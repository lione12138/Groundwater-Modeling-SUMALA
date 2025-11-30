#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MODFLOW6 Transient Model Calibration - Latin Hypercube Sampling (LHS) Optimization
===================================================================================
Automated calibration using LHS to sample parameter space for Specific Yield (Sy) 
and Specific Storage (Ss) across different lithologies.

Purpose:
    - Use Latin Hypercube Sampling (LHS) to efficiently explore parameter space
    - Sample in log-space for Ss (spans multiple orders of magnitude)
    - Run multiple MODFLOW6 simulations with different parameter combinations
    - Find optimal Sy and Ss values that minimize RMSE and maximize R²
    - Save best parameters and all results for analysis

Dependencies: numpy, pandas, flopy, scipy, fiona, shapely
"""

from __future__ import annotations
import os
import sys
import shutil
import glob
import time
import json
from pathlib import Path
from typing import Dict, List, Tuple
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.stats import qmc  # Latin Hypercube Sampling

import flopy
from flopy.mf6 import MFSimulation

# Optional: For shapefile processing
try:
    import fiona
    from shapely.geometry import shape, Point
    from flopy.utils import GridIntersect
    HAS_SPATIAL_LIBS = True
except ImportError:
    HAS_SPATIAL_LIBS = False
    print("Warning: fiona and/or shapely not available. Lithology mapping will be limited.")

# ------------------------
# User Configuration
# ------------------------
ROOT = Path(__file__).resolve().parent
CONFIG = {
    "base_model_ws": str((ROOT.parent / "Model").resolve()),
    "mf6_exe": r"D:\Program Files\USGS\mf6.6.2_win64\bin\mf6.exe",
    "runs_root": str((ROOT.parent / "Output" / "Calibrate").resolve()),
    "outputs_dir": str((ROOT.parent / "Output").resolve()),
    "obs_csv_name_hint": "Transient_cali.ob_gw_out_head",
    "observed_csv_path": str((ROOT.parent.parent / "Instructions" / "calibration_transient.csv").resolve()),
    "obs_points": ["P1", "Pz12"],
    # Mapping from MODFLOW observation names to user-friendly names
    "obs_name_mapping": {
        "HD__1_63_61": "P1",   # Layer 1, Row 63, Col 61
        "HD__2_63_61": "Pz12",  # Layer 2, Row 63, Col 61
    },
    "lithology_shapes": {
        "path": str((ROOT.parent.parent / "01_Steady_Calibrate" / "Lithology" / "KX-new.shp").resolve()),
        "name_field": "OBJECTNAME",
        "value_to_key": {
            "Clean_Gravel": "clean_gravel",
            "Clayey_Sand": "clayey_sand",
            "Loamy_Sand": "loamy_sand",
            "Sand_Stone": "sandstone"
        },
        "layer_by_key": {
            "clean_gravel": 0,
            "clayey_sand": 0,
            "loamy_sand": 0,
            "sandstone": 1
        }
    },
}

# ====================
# 👇 LHS Calibration Settings
# ====================
LHS_CONFIG = {
    "n_samples": 100,  # Number of parameter combinations to test (recommended: 50-200)
    "random_seed": 999,  # For reproducibility
    
    # Parameter bounds: [min, max]
    # Sy bounds (dimensionless): 0.01 - 0.35
    "sy_bounds": {
        "sy_clean_gravel": [0.1, 0.4],    # High permeability → high Sy
        "sy_clayey_sand": [0, 0.2],     # Clay reduces Sy
        "sy_loamy_sand": [0.1, 0.25],      # Medium Sy
        "sy_sandy_gravel": [0, 0.4],    # High Sy
        "sy_sandstone": [0.001, 0.001],     # Confined layer, Sy irrelevant (fixed)
        "sy_limestone": [0.001, 0.001],     # Confined layer, Sy irrelevant (fixed)
    },
    
    # Ss bounds (1/m): 1e-6 - 1e-3 (LOG-SCALE sampling)
    # Will be sampled in log10 space: log10(1e-6) to log10(1e-3) = -6 to -3
    "ss_log_bounds": {
        "ss_clean_gravel": [-7, -3],    # 10^-5.5 to 10^-4.5 (3.16e-6 to 3.16e-5)
        "ss_clayey_sand": [-7, -3],     # Higher due to clay compressibility
        "ss_loamy_sand": [-7, -3],      # Medium range
        "ss_sandy_gravel": [-7, -3],    # Similar to clean gravel
        "ss_sandstone": [-5, -4],       # Layer 2: critical parameter (1e-5 to 1e-4)
        "ss_limestone": [-7, -3],       # Lower compressibility (1e-6 to 1e-5)
    },
    
    # Optimization target
    "target_metric": "Overall_R2",  # Options: "Overall_R2", "Overall_RMSE", "P1_R2", "Pz12_R2"
    "maximize": True,  # True for R², False for RMSE
}

LITHOLOGY_TARGETS = [
    "clean_gravel",
    "clayey_sand",
    "loamy_sand",
    "sandstone",
    "sandy_gravel",
    "limestone"
]

# ------------------------
# Helper Functions (reused from Transient_test.py)
# ------------------------

def load_cellids_by_lithology(sim: MFSimulation, shapefile_path: str,
                               name_field: str, value_to_key: dict,
                               layer_by_key: dict) -> Dict[str, List[Tuple[int, int, int]]]:
    """Load shapefile and map cells to lithologies."""
    if not HAS_SPATIAL_LIBS:
        return {}
    
    gwf = sim.get_model()
    modelgrid = gwf.modelgrid
    gi = GridIntersect(modelgrid, method='structured')
    
    cellids_by_lith = {key: [] for key in value_to_key.values()}
    
    with fiona.open(shapefile_path, 'r') as shp:
        for feature in shp:
            lith_name = feature['properties'].get(name_field)
            if not lith_name or lith_name not in value_to_key:
                continue
            
            lith_key = value_to_key[lith_name]
            layer_idx = layer_by_key.get(lith_key, 0)
            
            geom = shape(feature['geometry'])
            result = gi.intersect(geom)
            if result is not None and len(result) > 0:
                if hasattr(result, 'cellids'):
                    # Structured result with cellids attribute
                    for cellid in result.cellids:
                        if isinstance(cellid, tuple) and len(cellid) == 2:
                            row, col = cellid
                            cellids_by_lith[lith_key].append((layer_idx, row, col))
                elif isinstance(result, np.ndarray) and 'cellids' in result.dtype.names:
                    # Numpy recarray with cellids field
                    for record in result:
                        cellid = record['cellids']
                        if isinstance(cellid, tuple) and len(cellid) == 2:
                            row, col = cellid
                            cellids_by_lith[lith_key].append((layer_idx, row, col))
    
    return cellids_by_lith


def apply_storage_parameters(sim: MFSimulation, params: dict,
                              cellids_by_lith: Dict[str, List[Tuple[int, int, int]]]) -> None:
    """Apply Sy and Ss parameters to storage package."""
    gwf = sim.get_model()
    sto = gwf.sto
    dis = gwf.dis
    
    nlay, nrow, ncol = dis.nlay.array, dis.nrow.array, dis.ncol.array
    
    # CRITICAL: Create NEW arrays each time, don't reuse old ones
    sy_arrays = {}
    ss_arrays = {}
    
    for lay in range(nlay):
        sy_arrays[lay] = np.full((nrow, ncol), 0.2)
        ss_arrays[lay] = np.full((nrow, ncol), 1e-5)
    
    # Apply shapefile lithologies
    for lith_key in ["clean_gravel", "clayey_sand", "loamy_sand", "sandstone"]:
        if lith_key not in cellids_by_lith or not cellids_by_lith[lith_key]:
            continue
        
        sy_val = params.get(f"sy_{lith_key}", 0.2)
        ss_val = params.get(f"ss_{lith_key}", 1e-5)
        
        for lay, row, col in cellids_by_lith[lith_key]:
            if lay in sy_arrays:
                sy_arrays[lay][row, col] = sy_val
            if lay in ss_arrays:
                ss_arrays[lay][row, col] = ss_val
    
    # Apply defaults to uncovered cells
    sy_sandy_gravel = params.get("sy_sandy_gravel", 0.25)
    ss_sandy_gravel = params.get("ss_sandy_gravel", 1e-5)
    
    for row in range(nrow):
        for col in range(ncol):
            is_covered = any((0, row, col) in cellids_by_lith[key] 
                           for key in ["clean_gravel", "clayey_sand", "loamy_sand"])
            if not is_covered and 0 in sy_arrays:
                sy_arrays[0][row, col] = sy_sandy_gravel
            if not is_covered and 0 in ss_arrays:
                ss_arrays[0][row, col] = ss_sandy_gravel
    
    sy_limestone = params.get("sy_limestone", 0.001)
    ss_limestone = params.get("ss_limestone", 1e-5)
    
    for row in range(nrow):
        for col in range(ncol):
            is_covered = (1, row, col) in cellids_by_lith.get("sandstone", [])
            if not is_covered and 1 in sy_arrays:
                sy_arrays[1][row, col] = sy_limestone
            if not is_covered and 1 in ss_arrays:
                ss_arrays[1][row, col] = ss_limestone
    
    # Write back to model
    sto.sy.set_data(sy_arrays)
    sto.ss.set_data(ss_arrays)
    
    # CRITICAL: Write simulation files so changes take effect
    sim.write_simulation()


def run_model(sim: MFSimulation, mf6_exe: str, verbose: bool = False) -> bool:
    """Run MODFLOW6 simulation."""
    # Don't call write_simulation() here - already done in apply_storage_parameters
    success, _ = sim.run_simulation(silent=not verbose)
    return success


def extract_drawdown(sim: MFSimulation, obs_csv_hint: str) -> pd.DataFrame:
    """Extract simulation results and convert to drawdown."""
    ws = sim.simulation_data.mfpath.get_sim_path()
    pattern = os.path.join(ws, f"*{obs_csv_hint}*.csv")
    csv_files = glob.glob(pattern)
    
    if not csv_files:
        raise FileNotFoundError(f"No observation CSV found matching: {pattern}")
    
    obs_file = csv_files[0]
    sim_df = pd.read_csv(obs_file)
    
    if 'totim' in sim_df.columns:
        sim_df['time (s)'] = sim_df['totim']
    elif 'time' in sim_df.columns:
        sim_df['time (s)'] = sim_df['time']
    else:
        raise ValueError("No time column found")
    
    drawdown_df = pd.DataFrame()
    drawdown_df['time (s)'] = sim_df['time (s)']
    
    # Get observation name mapping
    obs_mapping = CONFIG.get("obs_name_mapping", {})
    
    # Convert MODFLOW column names to user-friendly names and calculate drawdown
    for mf_name, user_name in obs_mapping.items():
        if mf_name in sim_df.columns:
            initial_head = sim_df[mf_name].iloc[0]
            drawdown_df[user_name] = initial_head - sim_df[mf_name]
    
    return drawdown_df


def compare_to_observations(sim_df: pd.DataFrame, obs_csv_path: str) -> Dict[str, float]:
    """Compare simulation to observed data and calculate metrics."""
    # Read observation data - skip first column which has semicolon-separated duplicate data
    obs_df = pd.read_csv(obs_csv_path)
    
    # The CSV has duplicate columns due to mixed delimiters - keep only the comma-separated ones
    if obs_df.columns[0].count(';') > 0:  # First column has semicolons
        obs_df = obs_df.iloc[:, 1:]  # Drop first column
    
    if 'time (s)' not in obs_df.columns and obs_df.columns[0] != 'time (s)':
        obs_df.rename(columns={obs_df.columns[0]: 'time (s)'}, inplace=True)
    
    for col in obs_df.columns:
        obs_df[col] = pd.to_numeric(obs_df[col], errors='coerce')
    
    obs_df = obs_df.dropna()
    
    # Ensure time column has consistent dtype (float)
    obs_df['time (s)'] = obs_df['time (s)'].astype(float)
    sim_df['time (s)'] = sim_df['time (s)'].astype(float)
    
    metrics = {}
    
    for obs_point in CONFIG["obs_points"]:
        if obs_point not in obs_df.columns or obs_point not in sim_df.columns:
            continue
        
        # Nearest neighbor interpolation
        merged = pd.merge_asof(
            obs_df[['time (s)', obs_point]].rename(columns={obs_point: 'observed'}),
            sim_df[['time (s)', obs_point]].rename(columns={obs_point: 'simulated'}),
            on='time (s)',
            direction='nearest'
        ).dropna()
        
        if len(merged) < 2:
            continue
        
        observed = merged['observed'].values
        simulated = merged['simulated'].values
        
        rmse = np.sqrt(np.mean((observed - simulated)**2))
        ss_res = np.sum((observed - simulated)**2)
        ss_tot = np.sum((observed - np.mean(observed))**2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        metrics[f"{obs_point}_RMSE"] = rmse
        metrics[f"{obs_point}_R2"] = r2
    
    # Overall metrics
    rmse_values = [v for k, v in metrics.items() if 'RMSE' in k]
    r2_values = [v for k, v in metrics.items() if 'R2' in k]
    
    if rmse_values:
        metrics['Overall_RMSE'] = np.mean(rmse_values)
        metrics['Overall_R2'] = np.mean(r2_values)
    
    return metrics


# ------------------------
# LHS Sampling Functions
# ------------------------

def generate_lhs_samples(n_samples: int, param_bounds: dict, 
                         log_bounds: dict, seed: int = 42) -> pd.DataFrame:
    """
    Generate Latin Hypercube Samples for calibration parameters.
    
    Parameters
    ----------
    n_samples : int
        Number of samples to generate
    param_bounds : dict
        Linear-scale parameter bounds {param_name: [min, max]}
    log_bounds : dict
        Log-scale parameter bounds {param_name: [log_min, log_max]}
    seed : int
        Random seed for reproducibility
    
    Returns
    -------
    pd.DataFrame
        Sampled parameters with columns for each parameter
    """
    n_params = len(param_bounds) + len(log_bounds)
    
    # Create LHS sampler
    sampler = qmc.LatinHypercube(d=n_params, seed=seed)
    samples = sampler.random(n=n_samples)  # Returns values in [0, 1]
    
    # Convert to parameter space
    param_names = list(param_bounds.keys()) + list(log_bounds.keys())
    param_df = pd.DataFrame(samples, columns=param_names)
    
    # Scale linear parameters
    for i, (param, bounds) in enumerate(param_bounds.items()):
        min_val, max_val = bounds
        param_df[param] = min_val + param_df[param] * (max_val - min_val)
        # Round to 3 significant figures
        param_df[param] = param_df[param].apply(lambda x: float(f'{x:.3g}'))
    
    # Scale log parameters (convert from [0,1] to log space, then to linear)
    for i, (param, log_bounds_val) in enumerate(log_bounds.items(), start=len(param_bounds)):
        log_min, log_max = log_bounds_val
        # Map [0,1] to [log_min, log_max] then convert to linear space
        log_val = log_min + param_df[param] * (log_max - log_min)
        param_df[param] = 10 ** log_val
        # Round to 3 significant figures
        param_df[param] = param_df[param].apply(lambda x: float(f'{x:.3g}'))
    
    return param_df


def run_single_iteration(iteration: int, params: dict, base_ws: str, 
                         mf6_exe: str, runs_root: str, obs_path: str,
                         cellids_by_lith: dict) -> Dict:
    """
    Run a single MODFLOW6 iteration with given parameters.
    
    Returns dictionary with iteration results including metrics.
    """
    iter_ws = Path(runs_root) / f"run_{iteration+1:03d}"
    
    # Clean and copy model
    if iter_ws.exists():
        import stat
        def remove_readonly(func, path, _):
            try:
                os.chmod(path, stat.S_IWRITE)
                func(path)
            except Exception:
                pass
        try:
            shutil.rmtree(iter_ws, onerror=remove_readonly)
        except Exception:
            pass  # Ignore errors in iteration cleanup
        time.sleep(0.1)
    
    shutil.copytree(base_ws, iter_ws)
    
    try:
        # Load simulation
        sim = MFSimulation.load(
            sim_ws=str(iter_ws),
            exe_name=mf6_exe,
            verbosity_level=0
        )
        
        # Apply parameters
        apply_storage_parameters(sim, params, cellids_by_lith)
        
        # Run model
        success = run_model(sim, mf6_exe, verbose=False)
        
        if not success:
            return {
                "iteration": iteration,
                "success": False,
                "error": "Model run failed",
                **{k: np.nan for k in params.keys()},
                "Overall_RMSE": np.nan,
                "Overall_R2": np.nan
            }
        
        # Extract and compare results
        try:
            sim_drawdown = extract_drawdown(sim, CONFIG["obs_csv_name_hint"])
            if len(sim_drawdown.columns) <= 1:  # Only has time column
                return {
                    "iteration": iteration,
                    "success": False,
                    "error": "No observation points found in simulation output",
                    **{k: np.nan for k in params.keys()},
                    "Overall_RMSE": np.nan,
                    "Overall_R2": np.nan
                }
            
            metrics = compare_to_observations(sim_drawdown, obs_path)
            
            if not metrics:  # Empty metrics dictionary
                return {
                    "iteration": iteration,
                    "success": False,
                    "error": "No metrics calculated - check observation point mapping",
                    **{k: np.nan for k in params.keys()},
                    "Overall_RMSE": np.nan,
                    "Overall_R2": np.nan
                }
        except Exception as e:
            return {
                "iteration": iteration,
                "success": False,
                "error": f"Metrics calculation failed: {str(e)}",
                **{k: np.nan for k in params.keys()},
                "Overall_RMSE": np.nan,
                "Overall_R2": np.nan
            }
        
        # Combine results
        result = {
            "iteration": iteration,
            "success": True,
            **params,
            **metrics
        }
        
        # Cleanup: Keep only observation output file
        try:
            obs_output_file = None
            for file in iter_ws.glob("*.ob_gw_out_head.csv"):
                obs_output_file = file
                break
            
            if obs_output_file and obs_output_file.exists():
                # Delete all files except the observation output
                for item in iter_ws.iterdir():
                    if item.is_file() and item != obs_output_file:
                        try:
                            os.chmod(item, 0o777)
                            item.unlink()
                        except Exception:
                            pass
                    elif item.is_dir():
                        # Remove subdirectories
                        try:
                            shutil.rmtree(item, ignore_errors=True)
                        except Exception:
                            pass
        except Exception as e:
            print(f"  Warning: Cleanup failed for run_{iteration+1:03d}: {e}")
        
        return result
        
    except Exception as e:
        return {
            "iteration": iteration,
            "success": False,
            "error": str(e),
            **{k: np.nan for k in params.keys()},
            "Overall_RMSE": np.nan,
            "Overall_R2": np.nan
        }


def main():
    """Main LHS calibration execution."""
    print("\n" + "="*80)
    print("MODFLOW6 Transient Calibration - Latin Hypercube Sampling (LHS)")
    print("="*80 + "\n")
    
    start_time = time.time()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Display configuration
    print("=== LHS Configuration ===")
    print(f"  Number of samples: {LHS_CONFIG['n_samples']}")
    print(f"  Random seed: {LHS_CONFIG['random_seed']}")
    print(f"  Target metric: {LHS_CONFIG['target_metric']} ({'maximize' if LHS_CONFIG['maximize'] else 'minimize'})")
    print()
    
    print("=== Parameter Bounds ===")
    print("\n--- Specific Yield (Sy) [Linear Scale] ---")
    for param, bounds in LHS_CONFIG['sy_bounds'].items():
        if bounds[0] != bounds[1]:  # Skip fixed parameters
            print(f"  {param:20s}: [{bounds[0]:.3f}, {bounds[1]:.3f}]")
    
    print("\n--- Specific Storage (Ss) [Log Scale] ---")
    for param, log_bounds in LHS_CONFIG['ss_log_bounds'].items():
        min_val = 10 ** log_bounds[0]
        max_val = 10 ** log_bounds[1]
        print(f"  {param:20s}: [{min_val:.2e}, {max_val:.2e}]  (log: [{log_bounds[0]:.1f}, {log_bounds[1]:.1f}])")
    print()
    
    # Generate LHS samples
    print("=== Generating LHS Samples ===")
    samples_df = generate_lhs_samples(
        LHS_CONFIG['n_samples'],
        LHS_CONFIG['sy_bounds'],
        LHS_CONFIG['ss_log_bounds'],
        LHS_CONFIG['random_seed']
    )
    print(f"✓ Generated {len(samples_df)} parameter combinations\n")
    
    # Prepare run directory
    runs_root = Path(CONFIG["runs_root"])
    if runs_root.exists():
        import stat
        def remove_readonly(func, path, _):
            try:
                os.chmod(path, stat.S_IWRITE)
                func(path)
            except Exception:
                pass
        
        try:
            shutil.rmtree(runs_root, onerror=remove_readonly)
        except Exception as e:
            print(f"⚠️ Cannot delete existing folder: {e}")
            backup = str(runs_root) + "_old_" + str(int(time.time()))
            try:
                os.rename(runs_root, backup)
                print(f"✓ Renamed old folder to: {backup}")
            except Exception as e2:
                raise RuntimeError(f"❌ Cannot handle target folder: {e2}") from e2
        time.sleep(0.2)
    
    runs_root.mkdir(parents=True, exist_ok=True)
    
    # Load lithology mapping once (reuse for all iterations)
    print("=== Loading Lithology Mapping ===")
    temp_ws = runs_root / "temp_mapping"
    shutil.copytree(CONFIG["base_model_ws"], temp_ws)
    
    sim_temp = MFSimulation.load(
        sim_ws=str(temp_ws),
        exe_name=CONFIG["mf6_exe"],
        verbosity_level=0
    )
    
    cellids_by_lith = {}
    if HAS_SPATIAL_LIBS and Path(CONFIG["lithology_shapes"]["path"]).exists():
        cellids_by_lith = load_cellids_by_lithology(
            sim_temp,
            CONFIG["lithology_shapes"]["path"],
            CONFIG["lithology_shapes"]["name_field"],
            CONFIG["lithology_shapes"]["value_to_key"],
            CONFIG["lithology_shapes"]["layer_by_key"]
        )
    
    # Clean temp directory
    import stat
    def remove_readonly(func, path, _):
        try:
            os.chmod(path, stat.S_IWRITE)
            func(path)
        except Exception:
            pass
    
    try:
        shutil.rmtree(temp_ws, onerror=remove_readonly)
    except Exception as e:
        print(f"⚠️ Warning: Could not delete temp folder: {e}")
    print()
    
    # Run iterations
    print("=== Running LHS Iterations ===")
    results = []
    
    for i in range(len(samples_df)):
        params = samples_df.iloc[i].to_dict()
        
        print(f"[{i+1}/{len(samples_df)}] Running iteration {i}...", end=" ", flush=True)
        
        result = run_single_iteration(
            i, params,
            CONFIG["base_model_ws"],
            CONFIG["mf6_exe"],
            str(runs_root),
            CONFIG["observed_csv_path"],
            cellids_by_lith
        )
        
        results.append(result)
        
        if result['success']:
            metric_val = result.get(LHS_CONFIG['target_metric'], np.nan)
            print(f"✓ {LHS_CONFIG['target_metric']} = {metric_val:.4f}")
        else:
            print(f"✗ {result.get('error', 'Unknown error')}")
    
    # Compile results
    results_df = pd.DataFrame(results)
    
    # Save all results
    output_dir = Path(CONFIG["outputs_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results_file = output_dir / "calibration_results.csv"
    results_df.to_csv(results_file, index=False)
    print(f"\n✓ All results saved to: {results_file}")
    
    # Find best result
    successful_results = results_df[results_df['success'] == True].copy()
    
    if len(successful_results) == 0:
        print("\n✗ No successful runs. Check model configuration and parameters.")
        return
    
    # Check if target metric exists
    if LHS_CONFIG['target_metric'] not in successful_results.columns:
        print(f"\n✗ Target metric '{LHS_CONFIG['target_metric']}' not found in results.")
        print(f"Available columns: {list(successful_results.columns)}")
        return
    
    # Filter out NaN values in target metric
    valid_results = successful_results[successful_results[LHS_CONFIG['target_metric']].notna()].copy()
    
    if len(valid_results) == 0:
        print(f"\n✗ No valid (non-NaN) values for metric '{LHS_CONFIG['target_metric']}'.")
        return
    
    if LHS_CONFIG['maximize']:
        best_idx = valid_results[LHS_CONFIG['target_metric']].idxmax()
    else:
        best_idx = valid_results[LHS_CONFIG['target_metric']].idxmin()
    
    best_result = valid_results.loc[best_idx]
    
    best_result = successful_results.loc[best_idx]
    
    # Display best result
    print("\n" + "="*80)
    print("=== BEST RESULT ===")
    print("="*80)
    print(f"\nIteration: {int(best_result['iteration'])}")
    print(f"Target Metric ({LHS_CONFIG['target_metric']}): {best_result[LHS_CONFIG['target_metric']]:.4f}")
    
    if 'Overall_RMSE' in best_result:
        print(f"Overall RMSE: {best_result['Overall_RMSE']:.4f} m")
    if 'Overall_R2' in best_result:
        print(f"Overall R²: {best_result['Overall_R2']:.4f}")
    
    if 'P1_RMSE' in best_result and 'P1_R2' in best_result:
        print(f"P1:  RMSE = {best_result['P1_RMSE']:.4f} m, R² = {best_result['P1_R2']:.4f}")
    if 'Pz12_RMSE' in best_result and 'Pz12_R2' in best_result:
        print(f"Pz12: RMSE = {best_result['Pz12_RMSE']:.4f} m, R² = {best_result['Pz12_R2']:.4f}")
    
    print("\n--- Best Parameters ---")
    print("\nSpecific Yield (Sy):")
    for lith in LITHOLOGY_TARGETS:
        param = f"sy_{lith}"
        if param in best_result:
            print(f"  {lith:20s}: {best_result[param]:.4f}")
    
    print("\nSpecific Storage (Ss):")
    for lith in LITHOLOGY_TARGETS:
        param = f"ss_{lith}"
        if param in best_result:
            print(f"  {lith:20s}: {best_result[param]:.2e} (1/m)")
    
    # Save best parameters
    best_params = {lith: {} for lith in LITHOLOGY_TARGETS}
    for lith in LITHOLOGY_TARGETS:
        sy_key = f"sy_{lith}"
        ss_key = f"ss_{lith}"
        if sy_key in best_result:
            best_params[lith]['Sy'] = float(best_result[sy_key])
        if ss_key in best_result:
            best_params[lith]['Ss'] = float(best_result[ss_key])
    
    best_params_file = output_dir / "best_parameters.json"
    with open(best_params_file, 'w') as f:
        json.dump({
            "timestamp": timestamp,
            "target_metric": LHS_CONFIG['target_metric'],
            "metric_value": float(best_result[LHS_CONFIG['target_metric']]),
            "iteration": int(best_result['iteration']),
            "parameters": best_params
        }, f, indent=2)
    
    print(f"\n✓ Best parameters saved to: {best_params_file}")
    
    # Summary statistics
    print("\n" + "="*80)
    print("=== SUMMARY STATISTICS ===")
    print("="*80)
    
    if LHS_CONFIG['target_metric'] in successful_results.columns:
        metric_values = successful_results[LHS_CONFIG['target_metric']]
        print(f"\n{LHS_CONFIG['target_metric']}:")
        print(f"  Best:   {metric_values.max():.4f}")
        print(f"  Worst:  {metric_values.min():.4f}")
        print(f"  Mean:   {metric_values.mean():.4f}")
        print(f"  Median: {metric_values.median():.4f}")
        print(f"  Std:    {metric_values.std():.4f}")
    
    success_rate = len(successful_results) / len(results_df) * 100
    print(f"\nSuccess Rate: {success_rate:.1f}% ({len(successful_results)}/{len(results_df)})")
    
    elapsed = time.time() - start_time
    print(f"\nTotal Time: {elapsed/60:.1f} minutes ({elapsed:.0f} seconds)")
    print(f"Average Time per Iteration: {elapsed/len(results_df):.1f} seconds")
    
    print("\n" + "="*80)
    print("✓ LHS Calibration Completed Successfully!")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()
