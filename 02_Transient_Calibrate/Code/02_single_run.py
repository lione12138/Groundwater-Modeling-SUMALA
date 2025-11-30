#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MODFLOW6 Transient Model Calibration - Manual Parameter Testing
================================================================
Calibrate transient model by adjusting Specific Yield (Sy) and Specific Storage (Ss)
for different lithologies while keeping all boundary conditions unchanged.

Purpose:
    - Modify storage parameters (Sy and Ss) for 6 lithologies (4 from shapefile + 2 defaults)
    - Run MODFLOW6 transient simulation
    - Extract drawdown results for observation points (P1 and Pz12)
    - Calculate RMSE and R² to evaluate model fit

Dependencies: numpy, pandas, flopy, fiona, shapely
"""

from __future__ import annotations
import os
import sys
import shutil
import glob
import time
import warnings
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

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
    "test_run_ws": str((ROOT.parent / "Output" / "run_test").resolve()),
    "outputs_dir": str((ROOT.parent / "Output").resolve()),
    "obs_csv_name_hint": "Transient_cali.ob_gw_out_head",
    "observed_csv_path": str((ROOT.parent.parent / "Instructions" / "calibration_transient.csv").resolve()),
    "obs_points": ["P1", "Pz12"],  # Observation points to calibrate
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
            "clean_gravel": 0,    # Layer 1 (unconfined)
            "clayey_sand": 0,     # Layer 1
            "loamy_sand": 0,      # Layer 1
            "sandstone": 1        # Layer 2 (confined)
        }
    },
}

# ====================
# 👇 MANUAL PARAMETERS - Adjust these values for calibration
# ====================
MANUAL_PARAMS = {
    # ===== Specific Yield (Sy) - Only affects Layer 1 (unconfined, ICONVERT=1) =====
    # Dimensionless, typical range: 0.01–0.35
    # Represents drainable porosity for unconfined aquifer

    "sy_clean_gravel": 0.30,   # Clean gravel
    "sy_clayey_sand": 0.15,    # Clayey sand
    "sy_loamy_sand": 0.19,      # Loamy sand
    "sy_sandy_gravel": 0.22,   # Sandy gravel

    # Note: Layer 2 is confined (ICONVERT=0), Sy irrelevant but must be provided
    "sy_sandstone": 0.001,                    # Sandstone
    "sy_limestone": 0.001,                    # Limestone

    # ===== Specific Storage (Ss) - Affects both Layer 1 and Layer 2 =====
    # Unit: 1/m, typical range: 1e-6 to 1e-3

    "ss_clean_gravel": 3.37e-05,    # Clean gravel
    "ss_clayey_sand": 5.37e-03,     # Clayey sand
    "ss_loamy_sand": 9.7e-03,      # Loamy sand
    "ss_sandy_gravel":2.66e-05,    # Sandy gravel
    "ss_sandstone": 8.58e-04,       # Sandstone
    "ss_limestone": 9.89e-03,       # Limestone
}


# MANUAL_PARAMS = {
#     # ===== Specific Yield (Sy) - Only affects Layer 1 (unconfined, ICONVERT=1) =====
#     # Dimensionless, typical range: 0.01–0.35
#     # Represents drainable porosity for unconfined aquifer

#     "sy_clean_gravel": 0.842,   # Clean gravel
#     "sy_clayey_sand": 0.943,    # Clayey sand
#     "sy_loamy_sand": 0.675,      # Loamy sand
#     "sy_sandy_gravel": 0.156,   # Sandy gravel

#     # Note: Layer 2 is confined (ICONVERT=0), Sy irrelevant but must be provided
#     "sy_sandstone": 0.0031,                    # Sandstone
#     "sy_limestone": 0.0079,                    # Limestone

#     # ===== Specific Storage (Ss) - Affects both Layer 1 and Layer 2 =====
#     # Unit: 1/m, typical range: 1e-6 to 1e-3

#     "ss_clean_gravel": 7.77e-02,    # Clean gravel
#     "ss_clayey_sand": 7.69e-04,     # Clayey sand
#     "ss_loamy_sand": 2.61e-04,      # Loamy sand
#     "ss_sandy_gravel":3.42e-02,    # Sandy gravel
#     "ss_sandstone": 7.79e-05,       # Sandstone
#     "ss_limestone": 3.83e-04,       # Limestone
# }


LITHOLOGY_TARGETS = [
    "clean_gravel",
    "clayey_sand",
    "loamy_sand",
    "sandstone",
    "sandy_gravel",
    "limestone"
]

# ------------------------
# Helper Functions
# ------------------------

def load_cellids_by_lithology(sim: MFSimulation, shapefile_path: str,
                               name_field: str, value_to_key: dict,
                               layer_by_key: dict) -> Dict[str, List[Tuple[int, int, int]]]:
    """
    Load shapefile and identify which model cells belong to each lithology.
    
    Parameters
    ----------
    sim : MFSimulation
        The FloPy MFSimulation object
    shapefile_path : str
        Path to the lithology shapefile
    name_field : str
        Field name in shapefile containing lithology names
    value_to_key : dict
        Mapping from shapefile values to parameter keys
    layer_by_key : dict
        Mapping from parameter keys to layer indices
    
    Returns
    -------
    dict
        Dictionary mapping lithology keys to lists of (layer, row, col) tuples
    """
    if not HAS_SPATIAL_LIBS:
        print("WARNING: Cannot load shapefile without fiona and shapely")
        return {}
    
    gwf = sim.get_model()
    dis = gwf.dis
    modelgrid = gwf.modelgrid
    
    # Create grid intersect object for efficient spatial queries
    # Suppress deprecation warning for structured method
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=DeprecationWarning)
        gi = GridIntersect(modelgrid, method='structured')
    
    cellids_by_lith = {key: [] for key in value_to_key.values()}
    
    # Read shapefile features
    with fiona.open(shapefile_path, 'r') as shp:
        for feature in shp:
            lith_name = feature['properties'].get(name_field)
            if not lith_name or lith_name not in value_to_key:
                continue
            
            lith_key = value_to_key[lith_name]
            layer_idx = layer_by_key.get(lith_key, 0)
            
            geom = shape(feature['geometry'])
            
            # Find all cells intersecting this feature
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
    
    # Print summary
    print("\n=== Lithology Cell Mapping Summary ===")
    for key, cellids in cellids_by_lith.items():
        if cellids:
            layer = cellids[0][0]
            print(f"  {key:20s}: {len(cellids):5d} cells in Layer {layer+1}")
    print()
    
    return cellids_by_lith


def apply_storage_parameters(sim: MFSimulation, params: dict,
                              cellids_by_lith: Dict[str, List[Tuple[int, int, int]]]) -> None:
    """
    Apply Sy and Ss parameters to the storage package based on lithology mapping.
    
    Parameters
    ----------
    sim : MFSimulation
        The FloPy MFSimulation object
    params : dict
        Dictionary of parameters (sy_* and ss_* keys)
    cellids_by_lith : dict
        Mapping from lithology keys to cell lists
    """
    gwf = sim.get_model()
    sto = gwf.sto
    dis = gwf.dis
    
    nlay, nrow, ncol = dis.nlay.array, dis.nrow.array, dis.ncol.array
    
    # CRITICAL: Create NEW arrays each time, don't reuse old ones from get_data()
    # This ensures parameters actually change between runs
    sy_arrays = {}
    ss_arrays = {}
    
    for lay in range(nlay):
        sy_arrays[lay] = np.full((nrow, ncol), 0.2)
        ss_arrays[lay] = np.full((nrow, ncol), 1e-5)
    
    # Apply shapefile-based lithologies
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
    
    # Apply defaults to cells not covered by shapefile
    # Layer 1: sandy_gravel (default for uncovered cells)
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
    
    # Layer 2: limestone (default for uncovered cells)
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
    
    print("=== Storage Parameters Applied ===")
    print(f"  Layer 1 Sy range: {np.min(sy_arrays[0]):.4f} - {np.max(sy_arrays[0]):.4f}")
    print(f"  Layer 1 Ss range: {np.min(ss_arrays[0]):.4e} - {np.max(ss_arrays[0]):.4e}")
    print(f"  Layer 2 Sy range: {np.min(sy_arrays[1]):.4f} - {np.max(sy_arrays[1]):.4f}")
    print(f"  Layer 2 Ss range: {np.min(ss_arrays[1]):.4e} - {np.max(ss_arrays[1]):.4e}")
    print()


def run_model(sim: MFSimulation, mf6_exe: str, verbose: bool = True) -> bool:
    """
    Run the MODFLOW6 simulation.
    
    Parameters
    ----------
    sim : MFSimulation
        The simulation object
    mf6_exe : str
        Path to mf6 executable
    verbose : bool
        Whether to print output
    
    Returns
    -------
    bool
        True if run was successful
    """
    print("=== Running MODFLOW6 ===")
    sim.write_simulation()
    success, _ = sim.run_simulation(silent=not verbose)
    
    if success:
        print("✓ MODFLOW6 run completed successfully\n")
    else:
        print("✗ MODFLOW6 run failed\n")
    
    return success


def extract_drawdown(sim: MFSimulation, obs_csv_hint: str) -> pd.DataFrame:
    """
    Extract simulation results and convert to drawdown.
    
    Drawdown = initial_head (t=1s) - current_head
    
    Parameters
    ----------
    sim : MFSimulation
        The simulation object
    obs_csv_hint : str
        Hint for finding observation CSV file
    
    Returns
    -------
    pd.DataFrame
        Drawdown data with columns: time (s), P1, Pz12
    """
    ws = sim.simulation_data.mfpath.get_sim_path()
    pattern = os.path.join(ws, f"*{obs_csv_hint}*.csv")
    csv_files = glob.glob(pattern)
    
    if not csv_files:
        raise FileNotFoundError(f"No observation CSV found matching: {pattern}")
    
    obs_file = csv_files[0]
    print(f"=== Reading simulation results ===")
    print(f"  File: {Path(obs_file).name}")
    
    # Read simulation output
    sim_df = pd.read_csv(obs_file)
    
    # Convert totim (simulation time in seconds) to time column
    if 'totim' in sim_df.columns:
        sim_df['time (s)'] = sim_df['totim']
    elif 'time' in sim_df.columns:
        sim_df['time (s)'] = sim_df['time']
    else:
        raise ValueError("No time column found in simulation output")
    
    # Calculate drawdown: initial_head - current_head
    # Initial head = head at first timestep (t=1s)
    drawdown_df = pd.DataFrame()
    drawdown_df['time (s)'] = sim_df['time (s)']
    
    # Get observation name mapping
    obs_mapping = CONFIG.get("obs_name_mapping", {})
    
    # Convert MODFLOW column names to user-friendly names and calculate drawdown
    for mf_name, user_name in obs_mapping.items():
        if mf_name in sim_df.columns:
            initial_head = sim_df[mf_name].iloc[0]  # Head at t=1s
            drawdown_df[user_name] = initial_head - sim_df[mf_name]
            print(f"  {user_name}: Initial head = {initial_head:.3f} m")
    
    print(f"  Extracted {len(drawdown_df)} timesteps\n")
    
    return drawdown_df


def compare_to_observations(sim_df: pd.DataFrame, obs_csv_path: str) -> Dict[str, float]:
    """
    Compare simulation results to observed drawdown and calculate metrics.
    
    Parameters
    ----------
    sim_df : pd.DataFrame
        Simulation drawdown data
    obs_csv_path : str
        Path to observed data CSV
    
    Returns
    -------
    dict
        Dictionary of metrics (RMSE, R²) for each observation point
    """
    # Read observed data
    obs_df = pd.read_csv(obs_csv_path)
    
    # The CSV has duplicate columns due to mixed delimiters - keep only the comma-separated ones
    if obs_df.columns[0].count(';') > 0:  # First column has semicolons
        obs_df = obs_df.iloc[:, 1:]  # Drop first column
    
    # Ensure time column
    if 'time (s)' not in obs_df.columns and obs_df.columns[0] != 'time (s)':
        obs_df.rename(columns={obs_df.columns[0]: 'time (s)'}, inplace=True)
    
    # Convert to numeric
    for col in obs_df.columns:
        obs_df[col] = pd.to_numeric(obs_df[col], errors='coerce')
    
    # Drop rows where ALL values are NA, but keep rows with some valid data
    obs_df = obs_df.dropna(how='all')
    
    # Ensure time column has consistent dtype (float)
    obs_df['time (s)'] = obs_df['time (s)'].astype(float)
    sim_df['time (s)'] = sim_df['time (s)'].astype(float)
    
    print("=== Comparison with Observations ===")
    metrics = {}
    
    for obs_point in CONFIG["obs_points"]:
        if obs_point not in obs_df.columns or obs_point not in sim_df.columns:
            print(f"  {obs_point}: SKIPPED (missing in data)")
            continue
        
        # Use nearest neighbor interpolation (as determined optimal)
        merged = pd.merge_asof(
            obs_df[['time (s)', obs_point]].rename(columns={obs_point: 'observed'}),
            sim_df[['time (s)', obs_point]].rename(columns={obs_point: 'simulated'}),
            on='time (s)',
            direction='nearest'
        ).dropna()
        
        if len(merged) < 2:
            print(f"  {obs_point}: SKIPPED (insufficient data)")
            continue
        
        observed = merged['observed'].values
        simulated = merged['simulated'].values
        
        # Calculate RMSE
        rmse = np.sqrt(np.mean((observed - simulated)**2))
        
        # Calculate R²
        ss_res = np.sum((observed - simulated)**2)
        ss_tot = np.sum((observed - np.mean(observed))**2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        metrics[f"{obs_point}_RMSE"] = rmse
        metrics[f"{obs_point}_R2"] = r2
        
        print(f"  {obs_point}: RMSE = {rmse:.4f} m, R² = {r2:.4f}")
    
    # Overall metrics
    rmse_values = [v for k, v in metrics.items() if 'RMSE' in k]
    r2_values = [v for k, v in metrics.items() if 'R2' in k]
    
    if rmse_values:
        metrics['Overall_RMSE'] = np.mean(rmse_values)
        metrics['Overall_R2'] = np.mean(r2_values)
        print(f"\n  Overall: RMSE = {metrics['Overall_RMSE']:.4f} m, R² = {metrics['Overall_R2']:.4f}")
    
    print()
    return metrics


def main():
    """Main execution function."""
    print("\n" + "="*70)
    print("MODFLOW6 Transient Model - Manual Parameter Testing")
    print("="*70 + "\n")
    
    # Display parameters
    print("=== Current Parameters ===")
    print("\n--- Specific Yield (Sy) ---")
    for key in LITHOLOGY_TARGETS:
        sy_key = f"sy_{key}"
        if sy_key in MANUAL_PARAMS:
            print(f"  {key:20s}: {MANUAL_PARAMS[sy_key]:.4f}")
    
    print("\n--- Specific Storage (Ss) ---")
    for key in LITHOLOGY_TARGETS:
        ss_key = f"ss_{key}"
        if ss_key in MANUAL_PARAMS:
            print(f"  {key:20s}: {MANUAL_PARAMS[ss_key]:.2e} (1/m)")
    print()
    
    # Prepare run directory
    test_ws = Path(CONFIG["test_run_ws"])
    if test_ws.exists():
        import stat
        def remove_readonly(func, path, _):
            """Clear readonly attribute and retry deletion"""
            try:
                os.chmod(path, stat.S_IWRITE)
                func(path)
            except Exception:
                pass
        
        try:
            shutil.rmtree(test_ws, onerror=remove_readonly)
        except Exception as e:
            print(f"⚠️ Cannot delete existing folder: {e}")
            # Try renaming instead
            backup = str(test_ws) + "_old_" + str(int(time.time()))
            try:
                os.rename(test_ws, backup)
                print(f"✓ Renamed old folder to: {backup}")
            except Exception as e2:
                raise RuntimeError(
                    f"❌ Cannot handle target folder: {e2}\n"
                    f"Suggestions:\n"
                    f"1. Close any programs using the folder\n"
                    f"2. Pause OneDrive sync temporarily\n"
                    f"3. Manually delete: {test_ws}"
                ) from e2
        
        time.sleep(0.2)  # Wait for filesystem
    
    shutil.copytree(CONFIG["base_model_ws"], test_ws)
    print(f"✓ Copied base model to: {test_ws}\n")
    
    # Load simulation
    print("=== Loading MODFLOW6 Model ===")
    sim = MFSimulation.load(
        sim_ws=str(test_ws),
        exe_name=CONFIG["mf6_exe"],
        verbosity_level=0
    )
    print(f"✓ Loaded simulation from: {test_ws}\n")
    
    # Load lithology mapping
    cellids_by_lith = {}
    if HAS_SPATIAL_LIBS and Path(CONFIG["lithology_shapes"]["path"]).exists():
        cellids_by_lith = load_cellids_by_lithology(
            sim,
            CONFIG["lithology_shapes"]["path"],
            CONFIG["lithology_shapes"]["name_field"],
            CONFIG["lithology_shapes"]["value_to_key"],
            CONFIG["lithology_shapes"]["layer_by_key"]
        )
    
    # Apply storage parameters
    apply_storage_parameters(sim, MANUAL_PARAMS, cellids_by_lith)
    
    # Run model
    success = run_model(sim, CONFIG["mf6_exe"], verbose=False)
    
    if not success:
        print("✗ Model run failed. Exiting.")
        return
    
    # Extract results
    try:
        sim_drawdown = extract_drawdown(sim, CONFIG["obs_csv_name_hint"])
        
        # Compare to observations
        metrics = compare_to_observations(sim_drawdown, CONFIG["observed_csv_path"])
        
        # Save results
        output_file = test_ws / "test_results.txt"
        with open(output_file, 'w') as f:
            f.write("MODFLOW6 Transient Model - Test Run Results\n")
            f.write("="*50 + "\n\n")
            f.write("Parameters:\n")
            for key in LITHOLOGY_TARGETS:
                sy_val = MANUAL_PARAMS.get(f"sy_{key}", "N/A")
                ss_val = MANUAL_PARAMS.get(f"ss_{key}", "N/A")
                f.write(f"  {key}:\n")
                f.write(f"    Sy = {sy_val}\n")
                f.write(f"    Ss = {ss_val}\n")
            f.write("\nMetrics:\n")
            for key, val in metrics.items():
                f.write(f"  {key} = {val:.4f}\n")
        
        print(f"✓ Results saved to: {output_file}")
        
        print("\n" + "="*70)
        print("✓ Test run completed successfully!")
        print("="*70 + "\n")
        
    except Exception as e:
        print(f"✗ Error during result extraction: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
