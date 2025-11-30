#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Transient Model Sensitivity Analysis
====================================
Performs One-at-a-Time (OAT) sensitivity analysis on storage parameters (Sy and Ss).
Each parameter is perturbed by ±10% and ±20% while others remain at base values.

Output:
- Sensitivity coefficients for each parameter
- Tornado diagrams
- Spider plots
- Summary tables
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
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('TkAgg')
warnings.filterwarnings('ignore')

import flopy
from flopy.mf6 import MFSimulation

# Optional: For shapefile processing
try:
    import fiona
    from shapely.geometry import shape
    from flopy.utils import GridIntersect
    HAS_SPATIAL_LIBS = True
except ImportError:
    HAS_SPATIAL_LIBS = False

# ============================================================================
# Configuration
# ============================================================================
ROOT = Path(__file__).resolve().parent
CONFIG = {
    "base_model_ws": str((ROOT.parent / "Model").resolve()),
    "mf6_exe": r"D:\Program Files\USGS\mf6.6.2_win64\bin\mf6.exe",
    "sensitivity_ws": str((ROOT.parent / "Output" / "Sensitivity Ana").resolve()),
    "outputs_dir": str((ROOT.parent / "Output" / "Sensitivity Ana").resolve()),
    "obs_csv_name_hint": "Transient_cali.ob_gw_out_head",
    "observed_csv_path": str((ROOT.parent.parent / "Instructions" / "calibration_transient.csv").resolve()),
    "obs_points": ["P1", "Pz12"],
    "obs_name_mapping": {
        "HD__1_63_61": "P1",
        "HD__2_63_61": "Pz12",
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

# ============================================================================
# Base (Calibrated) Parameters
# ============================================================================
BASE_PARAMS = {
    # Specific Yield (Sy) - Layer 1 (unconfined)
    "sy_clean_gravel": 0.30,
    "sy_clayey_sand": 0.15,
    "sy_loamy_sand": 0.19,
    "sy_sandy_gravel": 0.22,
    "sy_sandstone": 0.001,
    "sy_limestone": 0.001,

    # Specific Storage (Ss) - Both layers
    "ss_clean_gravel": 3.37e-05,
    "ss_clayey_sand": 5.37e-03,
    "ss_loamy_sand": 9.7e-03,
    "ss_sandy_gravel": 2.66e-05,
    "ss_sandstone": 8.58e-04,
    "ss_limestone": 9.89e-03,
}

# Parameters to analyze (12 storage parameters)
PARAMS_TO_ANALYZE = [
    # Sy parameters (Layer 1 relevant)
    "sy_clean_gravel",
    "sy_clayey_sand",
    "sy_loamy_sand",
    "sy_sandy_gravel",
    # Ss parameters (both layers)
    "ss_clean_gravel",
    "ss_clayey_sand",
    "ss_loamy_sand",
    "ss_sandy_gravel",
    "ss_sandstone",
    "ss_limestone",
]

# Perturbation levels
PERTURBATIONS = [-0.20, -0.10, 0.10, 0.20]  # ±10%, ±20%

# Lithology targets
LITHOLOGY_TARGETS = [
    "clean_gravel",
    "clayey_sand",
    "loamy_sand",
    "sandstone",
    "sandy_gravel",
    "limestone"
]

# Display names for plots
PARAM_DISPLAY_NAMES = {
    "sy_clean_gravel": "Sy (Clean Gravel)",
    "sy_clayey_sand": "Sy (Clayey Sand)",
    "sy_loamy_sand": "Sy (Loamy Sand)",
    "sy_sandy_gravel": "Sy (Sandy Gravel)",
    "sy_sandstone": "Sy (Sandstone)",
    "sy_limestone": "Sy (Limestone)",
    "ss_clean_gravel": "Ss (Clean Gravel)",
    "ss_clayey_sand": "Ss (Clayey Sand)",
    "ss_loamy_sand": "Ss (Loamy Sand)",
    "ss_sandy_gravel": "Ss (Sandy Gravel)",
    "ss_sandstone": "Ss (Sandstone)",
    "ss_limestone": "Ss (Limestone)",
}

# ============================================================================
# Model Functions
# ============================================================================

def load_cellids_by_lithology(sim: MFSimulation) -> Dict[str, List[Tuple[int, int, int]]]:
    """Load cell IDs for each lithology from shapefile."""
    if not HAS_SPATIAL_LIBS:
        return {}
    
    gwf = sim.get_model()
    modelgrid = gwf.modelgrid
    
    shp_cfg = CONFIG.get("lithology_shapes", {})
    shp_path = shp_cfg.get("path")
    name_field = shp_cfg.get("name_field", "OBJECTNAME")
    value_to_key = shp_cfg.get("value_to_key", {})
    layer_by_key = shp_cfg.get("layer_by_key", {})
    
    if not os.path.exists(shp_path):
        return {}
    
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=DeprecationWarning)
        gi = GridIntersect(modelgrid, method='structured')
    
    cellids_by_lith = {key: [] for key in value_to_key.values()}
    
    with fiona.open(shp_path, 'r') as shp:
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
                    for cellid in result.cellids:
                        if isinstance(cellid, tuple) and len(cellid) == 2:
                            row, col = cellid
                            cellids_by_lith[lith_key].append((layer_idx, row, col))
                elif isinstance(result, np.ndarray) and 'cellids' in result.dtype.names:
                    for record in result:
                        cellid = record['cellids']
                        if isinstance(cellid, tuple) and len(cellid) == 2:
                            row, col = cellid
                            cellids_by_lith[lith_key].append((layer_idx, row, col))
    
    return cellids_by_lith


def copy_model_to_run(base_ws: str, run_ws: str):
    """Copy model to run directory."""
    import stat
    
    if os.path.exists(run_ws):
        def remove_readonly(func, path, _):
            try:
                os.chmod(path, stat.S_IWRITE)
                func(path)
            except Exception:
                pass
        try:
            shutil.rmtree(run_ws, onerror=remove_readonly)
        except Exception:
            backup = run_ws + "_old_" + str(int(time.time()))
            os.rename(run_ws, backup)
    
    time.sleep(0.2)
    shutil.copytree(base_ws, run_ws, dirs_exist_ok=True)


def apply_storage_parameters(sim: MFSimulation, params: dict,
                              cellids_by_lith: Dict[str, List[Tuple[int, int, int]]]) -> None:
    """Apply Sy and Ss parameters to the storage package."""
    gwf = sim.get_model()
    sto = gwf.sto
    dis = gwf.dis
    
    nlay, nrow, ncol = dis.nlay.array, dis.nrow.array, dis.ncol.array
    
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
    
    # Apply defaults for uncovered cells
    sy_sandy_gravel = params.get("sy_sandy_gravel", 0.25)
    ss_sandy_gravel = params.get("ss_sandy_gravel", 1e-5)
    
    for row in range(nrow):
        for col in range(ncol):
            is_covered = any((0, row, col) in cellids_by_lith.get(key, []) 
                           for key in ["clean_gravel", "clayey_sand", "loamy_sand"])
            if not is_covered and 0 in sy_arrays:
                sy_arrays[0][row, col] = sy_sandy_gravel
            if not is_covered and 0 in ss_arrays:
                ss_arrays[0][row, col] = ss_sandy_gravel
    
    # Layer 2: limestone default
    sy_limestone = params.get("sy_limestone", 0.001)
    ss_limestone = params.get("ss_limestone", 1e-5)
    
    for row in range(nrow):
        for col in range(ncol):
            is_covered = (1, row, col) in cellids_by_lith.get("sandstone", [])
            if not is_covered and 1 in sy_arrays:
                sy_arrays[1][row, col] = sy_limestone
            if not is_covered and 1 in ss_arrays:
                ss_arrays[1][row, col] = ss_limestone
    
    sto.sy.set_data(sy_arrays)
    sto.ss.set_data(ss_arrays)
    sim.write_simulation()


def run_simulation(sim: MFSimulation) -> bool:
    """Run MODFLOW 6 simulation."""
    success, _ = sim.run_simulation(report=False, silent=True)
    return bool(success)


def extract_drawdown(sim: MFSimulation) -> pd.DataFrame:
    """Extract simulation results and convert to drawdown."""
    ws = sim.simulation_data.mfpath.get_sim_path()
    pattern = os.path.join(ws, f"*{CONFIG['obs_csv_name_hint']}*.csv")
    csv_files = glob.glob(pattern)
    
    if not csv_files:
        return None
    
    obs_file = csv_files[0]
    sim_df = pd.read_csv(obs_file)
    
    if 'totim' in sim_df.columns:
        sim_df['time (s)'] = sim_df['totim']
    elif 'time' in sim_df.columns:
        sim_df['time (s)'] = sim_df['time']
    else:
        return None
    
    drawdown_df = pd.DataFrame()
    drawdown_df['time (s)'] = sim_df['time (s)']
    
    obs_mapping = CONFIG.get("obs_name_mapping", {})
    
    for mf_name, user_name in obs_mapping.items():
        if mf_name in sim_df.columns:
            initial_head = sim_df[mf_name].iloc[0]
            drawdown_df[user_name] = initial_head - sim_df[mf_name]
    
    return drawdown_df


def calculate_metrics(sim_df: pd.DataFrame) -> Dict[str, Tuple[float, float]]:
    """Calculate RMSE and R² for each observation point."""
    obs_csv_path = CONFIG["observed_csv_path"]
    
    obs_df = pd.read_csv(obs_csv_path)
    
    if obs_df.columns[0].count(';') > 0:
        obs_df = obs_df.iloc[:, 1:]
    
    if 'time (s)' not in obs_df.columns and obs_df.columns[0] != 'time (s)':
        obs_df.rename(columns={obs_df.columns[0]: 'time (s)'}, inplace=True)
    
    for col in obs_df.columns:
        obs_df[col] = pd.to_numeric(obs_df[col], errors='coerce')
    
    obs_df = obs_df.dropna(how='all')
    obs_df['time (s)'] = obs_df['time (s)'].astype(float)
    sim_df['time (s)'] = sim_df['time (s)'].astype(float)
    
    metrics = {}
    
    for obs_point in CONFIG["obs_points"]:
        if obs_point not in obs_df.columns or obs_point not in sim_df.columns:
            continue
        
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
        
        metrics[obs_point] = (rmse, r2)
    
    # Overall metrics
    if metrics:
        rmse_vals = [v[0] for v in metrics.values()]
        r2_vals = [v[1] for v in metrics.values()]
        metrics['Overall'] = (np.mean(rmse_vals), np.mean(r2_vals))
    
    return metrics


# ============================================================================
# Sensitivity Analysis Functions
# ============================================================================

def run_single_model(params: Dict[str, float], run_name: str, cellids: Dict) -> Dict[str, Tuple[float, float]]:
    """Run a single model with given parameters and return metrics."""
    base_ws = CONFIG["base_model_ws"]
    sensitivity_ws = CONFIG["sensitivity_ws"]
    mf6_exe = CONFIG["mf6_exe"]
    
    run_ws = os.path.join(sensitivity_ws, "runs", run_name)
    
    copy_model_to_run(base_ws, run_ws)
    
    sim = MFSimulation.load(sim_ws=run_ws, exe_name=mf6_exe, verbosity_level=0)
    apply_storage_parameters(sim, params, cellids)
    
    success = run_simulation(sim)
    
    if not success:
        return {"P1": (float("inf"), float("nan")), "Pz12": (float("inf"), float("nan")), "Overall": (float("inf"), float("nan"))}
    
    sim_df = extract_drawdown(sim)
    if sim_df is None:
        return {"P1": (float("inf"), float("nan")), "Pz12": (float("inf"), float("nan")), "Overall": (float("inf"), float("nan"))}
    
    return calculate_metrics(sim_df)


def perform_sensitivity_analysis():
    """Perform OAT sensitivity analysis."""
    print("\n" + "="*70)
    print("TRANSIENT MODEL SENSITIVITY ANALYSIS")
    print("="*70)
    
    output_dir = Path(CONFIG["outputs_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "runs").mkdir(exist_ok=True)
    
    # Load cellids once
    print("\nLoading lithology cell mappings...")
    base_ws = CONFIG["base_model_ws"]
    mf6_exe = CONFIG["mf6_exe"]
    
    temp_ws = str(output_dir / "runs" / "temp_load")
    copy_model_to_run(base_ws, temp_ws)
    sim = MFSimulation.load(sim_ws=temp_ws, exe_name=mf6_exe, verbosity_level=0)
    cellids = load_cellids_by_lithology(sim)
    
    # Run base case
    print("\n" + "-"*70)
    print("Running BASE CASE...")
    print("-"*70)
    base_metrics = run_single_model(BASE_PARAMS.copy(), "base", cellids)
    
    base_rmse_p1 = base_metrics.get("P1", (np.nan, np.nan))[0]
    base_r2_p1 = base_metrics.get("P1", (np.nan, np.nan))[1]
    base_rmse_pz12 = base_metrics.get("Pz12", (np.nan, np.nan))[0]
    base_r2_pz12 = base_metrics.get("Pz12", (np.nan, np.nan))[1]
    base_rmse_overall = base_metrics.get("Overall", (np.nan, np.nan))[0]
    base_r2_overall = base_metrics.get("Overall", (np.nan, np.nan))[1]
    
    print(f"  P1:      RMSE = {base_rmse_p1:.4f} m, R² = {base_r2_p1:.4f}")
    print(f"  Pz12:    RMSE = {base_rmse_pz12:.4f} m, R² = {base_r2_pz12:.4f}")
    print(f"  Overall: RMSE = {base_rmse_overall:.4f} m, R² = {base_r2_overall:.4f}")
    
    # Store results
    results = []
    
    total_runs = len(PARAMS_TO_ANALYZE) * len(PERTURBATIONS)
    run_count = 0
    
    print("\n" + "-"*70)
    print(f"Running {total_runs} PERTURBATION RUNS...")
    print("-"*70)
    
    for param in PARAMS_TO_ANALYZE:
        for delta in PERTURBATIONS:
            run_count += 1
            delta_pct = int(delta * 100)
            run_name = f"{param}_{delta_pct:+d}pct"
            
            perturbed_params = BASE_PARAMS.copy()
            perturbed_params[param] = BASE_PARAMS[param] * (1 + delta)
            
            print(f"  [{run_count}/{total_runs}] {run_name}...", end=" ", flush=True)
            
            metrics = run_single_model(perturbed_params, run_name, cellids)
            
            rmse_overall = metrics.get("Overall", (np.nan, np.nan))[0]
            r2_overall = metrics.get("Overall", (np.nan, np.nan))[1]
            rmse_p1 = metrics.get("P1", (np.nan, np.nan))[0]
            r2_p1 = metrics.get("P1", (np.nan, np.nan))[1]
            rmse_pz12 = metrics.get("Pz12", (np.nan, np.nan))[0]
            r2_pz12 = metrics.get("Pz12", (np.nan, np.nan))[1]
            
            # Calculate sensitivity coefficients
            def calc_sensitivity(new_val, base_val, delta):
                if base_val > 0 and not np.isinf(new_val) and not np.isnan(new_val):
                    return ((new_val - base_val) / base_val) / delta
                return np.nan
            
            s_rmse_overall = calc_sensitivity(rmse_overall, base_rmse_overall, delta)
            s_r2_overall = calc_sensitivity(r2_overall, base_r2_overall, delta)
            s_rmse_p1 = calc_sensitivity(rmse_p1, base_rmse_p1, delta)
            s_r2_p1 = calc_sensitivity(r2_p1, base_r2_p1, delta)
            
            def calc_pct_change(new_val, base_val):
                if base_val > 0 and not np.isinf(new_val) and not np.isnan(new_val):
                    return (new_val - base_val) / base_val * 100
                return np.nan
            
            results.append({
                "parameter": param,
                "perturbation": delta,
                "perturbed_value": perturbed_params[param],
                "rmse_p1": rmse_p1,
                "r2_p1": r2_p1,
                "rmse_pz12": rmse_pz12,
                "r2_pz12": r2_pz12,
                "rmse_overall": rmse_overall,
                "r2_overall": r2_overall,
                "delta_rmse_overall_pct": calc_pct_change(rmse_overall, base_rmse_overall),
                "delta_r2_overall_pct": calc_pct_change(r2_overall, base_r2_overall),
                "delta_rmse_p1_pct": calc_pct_change(rmse_p1, base_rmse_p1),
                "delta_r2_p1_pct": calc_pct_change(r2_p1, base_r2_p1),
                "sensitivity_rmse_overall": s_rmse_overall,
                "sensitivity_r2_overall": s_r2_overall,
                "sensitivity_rmse_p1": s_rmse_p1,
                "sensitivity_r2_p1": s_r2_p1,
            })
            
            print(f"RMSE={rmse_overall:.4f}, R²={r2_overall:.4f}")
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    
    # Add base case
    base_row = {
        "parameter": "BASE",
        "perturbation": 0,
        "perturbed_value": np.nan,
        "rmse_p1": base_rmse_p1,
        "r2_p1": base_r2_p1,
        "rmse_pz12": base_rmse_pz12,
        "r2_pz12": base_r2_pz12,
        "rmse_overall": base_rmse_overall,
        "r2_overall": base_r2_overall,
        "delta_rmse_overall_pct": 0,
        "delta_r2_overall_pct": 0,
        "delta_rmse_p1_pct": 0,
        "delta_r2_p1_pct": 0,
        "sensitivity_rmse_overall": np.nan,
        "sensitivity_r2_overall": np.nan,
        "sensitivity_rmse_p1": np.nan,
        "sensitivity_r2_p1": np.nan,
    }
    results_df = pd.concat([pd.DataFrame([base_row]), results_df], ignore_index=True)
    
    # Save results
    results_csv = output_dir / "sensitivity_results.csv"
    results_df.to_csv(results_csv, index=False)
    print(f"\n✓ Results saved to: {results_csv}")
    
    # Calculate average sensitivity per parameter
    print("\n" + "="*70)
    print("SENSITIVITY SUMMARY (Average |S| across perturbations)")
    print("="*70)
    
    summary_data = []
    for param in PARAMS_TO_ANALYZE:
        param_df = results_df[results_df["parameter"] == param]
        avg_s_rmse_overall = np.nanmean(np.abs(param_df["sensitivity_rmse_overall"]))
        avg_s_r2_overall = np.nanmean(np.abs(param_df["sensitivity_r2_overall"]))
        avg_s_rmse_p1 = np.nanmean(np.abs(param_df["sensitivity_rmse_p1"]))
        avg_s_r2_p1 = np.nanmean(np.abs(param_df["sensitivity_r2_p1"]))
        
        summary_data.append({
            "Parameter": PARAM_DISPLAY_NAMES.get(param, param),
            "param_key": param,
            "|S| RMSE Overall": avg_s_rmse_overall,
            "|S| R² Overall": avg_s_r2_overall,
            "|S| RMSE P1": avg_s_rmse_p1,
            "|S| R² P1": avg_s_r2_p1,
            "Avg |S|": np.nanmean([avg_s_rmse_overall, avg_s_r2_overall, avg_s_rmse_p1, avg_s_r2_p1]),
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values("Avg |S|", ascending=False)
    
    print(f"\n{'Parameter':<25} {'|S| RMSE':>12} {'|S| R²':>12} {'|S| RMSE P1':>14} {'|S| R² P1':>12} {'Avg |S|':>10}")
    print("-"*90)
    for _, row in summary_df.iterrows():
        print(f"{row['Parameter']:<25} {row['|S| RMSE Overall']:>12.3f} {row['|S| R² Overall']:>12.3f} {row['|S| RMSE P1']:>14.3f} {row['|S| R² P1']:>12.3f} {row['Avg |S|']:>10.3f}")
    
    summary_csv = output_dir / "sensitivity_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"\n✓ Summary saved to: {summary_csv}")
    
    return results_df, summary_df, base_metrics


def create_tornado_diagram(results_df: pd.DataFrame, summary_df: pd.DataFrame, output_dir: Path):
    """Create tornado diagram showing sensitivity."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 8))
    
    sorted_params = summary_df.sort_values("Avg |S|", ascending=True)["param_key"].tolist()
    
    colors_neg = '#3498db'
    colors_pos = '#e74c3c'
    
    for idx, (ax, title, metric_col) in enumerate([
        (axes[0], "Overall RMSE", "delta_rmse_overall_pct"),
        (axes[1], "P1 RMSE", "delta_rmse_p1_pct"),
    ]):
        y_pos = np.arange(len(sorted_params))
        
        neg_20 = []
        neg_10 = []
        pos_10 = []
        pos_20 = []
        
        for param in sorted_params:
            param_df = results_df[results_df["parameter"] == param]
            
            for delta in [-0.20, -0.10, 0.10, 0.20]:
                val = param_df[param_df["perturbation"] == delta][metric_col].values
                val = val[0] if len(val) > 0 and not np.isnan(val[0]) and not np.isinf(val[0]) else 0
                
                if delta == -0.20:
                    neg_20.append(val)
                elif delta == -0.10:
                    neg_10.append(val)
                elif delta == 0.10:
                    pos_10.append(val)
                elif delta == 0.20:
                    pos_20.append(val)
        
        bar_height = 0.35
        
        ax.barh(y_pos - bar_height/2, neg_20, height=bar_height, color=colors_neg, alpha=0.9, label='-20%')
        ax.barh(y_pos - bar_height/2, neg_10, height=bar_height, color=colors_neg, alpha=0.5, label='-10%')
        ax.barh(y_pos + bar_height/2, pos_10, height=bar_height, color=colors_pos, alpha=0.5, label='+10%')
        ax.barh(y_pos + bar_height/2, pos_20, height=bar_height, color=colors_pos, alpha=0.9, label='+20%')
        
        ax.axvline(x=0, color='black', linewidth=0.8)
        ax.set_yticks(y_pos)
        ax.set_yticklabels([PARAM_DISPLAY_NAMES.get(p, p) for p in sorted_params])
        ax.set_xlabel('Change in RMSE (%)')
        ax.set_title(f'{title} Sensitivity', fontweight='bold')
        ax.legend(loc='lower right', fontsize=9)
        ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    
    output_path = output_dir / "tornado_diagram.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ Tornado diagram saved to: {output_path}")
    
    plt.show()


def create_spider_plot(results_df: pd.DataFrame, output_dir: Path):
    """Create spider plot showing parameter sensitivity."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    perturbation_levels = [-20, -10, 10, 20]
    
    for idx, (ax, title, metric_col) in enumerate([
        (axes[0], "Overall RMSE Response", "delta_rmse_overall_pct"),
        (axes[1], "P1 RMSE Response", "delta_rmse_p1_pct"),
    ]):
        for param in PARAMS_TO_ANALYZE:
            param_df = results_df[results_df["parameter"] == param].sort_values("perturbation")
            x_vals = [int(p * 100) for p in param_df["perturbation"]]
            y_vals = param_df[metric_col].values
            
            # Filter out invalid values
            valid_mask = ~np.isnan(y_vals) & ~np.isinf(y_vals)
            if valid_mask.sum() > 0:
                x_vals = [x_vals[i] for i in range(len(x_vals)) if valid_mask[i]]
                y_vals = y_vals[valid_mask]
                
                ax.plot(x_vals, y_vals, marker='o', markersize=5, linewidth=1.5, 
                       label=PARAM_DISPLAY_NAMES.get(param, param), alpha=0.8)
        
        ax.axhline(y=0, color='black', linewidth=1, linestyle='--', alpha=0.5)
        ax.axvline(x=0, color='black', linewidth=1, linestyle='--', alpha=0.5)
        ax.set_xlabel('Parameter Perturbation (%)')
        ax.set_ylabel('Change in RMSE (%)')
        ax.set_title(title, fontweight='bold')
        ax.legend(loc='best', fontsize=7, ncol=2)
        ax.grid(True, alpha=0.3)
        ax.set_xticks(perturbation_levels)
    
    plt.tight_layout()
    
    output_path = output_dir / "spider_plot.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ Spider plot saved to: {output_path}")
    
    plt.show()


def create_sensitivity_ranking(summary_df: pd.DataFrame, output_dir: Path):
    """Create bar chart of sensitivity ranking."""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    sorted_df = summary_df.sort_values("Avg |S|", ascending=True)
    
    y_pos = np.arange(len(sorted_df))
    
    colors = plt.cm.RdYlGn_r(np.linspace(0.2, 0.8, len(sorted_df)))
    
    bars = ax.barh(y_pos, sorted_df["Avg |S|"], color=colors, edgecolor='black', linewidth=0.5)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sorted_df["Parameter"])
    ax.set_xlabel('Average Sensitivity Coefficient |S|', fontsize=11)
    ax.set_title('Storage Parameter Sensitivity Ranking (Transient)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    
    for bar, val in zip(bars, sorted_df["Avg |S|"]):
        ax.text(bar.get_width() + 0.005, bar.get_y() + bar.get_height()/2, 
               f'{val:.3f}', va='center', fontsize=9)
    
    plt.tight_layout()
    
    output_path = output_dir / "sensitivity_ranking.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ Sensitivity ranking saved to: {output_path}")
    
    plt.show()


# ============================================================================
# Main
# ============================================================================

def main():
    start_time = time.time()
    
    results_df, summary_df, base_metrics = perform_sensitivity_analysis()
    
    output_dir = Path(CONFIG["outputs_dir"])
    
    print("\n" + "="*70)
    print("CREATING PLOTS...")
    print("="*70)
    
    create_tornado_diagram(results_df, summary_df, output_dir)
    create_spider_plot(results_df, output_dir)
    create_sensitivity_ranking(summary_df, output_dir)
    
    elapsed = time.time() - start_time
    print(f"\n✓ Sensitivity analysis completed in {elapsed/60:.1f} minutes")
    print(f"✓ All outputs saved to: {output_dir}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
