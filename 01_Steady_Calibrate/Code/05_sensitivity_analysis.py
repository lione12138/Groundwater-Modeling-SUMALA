#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Steady-State Sensitivity Analysis
=================================
Performs One-at-a-Time (OAT) sensitivity analysis on calibrated parameters.
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
from pathlib import Path
from typing import Dict, List, Tuple
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('TkAgg')

import flopy
from flopy.mf6 import MFSimulation

# ============================================================================
# Configuration
# ============================================================================
ROOT = Path(__file__).resolve().parent
CONFIG = {
    "base_model_ws": str((ROOT.parent / "Model").resolve()),
    "mf6_exe": r"D:\Program Files\USGS\mf6.6.2_win64\bin\mf6.exe",
    "sensitivity_ws": str((ROOT.parent / "Out" / "Sensitivity Ana").resolve()),
    "outputs_dir": str((ROOT.parent / "Out" / "Sensitivity Ana").resolve()),
    "obs_csv_name_hint": "Steady_cali.ob_gw_out_head",
    "observed_csv_path": str((ROOT.parent / "steady_state_calibration.csv").resolve()),
    "observed_csv_delim": ";",
    "observed_index_col": "Unnamed: 0",
    "conditions": ["natural", "pumping"],
    "obs_index_labels": ["P1", "Pz2", "Pz3", "Pz4", "Pz5", "Pz6", "Pz7", "Pz8", "Pz9", "Pz10", "Pz11", "Pz12"],
    "lithology_shapes": {
        "path": str((ROOT.parent / "Lithology" / "KX-new.shp").resolve()),
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
    # Shapefile lithology parameters (m/s)
    "clean_gravel": 0.06,
    "clayey_sand": 0.0000137,
    "loamy_sand": 0.000015,
    "sandstone": 0.0000411,
    
    # Default lithology parameters (m/s)
    "sandy_gravel": 0.0009,
    "limestone": 0.00003,
    
    # Boundary condition parameters
    "riv_cond": 0.00025,
    "initial_head": 37.4,
    "ghb_north_head": 37.2,
    "ghb_south_head": 37.6,
}

# Parameters to analyze (all 10 calibration parameters)
PARAMS_TO_ANALYZE = [
    "clean_gravel",
    "clayey_sand", 
    "loamy_sand",
    "sandstone",
    "sandy_gravel",
    "limestone",
    "riv_cond",
    "initial_head",
    "ghb_north_head",
    "ghb_south_head",
]

# Perturbation levels
PERTURBATIONS = [-0.20, -0.10, 0.10, 0.20]  # ±10%, ±20%

# Parameter categories for grouping
LITHOLOGY_TARGETS = ["clean_gravel", "clayey_sand", "loamy_sand", "sandstone"]
DEFAULT_PARAMS = ["sandy_gravel", "limestone"]
BOUNDARY_PARAMS = ["riv_cond", "initial_head", "ghb_north_head", "ghb_south_head"]

# Display names for plots
PARAM_DISPLAY_NAMES = {
    "clean_gravel": "K (Clean Gravel)",
    "clayey_sand": "K (Clayey Sand)",
    "loamy_sand": "K (Loamy Sand)",
    "sandstone": "K (Sandstone)",
    "sandy_gravel": "K (Sandy Gravel)",
    "limestone": "K (Limestone)",
    "riv_cond": "River Conductance",
    "initial_head": "Initial Head",
    "ghb_north_head": "GHB North Head",
    "ghb_south_head": "GHB South Head",
}

# ============================================================================
# Model Functions (from Steady_test.py)
# ============================================================================

def load_cellids_by_lithology(sim: MFSimulation) -> Dict[str, List[Tuple[int, int, int]]]:
    """Load cell IDs for each lithology from shapefile."""
    from flopy.utils.gridintersect import GridIntersect
    import fiona
    from shapely.geometry import shape as shp_shape

    gwf = sim.get_model(list(sim.model_names)[0])
    mg = gwf.modelgrid

    shp_cfg = CONFIG.get("lithology_shapes", {})
    shp_path = shp_cfg.get("path")
    name_field = shp_cfg.get("name_field", "NAME")
    value_to_key = shp_cfg.get("value_to_key", {})
    layer_by_key = shp_cfg.get("layer_by_key", {})

    gi = GridIntersect(mg, method="vertex")
    out = {k: [] for k in LITHOLOGY_TARGETS}

    with fiona.open(shp_path, "r") as src:
        for feat in src:
            props = feat.get("properties") or {}
            geom = feat.get("geometry")
            if not geom:
                continue
            poly = shp_shape(geom)
            raw_val = props.get(name_field)
            if raw_val is None:
                continue
            key = value_to_key.get(str(raw_val))
            if key is None:
                key = str(raw_val)
            if key not in LITHOLOGY_TARGETS:
                continue
            lay = layer_by_key.get(key)
            if lay is None:
                continue

            res = gi.intersect(poly)
            cell_rcs = []
            try:
                cell_rcs = [tuple(rc) for rc in res['cellids']]
            except Exception:
                cell_rcs = [tuple(r.cellid) for r in res]

            out[key].extend([(lay, r, c) for (r, c) in cell_rcs])

    return out


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


def fix_nam_paths(workspace: str):
    """Fix relative paths in .nam files."""
    import re
    
    mfsim_path = os.path.join(workspace, "mfsim.nam")
    if os.path.exists(mfsim_path):
        with open(mfsim_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        content = re.sub(r'\.\.[\\/]\.\.[\\/]model[\\/]model\.\w+[\\/]', '', content)
        content = re.sub(r'\.\.[\\/]\.\.[\\/]output[\\/]output\.\w+[\\/]', '', content)
        with open(mfsim_path, 'w', encoding='utf-8') as f:
            f.write(content)
    
    for nam_file in glob.glob(os.path.join(workspace, "*.nam")):
        if os.path.basename(nam_file).lower() == "mfsim.nam":
            continue
        with open(nam_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        content = re.sub(r'\.\.[\\/]\.\.[\\/]model[\\/]model\.\w+[\\/]', '', content)
        content = re.sub(r'\.\.[\\/]\.\.[\\/]output[\\/]output\.\w+[\\/]', '', content)
        with open(nam_file, 'w', encoding='utf-8') as f:
            f.write(content)


def load_simulation(workspace: str, mf6_exe: str) -> MFSimulation:
    """Load MODFLOW 6 simulation."""
    fix_nam_paths(workspace)
    sim = flopy.mf6.MFSimulation.load(sim_ws=workspace, exe_name=mf6_exe, verbosity_level=0)
    return sim


def set_k_for_zones(sim: MFSimulation, kmap: Dict[str, float], cellids: Dict[str, List[Tuple[int,int,int]]]):
    """Set hydraulic conductivity for zones."""
    gwf = sim.get_model(list(sim.model_names)[0])
    npf = gwf.get_package("npf")

    k_arr = np.array(npf.k.array)
    k33_arr = np.array(npf.k33.array)
    
    # Set layer defaults
    if "sandy_gravel" in kmap:
        k_arr[0, :, :] = kmap["sandy_gravel"]
        k33_arr[0, :, :] = kmap["sandy_gravel"] / 10.0
    
    if "limestone" in kmap:
        k_arr[1, :, :] = kmap["limestone"]
        k33_arr[1, :, :] = kmap["limestone"] / 10.0

    # Override with shapefile zones
    for lith, kx in kmap.items():
        if lith in DEFAULT_PARAMS or lith in BOUNDARY_PARAMS:
            continue
        cells = cellids.get(lith, [])
        for lay, i, j in cells:
            k_arr[lay, i, j] = kx
            k33_arr[lay, i, j] = kx / 10.0

    npf.k.set_data(k_arr)
    npf.k33.set_data(k33_arr)
    sim.write_simulation()


def set_boundary_conditions(sim: MFSimulation, params: Dict[str, float]):
    """Set boundary condition parameters."""
    gwf = sim.get_model(list(sim.model_names)[0])
    
    # RIV conductance
    if "riv_cond" in params:
        riv = gwf.get_package("riv")
        if riv is not None:
            riv_data = riv.stress_period_data.get_data(0)
            for i in range(len(riv_data)):
                boundname = str(riv_data[i]['boundname']).lower().strip("'\"")
                if 'river_boundary_whole' in boundname:
                    riv_data[i]['cond'] = params["riv_cond"]
            riv.stress_period_data.set_data(riv_data, 0)
    
    # Initial head
    if "initial_head" in params:
        ic = gwf.get_package("ic")
        if ic is not None:
            strt_arr = np.array(ic.strt.array)
            strt_arr[:, :, :] = params["initial_head"]
            ic.strt.set_data(strt_arr)
    
    # GHB
    ghb = gwf.get_package("ghb")
    if ghb is not None:
        ghb_data = ghb.stress_period_data.get_data(0)
        sandy_gravel_kx = params.get("sandy_gravel", 0.0005)
        sandstone_kx = params.get("sandstone", 1e-5)
        limestone_kx = params.get("limestone", 3e-5)
        
        for i in range(len(ghb_data)):
            cellid = ghb_data[i]['cellid']
            layer, row, col = cellid[0], cellid[1], cellid[2]
            
            if layer == 0:
                ghb_data[i]['cond'] = sandy_gravel_kx
                if row == 0 and "ghb_north_head" in params:
                    ghb_data[i]['bhead'] = params["ghb_north_head"]
                elif row == 126 and col <= 91 and "ghb_south_head" in params:
                    ghb_data[i]['bhead'] = params["ghb_south_head"]
            elif layer == 1:
                if row == 0:
                    ghb_data[i]['cond'] = limestone_kx
                    if "ghb_north_head" in params:
                        ghb_data[i]['bhead'] = params["ghb_north_head"]
                elif row == 126 and col <= 91:
                    ghb_data[i]['cond'] = sandstone_kx
                    if "ghb_south_head" in params:
                        ghb_data[i]['bhead'] = params["ghb_south_head"]
        
        ghb.stress_period_data.set_data(ghb_data, 0)
    
    sim.write_simulation()


def run_simulation(sim: MFSimulation) -> bool:
    """Run MODFLOW 6 simulation."""
    success, _ = sim.run_simulation(report=False, silent=True)
    return bool(success)


def find_obs_csv(run_ws: str) -> str | None:
    """Find observation CSV file."""
    csvs = glob.glob(os.path.join(run_ws, "*.csv"))
    for p in csvs:
        if "ob_gw_out_head" in os.path.basename(p).lower():
            return p
    return None


def calculate_metrics(obs_csv: str) -> Dict[str, Tuple[float, float]]:
    """Calculate RMSE and R² for each condition."""
    from flopy.utils.observationfile import Mf6Obs
    
    obs_path = CONFIG["observed_csv_path"]
    delim = CONFIG["observed_csv_delim"]
    index_col = CONFIG["observed_index_col"]
    conds = CONFIG["conditions"]
    labels = CONFIG["obs_index_labels"]

    # Read observed data
    obs_df = pd.read_csv(obs_path, delimiter=delim)
    if index_col in obs_df.columns:
        obs_df = obs_df.set_index(index_col)
    
    # Read simulated data
    sim_out = Mf6Obs(obs_csv, isBinary=False)
    df = sim_out.get_dataframe()
    df_t = df.transpose()
    df_t = df_t.iloc[1:]
    df_t = df_t[~df_t.index.duplicated(keep='first')]
    
    if len(df_t.columns) >= len(conds):
        df_t.columns = conds[:len(df_t.columns)]
    df_t.index = labels[:len(df_t)]
    
    obs_df = obs_df.loc[df_t.index]

    metrics = {}
    for condition in conds:
        observed = obs_df[condition].values.astype(float)
        simulated = df_t[condition].values.astype(float)
        
        valid_mask = (simulated > -1e+20)
        if valid_mask.sum() == 0:
            metrics[condition] = (float("inf"), float("nan"))
            continue
        
        observed = observed[valid_mask]
        simulated = simulated[valid_mask]
        
        residuals = simulated - observed
        rmse = float(np.sqrt(np.mean(residuals ** 2)))
        
        ss_res = float(np.sum(residuals ** 2))
        ss_tot = float(np.sum((observed - np.mean(observed)) ** 2))
        
        r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 1e-10 else float("nan")
        
        metrics[condition] = (rmse, r2)
    
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
    
    # Copy model
    copy_model_to_run(base_ws, run_ws)
    
    # Load and modify
    sim = load_simulation(run_ws, mf6_exe)
    set_k_for_zones(sim, params, cellids)
    set_boundary_conditions(sim, params)
    
    # Run
    success = run_simulation(sim)
    
    if not success:
        return {"natural": (float("inf"), float("nan")), "pumping": (float("inf"), float("nan"))}
    
    # Get metrics
    obs_csv = find_obs_csv(run_ws)
    if obs_csv:
        return calculate_metrics(obs_csv)
    else:
        return {"natural": (float("inf"), float("nan")), "pumping": (float("inf"), float("nan"))}


def perform_sensitivity_analysis():
    """Perform OAT sensitivity analysis."""
    print("\n" + "="*70)
    print("STEADY-STATE SENSITIVITY ANALYSIS")
    print("="*70)
    
    # Create output directory
    output_dir = Path(CONFIG["outputs_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "runs").mkdir(exist_ok=True)
    
    # Load cellids once
    print("\nLoading lithology cell mappings...")
    base_ws = CONFIG["base_model_ws"]
    mf6_exe = CONFIG["mf6_exe"]
    
    temp_ws = str(output_dir / "runs" / "temp_load")
    copy_model_to_run(base_ws, temp_ws)
    sim = load_simulation(temp_ws, mf6_exe)
    cellids = load_cellids_by_lithology(sim)
    
    # Run base case
    print("\n" + "-"*70)
    print("Running BASE CASE...")
    print("-"*70)
    base_metrics = run_single_model(BASE_PARAMS.copy(), "base", cellids)
    
    base_rmse_nat = base_metrics["natural"][0]
    base_r2_nat = base_metrics["natural"][1]
    base_rmse_pump = base_metrics["pumping"][0]
    base_r2_pump = base_metrics["pumping"][1]
    
    print(f"  Natural:  RMSE = {base_rmse_nat:.4f} m, R² = {base_r2_nat:.4f}")
    print(f"  Pumping:  RMSE = {base_rmse_pump:.4f} m, R² = {base_r2_pump:.4f}")
    
    # Store results
    results = []
    
    # Run perturbations
    total_runs = len(PARAMS_TO_ANALYZE) * len(PERTURBATIONS)
    run_count = 0
    
    print("\n" + "-"*70)
    print(f"Running {total_runs} PERTURBATION RUNS...")
    print("-"*70)
    
    for param in PARAMS_TO_ANALYZE:
        param_results = {"parameter": param, "base_value": BASE_PARAMS[param]}
        
        for delta in PERTURBATIONS:
            run_count += 1
            delta_pct = int(delta * 100)
            run_name = f"{param}_{delta_pct:+d}pct"
            
            # Create perturbed parameters
            perturbed_params = BASE_PARAMS.copy()
            perturbed_params[param] = BASE_PARAMS[param] * (1 + delta)
            
            print(f"  [{run_count}/{total_runs}] {run_name}...", end=" ", flush=True)
            
            metrics = run_single_model(perturbed_params, run_name, cellids)
            
            rmse_nat = metrics["natural"][0]
            r2_nat = metrics["natural"][1]
            rmse_pump = metrics["pumping"][0]
            r2_pump = metrics["pumping"][1]
            
            # Calculate sensitivity coefficients
            # S = (ΔOutput/Output_base) / (ΔParam/Param_base) = (ΔOutput/Output_base) / delta
            if base_rmse_nat > 0 and not np.isinf(rmse_nat):
                s_rmse_nat = ((rmse_nat - base_rmse_nat) / base_rmse_nat) / delta
            else:
                s_rmse_nat = np.nan
                
            if base_r2_nat > 0 and not np.isnan(r2_nat):
                s_r2_nat = ((r2_nat - base_r2_nat) / base_r2_nat) / delta
            else:
                s_r2_nat = np.nan
                
            if base_rmse_pump > 0 and not np.isinf(rmse_pump):
                s_rmse_pump = ((rmse_pump - base_rmse_pump) / base_rmse_pump) / delta
            else:
                s_rmse_pump = np.nan
                
            if base_r2_pump > 0 and not np.isnan(r2_pump):
                s_r2_pump = ((r2_pump - base_r2_pump) / base_r2_pump) / delta
            else:
                s_r2_pump = np.nan
            
            results.append({
                "parameter": param,
                "perturbation": delta,
                "perturbed_value": perturbed_params[param],
                "rmse_natural": rmse_nat,
                "r2_natural": r2_nat,
                "rmse_pumping": rmse_pump,
                "r2_pumping": r2_pump,
                "delta_rmse_nat_pct": (rmse_nat - base_rmse_nat) / base_rmse_nat * 100 if base_rmse_nat > 0 else np.nan,
                "delta_r2_nat_pct": (r2_nat - base_r2_nat) / base_r2_nat * 100 if base_r2_nat > 0 else np.nan,
                "delta_rmse_pump_pct": (rmse_pump - base_rmse_pump) / base_rmse_pump * 100 if base_rmse_pump > 0 else np.nan,
                "delta_r2_pump_pct": (r2_pump - base_r2_pump) / base_r2_pump * 100 if base_r2_pump > 0 else np.nan,
                "sensitivity_rmse_nat": s_rmse_nat,
                "sensitivity_r2_nat": s_r2_nat,
                "sensitivity_rmse_pump": s_rmse_pump,
                "sensitivity_r2_pump": s_r2_pump,
            })
            
            print(f"RMSE_nat={rmse_nat:.4f}, R²_nat={r2_nat:.4f}")
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    
    # Add base case
    base_row = {
        "parameter": "BASE",
        "perturbation": 0,
        "perturbed_value": np.nan,
        "rmse_natural": base_rmse_nat,
        "r2_natural": base_r2_nat,
        "rmse_pumping": base_rmse_pump,
        "r2_pumping": base_r2_pump,
        "delta_rmse_nat_pct": 0,
        "delta_r2_nat_pct": 0,
        "delta_rmse_pump_pct": 0,
        "delta_r2_pump_pct": 0,
        "sensitivity_rmse_nat": np.nan,
        "sensitivity_r2_nat": np.nan,
        "sensitivity_rmse_pump": np.nan,
        "sensitivity_r2_pump": np.nan,
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
        avg_s_rmse_nat = np.nanmean(np.abs(param_df["sensitivity_rmse_nat"]))
        avg_s_r2_nat = np.nanmean(np.abs(param_df["sensitivity_r2_nat"]))
        avg_s_rmse_pump = np.nanmean(np.abs(param_df["sensitivity_rmse_pump"]))
        avg_s_r2_pump = np.nanmean(np.abs(param_df["sensitivity_r2_pump"]))
        
        summary_data.append({
            "Parameter": PARAM_DISPLAY_NAMES.get(param, param),
            "param_key": param,
            "|S| RMSE Natural": avg_s_rmse_nat,
            "|S| R² Natural": avg_s_r2_nat,
            "|S| RMSE Pumping": avg_s_rmse_pump,
            "|S| R² Pumping": avg_s_r2_pump,
            "Avg |S|": np.nanmean([avg_s_rmse_nat, avg_s_r2_nat, avg_s_rmse_pump, avg_s_r2_pump]),
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values("Avg |S|", ascending=False)
    
    print(f"\n{'Parameter':<25} {'|S| RMSE Nat':>12} {'|S| R² Nat':>12} {'|S| RMSE Pump':>14} {'|S| R² Pump':>12} {'Avg |S|':>10}")
    print("-"*90)
    for _, row in summary_df.iterrows():
        print(f"{row['Parameter']:<25} {row['|S| RMSE Natural']:>12.3f} {row['|S| R² Natural']:>12.3f} {row['|S| RMSE Pumping']:>14.3f} {row['|S| R² Pumping']:>12.3f} {row['Avg |S|']:>10.3f}")
    
    # Save summary
    summary_csv = output_dir / "sensitivity_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"\n✓ Summary saved to: {summary_csv}")
    
    return results_df, summary_df, base_metrics


def create_tornado_diagram(results_df: pd.DataFrame, summary_df: pd.DataFrame, output_dir: Path):
    """Create tornado diagram showing sensitivity."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 8))
    
    # Exclude GHB head parameters (too sensitive, causes model instability)
    params_to_plot = [p for p in PARAMS_TO_ANALYZE if p not in ["ghb_north_head", "ghb_south_head"]]
    
    # Sort by average sensitivity (excluding GHB heads)
    filtered_summary = summary_df[summary_df["param_key"].isin(params_to_plot)]
    sorted_params = filtered_summary.sort_values("Avg |S|", ascending=True)["param_key"].tolist()
    
    colors_neg = '#3498db'  # Blue for negative perturbation
    colors_pos = '#e74c3c'  # Red for positive perturbation
    
    for idx, (ax, condition, metric_col) in enumerate([
        (axes[0], "natural", "delta_rmse_nat_pct"),
        (axes[1], "pumping", "delta_rmse_pump_pct"),
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
        
        # Plot bars
        bar_height = 0.35
        
        ax.barh(y_pos - bar_height/2, neg_20, height=bar_height, color=colors_neg, alpha=0.9, label='-20%')
        ax.barh(y_pos - bar_height/2, neg_10, height=bar_height, color=colors_neg, alpha=0.5, label='-10%')
        ax.barh(y_pos + bar_height/2, pos_10, height=bar_height, color=colors_pos, alpha=0.5, label='+10%')
        ax.barh(y_pos + bar_height/2, pos_20, height=bar_height, color=colors_pos, alpha=0.9, label='+20%')
        
        ax.axvline(x=0, color='black', linewidth=0.8)
        ax.set_yticks(y_pos)
        ax.set_yticklabels([PARAM_DISPLAY_NAMES.get(p, p) for p in sorted_params])
        ax.set_xlabel('Change in RMSE (%)')
        ax.set_title(f'{condition.capitalize()} Condition - RMSE Sensitivity', fontweight='bold')
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
    
    # Exclude GHB head parameters
    params_to_plot = [p for p in PARAMS_TO_ANALYZE if p not in ["ghb_north_head", "ghb_south_head"]]
    
    for idx, (ax, condition, metric_col) in enumerate([
        (axes[0], "Natural", "delta_rmse_nat_pct"),
        (axes[1], "Pumping", "delta_rmse_pump_pct"),
    ]):
        for param in params_to_plot:
            param_df = results_df[results_df["parameter"] == param].sort_values("perturbation")
            x_vals = [int(p * 100) for p in param_df["perturbation"]]
            y_vals = param_df[metric_col].values
            
            ax.plot(x_vals, y_vals, marker='o', markersize=5, linewidth=1.5, 
                   label=PARAM_DISPLAY_NAMES.get(param, param), alpha=0.8)
        
        ax.axhline(y=0, color='black', linewidth=1, linestyle='--', alpha=0.5)
        ax.axvline(x=0, color='black', linewidth=1, linestyle='--', alpha=0.5)
        ax.set_xlabel('Parameter Perturbation (%)')
        ax.set_ylabel('Change in RMSE (%)')
        ax.set_title(f'{condition} Condition - RMSE Response', fontweight='bold')
        ax.legend(loc='best', fontsize=8, ncol=2)
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
    
    # Exclude GHB head parameters and sort by average sensitivity
    params_to_plot = [p for p in PARAMS_TO_ANALYZE if p not in ["ghb_north_head", "ghb_south_head"]]
    filtered_df = summary_df[summary_df["param_key"].isin(params_to_plot)]
    sorted_df = filtered_df.sort_values("Avg |S|", ascending=True)
    
    y_pos = np.arange(len(sorted_df))
    
    colors = plt.cm.RdYlGn_r(np.linspace(0.2, 0.8, len(sorted_df)))
    
    bars = ax.barh(y_pos, sorted_df["Avg |S|"], color=colors, edgecolor='black', linewidth=0.5)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sorted_df["Parameter"])
    ax.set_xlabel('Average Sensitivity Coefficient |S|', fontsize=11)
    ax.set_title('Parameter Sensitivity Ranking', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    
    # Add value labels
    for bar, val in zip(bars, sorted_df["Avg |S|"]):
        ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2, 
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
    
    # Run sensitivity analysis
    results_df, summary_df, base_metrics = perform_sensitivity_analysis()
    
    output_dir = Path(CONFIG["outputs_dir"])
    
    # Create plots
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
