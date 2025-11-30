#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MF6 Single Test Run (Manual Parameter Input)
================================================
Based on Steady_cali_06.py, but uses manual parameter input instead of LHS sampling.

Dependencies: numpy, pandas, flopy
"""

from __future__ import annotations
import os
import sys
import shutil
import glob
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

import flopy
from flopy.mf6 import MFSimulation

# ------------------------
# User Configuration
# ------------------------
ROOT = Path(__file__).resolve().parent
CONFIG = {
    "base_model_ws": str((ROOT.parent / "Model").resolve()),
    # "base_model_obs": str((ROOT.parent / "Out" / "run_modelmuse").resolve()),
    "base_model_obs": str((ROOT.parent / "Final Model").resolve()),
    "mf6_exe": r"D:\Program Files\USGS\mf6.6.2_win64\bin\mf6.exe",
    "test_run_ws": str((ROOT.parent / "Out" / "run_flopy_shp").resolve()),
    "outputs_dir": str((ROOT.parent / "Out").resolve()),
    "obs_csv_name_hint": "Steady.ob_gw_out_head",
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

# ====================
# 👇 Manual Parameter Settings
# ====================
# clean_gravel  >  sandy_gravel  >  loamy_sand  >  clayey_sand  >  limestone  >  sandstone
MANUAL_PARAMS = {   
    # Shapefile lithology parameters (m/s)
    "clean_gravel": 0.06,      # Layer 1: Clean_Gravel
    "clayey_sand": 0.0000137,     # Layer 1: Clayey_Sand
    "loamy_sand": 0.000015,      # Layer 1: Loamy_Sand
    "sandstone": 0.0000411,         # Layer 2: Sand_Stone
    
    # Default formula parameters (m/s)
    "sandy_gravel": 0.0009,     # Layer 1 default value
    "limestone": 0.00003,       # Layer 2 default value
    
    # Boundary condition parameters
    "riv_cond": 0.00025,       # River_Boundary_Whole RIV conductance (m2/s)
    "initial_head": 37.4,      # Initial head (m)
    "ghb_north_head": 37.2,   # North_Boundary GHB boundary head (m)
    "ghb_south_head": 37.6,   # South_Boundary GHB boundary head (m)
}
# 0.956, 0.866
# MANUAL_PARAMS = {   
#     # Shapefile lithology parameters (m/s)
#     "clean_gravel": 0.1,      # Layer 1: Clean_Gravel
#     "clayey_sand": 0.0000537,     # Layer 1: Clayey_Sand
#     "loamy_sand": 0.0000035,      # Layer 1: Loamy_Sand
#     "sandstone": 0.0000411,         # Layer 2: Sand_Stone
    
#     # Default formula parameters (m/s)
#     "sandy_gravel": 0.0009,     # Layer 1 default value
#     "limestone": 0.00003,       # Layer 2 default value
    
#     # Boundary condition parameters
#     "riv_cond": 0.00025,       # River_Boundary_Whole RIV conductance (m2/s)
#     "initial_head": 37.4,      # Initial head (m)
#     "ghb_north_head": 37.2,   # North_Boundary GHB boundary head (m)
#     "ghb_south_head": 37.6,   # South_Boundary GHB boundary head (m)
# }

LITHOLOGY_TARGETS = [
    "clean_gravel",
    "clayey_sand",
    "loamy_sand",
    "sandstone",
]

DEFAULT_PARAMS = [
    "sandy_gravel",
    "limestone",
]

BOUNDARY_PARAMS = [
    "riv_cond",
    "initial_head",
    "ghb_north_head",
    "ghb_south_head",
]

# ------------------------
# Helper Functions (copied from Steady_cali_06.py)
# ------------------------

def load_cellids_by_lithology(sim: MFSimulation) -> Dict[str, List[Tuple[int, int, int]]]:
    """Generate (lay,row,col) cellid list for each lithology based on Shapefile polygon intersection."""
    from flopy.utils.gridintersect import GridIntersect
    try:
        import fiona
        from shapely.geometry import shape as shp_shape
    except Exception as e:
        raise RuntimeError("fiona and shapely are required: pip install fiona shapely") from e

    gwf = sim.get_model(list(sim.model_names)[0])
    mg = gwf.modelgrid

    if mg.grid_type.lower() != "structured":
        raise NotImplementedError("Current implementation assumes DIS structured grid.")

    shp_cfg = CONFIG.get("lithology_shapes", {})
    shp_path = shp_cfg.get("path")
    name_field = shp_cfg.get("name_field", "NAME")
    value_to_key: Dict[str, str] = shp_cfg.get("value_to_key", {})
    layer_by_key: Dict[str, int] = shp_cfg.get("layer_by_key", {})

    if not shp_path or not os.path.exists(shp_path):
        raise FileNotFoundError(f"Shapefile not found: {shp_path}")

    gi = GridIntersect(mg, method="vertex")
    out: Dict[str, List[Tuple[int,int,int]]] = {k: [] for k in LITHOLOGY_TARGETS}

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
                raise ValueError(f"Missing layer mapping for lithology {key}.")

            res = gi.intersect(poly)
            cell_rcs = []
            try:
                cell_rcs = [tuple(rc) for rc in res['cellids']]
            except Exception:
                cell_rcs = [tuple(r.cellid) for r in res]

            out[key].extend([(lay, r, c) for (r, c) in cell_rcs])

    for k in LITHOLOGY_TARGETS:
        print(f"  - {k}: {len(out.get(k, []))} cells")

    return out


def compare_to_observations(sim_df: pd.DataFrame) -> Dict[str, Tuple[float, float]]:
    """Calculate RMSE and R² between simulated and observed values"""
    obs_path = CONFIG["observed_csv_path"]
    delim = CONFIG["observed_csv_delim"]
    index_col = CONFIG["observed_index_col"]
    conds = CONFIG["conditions"]

    obs_df = pd.read_csv(obs_path, delimiter=delim)
    if index_col in obs_df.columns:
        obs_df = obs_df.set_index(index_col)
    
    obs_df = obs_df.loc[sim_df.index]

    metrics: Dict[str, Tuple[float, float]] = {}
    for condition in conds:
        observed = obs_df[condition].values.astype(float)
        simulated = sim_df[condition].values.astype(float)
        
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
        
        if ss_tot > 1e-10:
            r2 = float(1.0 - ss_res / ss_tot)
        else:
            r2 = float("nan")
        
        metrics[condition] = (rmse, r2)
    
    return metrics


def copy_model_to_run(base_ws: str, run_ws: str):
    """Copy model to run directory"""
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
        except Exception as e:
            print(f"⚠️ Failed to delete old folder: {e}")
            backup = run_ws + "_old_" + str(int(time.time()))
            try:
                os.rename(run_ws, backup)
                print(f"Renamed old folder to: {backup}")
            except Exception as e2:
                raise RuntimeError(f"❌ Cannot process target folder: {e2}") from e2
    
    time.sleep(0.3)
    
    try:
        shutil.copytree(base_ws, run_ws, dirs_exist_ok=True)
    except Exception as e:
        raise RuntimeError(f"❌ Failed to copy model: {e}") from e


def fix_nam_paths(workspace: str):
    """Fix relative paths in .nam files"""
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
    """Load MODFLOW 6 simulation"""
    print(f"Loading simulation...sim_ws = {workspace}")
    fix_nam_paths(workspace)
    sim = flopy.mf6.MFSimulation.load(sim_ws=workspace, exe_name=mf6_exe)
    return sim


def set_k_for_zones(sim: MFSimulation, kmap: Dict[str, float], cellids: Dict[str, List[Tuple[int,int,int]]]):
    """Update kx by lithology (also ky=kx, kz=kx/10)"""
    gwf = sim.get_model(list(sim.model_names)[0])
    npf = gwf.get_package("npf")

    k_arr = np.array(npf.k.array)
    k33_arr = np.array(npf.k33.array)
    
    # Step 1: Set global default values (by layer)
    if "sandy_gravel" in kmap:
        kx_layer1 = kmap["sandy_gravel"]
        k_arr[0, :, :] = kx_layer1
        k33_arr[0, :, :] = kx_layer1 / 10.0
        print(f"  Set Layer 1 default (sandy_gravel): kx={kx_layer1:.3e}")
    
    if "limestone" in kmap:
        kx_layer2 = kmap["limestone"]
        k_arr[1, :, :] = kx_layer2
        k33_arr[1, :, :] = kx_layer2 / 10.0
        print(f"  Set Layer 2 default (limestone): kx={kx_layer2:.3e}")

    # Step 2: Override specific cells with shapefile-defined lithology values
    for lith, kx in kmap.items():
        if lith in DEFAULT_PARAMS or lith in BOUNDARY_PARAMS:
            continue
        cells = cellids.get(lith, [])
        if not cells:
            continue
        print(f"  Override {lith}: {len(cells)} cells, kx={kx:.3e}")
        for lay, i, j in cells:
            k_arr[lay, i, j] = kx
            k33_arr[lay, i, j] = kx / 10.0

    npf.k.set_data(k_arr)
    npf.k33.set_data(k33_arr)
    
    sim.write_simulation()


def set_boundary_conditions(sim: MFSimulation, params: Dict[str, float]):
    """Set boundary condition parameters"""
    gwf = sim.get_model(list(sim.model_names)[0])
    
    # 1. Set RIV conductance
    if "riv_cond" in params:
        riv = gwf.get_package("riv")
        if riv is not None:
            riv_data = riv.stress_period_data.get_data(0)
            for i in range(len(riv_data)):
                boundname = str(riv_data[i]['boundname']).lower().strip("'\"")
                if 'river_boundary_whole' in boundname:
                    riv_data[i]['cond'] = params["riv_cond"]
            riv.stress_period_data.set_data(riv_data, 0)
            print(f"  Set RIV conductance: {params['riv_cond']:.3e}")
    
    # 2. Set initial head
    if "initial_head" in params:
        ic = gwf.get_package("ic")
        if ic is not None:
            strt_arr = np.array(ic.strt.array)
            strt_arr[:, :, :] = params["initial_head"]
            ic.strt.set_data(strt_arr)
            print(f"  Set initial head: {params['initial_head']:.3f} m")
    
    # 3. Set GHB (North_Boundary and South_Boundary)
    ghb = gwf.get_package("ghb")
    if ghb is not None:
        ghb_data = ghb.stress_period_data.get_data(0)
        sandy_gravel_kx = params.get("sandy_gravel", 0.0005)
        sandstone_kx = params.get("sandstone", 1e-5)
        limestone_kx = params.get("limestone", 3e-5)
        
        layer1_north_count = 0
        layer1_south_count = 0
        layer2_north_count = 0
        layer2_south_count = 0
        
        for i in range(len(ghb_data)):
            cellid = ghb_data[i]['cellid']
            layer, row, col = cellid[0], cellid[1], cellid[2]
            
            if layer == 0:
                # Layer 1: Use sandy_gravel kx
                ghb_data[i]['cond'] = sandy_gravel_kx
                
                if row == 0 and "ghb_north_head" in params:
                    ghb_data[i]['bhead'] = params["ghb_north_head"]
                    layer1_north_count += 1
                elif row == 126 and col <= 91 and "ghb_south_head" in params:
                    ghb_data[i]['bhead'] = params["ghb_south_head"]
                    layer1_south_count += 1
            elif layer == 1:
                # Layer 2: South boundary uses sandstone, North boundary uses limestone
                if row == 0:
                    ghb_data[i]['cond'] = limestone_kx
                    if "ghb_north_head" in params:
                        ghb_data[i]['bhead'] = params["ghb_north_head"]
                    layer2_north_count += 1
                elif row == 126 and col <= 91:
                    ghb_data[i]['cond'] = sandstone_kx
                    if "ghb_south_head" in params:
                        ghb_data[i]['bhead'] = params["ghb_south_head"]
                    layer2_south_count += 1
        
        ghb.stress_period_data.set_data(ghb_data, 0)
        
        if layer1_north_count > 0:
            print(f"  Set Layer 1 North GHB (row=0): {layer1_north_count} cells, head={params.get('ghb_north_head', 'N/A'):.3f} m, cond={sandy_gravel_kx:.3e}")
        if layer1_south_count > 0:
            print(f"  Set Layer 1 South GHB (row=126, col<=91): {layer1_south_count} cells, head={params.get('ghb_south_head', 'N/A'):.3f} m, cond={sandy_gravel_kx:.3e}")
        if layer2_north_count > 0:
            print(f"  Set Layer 2 North GHB (row=0): {layer2_north_count} cells, head={params.get('ghb_north_head', 'N/A'):.3f} m, cond={limestone_kx:.3e}")
        if layer2_south_count > 0:
            print(f"  Set Layer 2 South GHB (row=126, col<=91): {layer2_south_count} cells, head={params.get('ghb_south_head', 'N/A'):.3f} m, cond={sandstone_kx:.3e}")
    
    sim.write_simulation()


def run_simulation(sim: MFSimulation, run_ws: str) -> bool:
    """Run MODFLOW 6 simulation"""
    success, buff = sim.run_simulation(report=True, silent=True)
    if not success:
        print("Run failed, output:\n" + "\n".join(buff))
    return bool(success)


def find_obs_csv(run_ws: str, hint: str | None = None) -> str | None:
    """Find observation CSV file"""
    csvs = glob.glob(os.path.join(run_ws, "**", "*.csv"), recursive=True)
    if not csvs:
        return None

    outputs_dir = CONFIG.get("outputs_dir", "")
    observed_csv_path = CONFIG.get("observed_csv_path", "")
    observed_csv_path = os.path.normpath(observed_csv_path) if observed_csv_path else ""

    candidates = []
    for p in csvs:
        npth = os.path.normpath(p)
        if outputs_dir:
            abs_p = os.path.abspath(npth)
            abs_out = os.path.abspath(outputs_dir)
            if abs_p.startswith(abs_out + os.sep):
                rel_path = os.path.relpath(abs_p, abs_out)
                if not rel_path.startswith("test_run" + os.sep):
                    continue
        if observed_csv_path and os.path.abspath(npth) == os.path.abspath(observed_csv_path):
            continue
        base = os.path.basename(npth).lower()
        if base in {"results_running.csv", "results.csv", "results.xlsx"}:
            continue
        candidates.append(npth)

    if not candidates:
        return None

    def latest(paths):
        return max(paths, key=lambda x: os.path.getmtime(x))

    if hint:
        hint_lower = hint.lower()
        exacts = [p for p in candidates if os.path.basename(p).lower() in [hint_lower, f"{hint_lower}.csv"]]
        if exacts:
            return latest(exacts)
        matches = [p for p in candidates if hint_lower in os.path.basename(p).lower()]
        if matches:
            return latest(matches)

    return latest(candidates)


def make_sim_df(csv_path: str) -> pd.DataFrame:
    """Read MF6 OBS CSV and construct DataFrame"""
    from flopy.utils.observationfile import Mf6Obs
    conds = CONFIG["conditions"]
    labels = CONFIG["obs_index_labels"]

    sim_out = Mf6Obs(csv_path, isBinary=False)
    df = sim_out.get_dataframe()
    
    df_t = df.transpose()
    df_t = df_t.iloc[1:]
    df_t = df_t[~df_t.index.duplicated(keep='first')]
    
    if len(df_t.columns) >= len(conds):
        df_t.columns = conds[:len(df_t.columns)]
    
    if len(df_t) != len(labels):
        print(f"⚠️ Warning: OBS data has {len(df_t)} observation points, but config requires {len(labels)}")
    
    df_t.index = labels[:len(df_t)]
    
    return df_t


def plot_comparison(obs_csv_path: str, model_csv_path: str, test_csv_path: str):
    """
    Plot comparison: 4 line plots (2 rows x 2 columns)
    First row: Original model (Model/) natural and pumping
    Second row: Test run (Out/test_run/) natural and pumping
    """
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    
    plt.rcParams['figure.dpi'] = 200
    
    # Read observation data
    obs_df = pd.read_csv(obs_csv_path, delimiter=CONFIG["observed_csv_delim"])
    if CONFIG["observed_index_col"] in obs_df.columns:
        obs_df = obs_df.set_index(CONFIG["observed_index_col"])
    
    # Read both simulated datasets
    model_df = make_sim_df(model_csv_path)
    test_df = make_sim_df(test_csv_path)
    
    # Ensure observation data aligns with simulated data observation point order
    obs_df = obs_df.loc[model_df.index]
    
    labels = model_df.index.tolist()
    x = np.arange(len(labels))
    
    # Extract data
    obs_nat = obs_df['natural'].values
    obs_pmp = obs_df['pumping'].values
    
    model_nat = model_df['natural'].values
    model_pmp = model_df['pumping'].values
    
    test_nat = test_df['natural'].values
    test_pmp = test_df['pumping'].values
    
    # Calculate metrics - Original model
    model_resid_nat = model_nat - obs_nat
    model_rmse_nat = np.sqrt(np.mean(model_resid_nat**2))
    model_ss_res_nat = np.sum(model_resid_nat**2)
    model_ss_tot_nat = np.sum((obs_nat - np.mean(obs_nat))**2)
    model_r2_nat = 1 - (model_ss_res_nat / model_ss_tot_nat)
    
    model_resid_pmp = model_pmp - obs_pmp
    model_rmse_pmp = np.sqrt(np.mean(model_resid_pmp**2))
    model_ss_res_pmp = np.sum(model_resid_pmp**2)
    model_ss_tot_pmp = np.sum((obs_pmp - np.mean(obs_pmp))**2)
    model_r2_pmp = 1 - (model_ss_res_pmp / model_ss_tot_pmp)
    
    # Calculate metrics - Test run
    test_resid_nat = test_nat - obs_nat
    test_rmse_nat = np.sqrt(np.mean(test_resid_nat**2))
    test_ss_res_nat = np.sum(test_resid_nat**2)
    test_ss_tot_nat = np.sum((obs_nat - np.mean(obs_nat))**2)
    test_r2_nat = 1 - (test_ss_res_nat / test_ss_tot_nat)
    
    test_resid_pmp = test_pmp - obs_pmp
    test_rmse_pmp = np.sqrt(np.mean(test_resid_pmp**2))
    test_ss_res_pmp = np.sum(test_resid_pmp**2)
    test_ss_tot_pmp = np.sum((obs_pmp - np.mean(obs_pmp))**2)
    test_r2_pmp = 1 - (test_ss_res_pmp / test_ss_tot_pmp)
    
    # Create figure
    fig = plt.figure(figsize=(14, 10))
    gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[1, 1], width_ratios=[1, 1])
    gs.update(hspace=0.3, wspace=0.25)
    
    # First row: Original model
    ax1 = fig.add_subplot(gs[0, 0])  # Top-left: Original model natural
    ax2 = fig.add_subplot(gs[0, 1])  # Top-right: Original model pumping
    
    # Second row: Test run
    ax3 = fig.add_subplot(gs[1, 0])  # Bottom-left: Test run natural
    ax4 = fig.add_subplot(gs[1, 1])  # Bottom-right: Test run pumping
    
    # ---------- Subplot 1: Original Model Natural ----------
    ax1.plot(x, obs_nat, marker='o', markersize=5, linewidth=1.5, label='Observed', color='#1f77b4')
    ax1.plot(x, model_nat, marker='s', markersize=5, linewidth=1.5, label='Simulated (Original)', color='#ff7f0e')
    ax1.set_title('GUI Model - Natural Condition', fontsize=11, fontweight='bold')
    ax1.set_xlabel('Observation Points', fontsize=10)
    ax1.set_ylabel('Head (m)', fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='best', fontsize=9)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=45, fontsize=8, ha='right')
    
    textstr1 = f'R² = {model_r2_nat:.3f}\nRMSE = {model_rmse_nat:.3f} m'
    ax1.text(0.02, 0.98, textstr1, transform=ax1.transAxes, fontsize=9,
             va='top', ha='left',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.75, linewidth=0.5))
    
    # ---------- Subplot 2: Original Model Pumping ----------
    ax2.plot(x, obs_pmp, marker='o', markersize=5, linewidth=1.5, label='Observed', color='#1f77b4')
    ax2.plot(x, model_pmp, marker='s', markersize=5, linewidth=1.5, label='Simulated (Original)', color='#ff7f0e')
    ax2.set_title('GUI - Pumping Condition', fontsize=11, fontweight='bold')
    ax2.set_xlabel('Observation Points', fontsize=10)
    ax2.set_ylabel('Head (m)', fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best', fontsize=9)
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, rotation=45, fontsize=8, ha='right')
    
    textstr2 = f'R² = {model_r2_pmp:.3f}\nRMSE = {model_rmse_pmp:.3f} m'
    ax2.text(0.02, 0.98, textstr2, transform=ax2.transAxes, fontsize=9,
             va='top', ha='left',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.75, linewidth=0.5))
    
    # ---------- Subplot 3: Test Run Natural ----------
    ax3.plot(x, obs_nat, marker='o', markersize=5, linewidth=1.5, label='Observed', color='#1f77b4')
    ax3.plot(x, test_nat, marker='s', markersize=5, linewidth=1.5, label='Simulated (Test)', color='#2ca02c')
    ax3.set_title('Test Run - Natural Condition', fontsize=11, fontweight='bold')
    ax3.set_xlabel('Observation Points', fontsize=10)
    ax3.set_ylabel('Head (m)', fontsize=10)
    ax3.grid(True, alpha=0.3)
    ax3.legend(loc='best', fontsize=9)
    ax3.set_xticks(x)
    ax3.set_xticklabels(labels, rotation=45, fontsize=8, ha='right')
    
    textstr3 = f'R² = {test_r2_nat:.3f}\nRMSE = {test_rmse_nat:.3f} m'
    ax3.text(0.02, 0.98, textstr3, transform=ax3.transAxes, fontsize=9,
             va='top', ha='left',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.75, linewidth=0.5))
    
    # ---------- Subplot 4: Test Run Pumping ----------
    ax4.plot(x, obs_pmp, marker='o', markersize=5, linewidth=1.5, label='Observed', color='#1f77b4')
    ax4.plot(x, test_pmp, marker='s', markersize=5, linewidth=1.5, label='Simulated (Test)', color='#2ca02c')
    ax4.set_title('Test Run - Pumping Condition', fontsize=11, fontweight='bold')
    ax4.set_xlabel('Observation Points', fontsize=10)
    ax4.set_ylabel('Head (m)', fontsize=10)
    ax4.grid(True, alpha=0.3)
    ax4.legend(loc='best', fontsize=9)
    ax4.set_xticks(x)
    ax4.set_xticklabels(labels, rotation=45, fontsize=8, ha='right')
    
    textstr4 = f'R² = {test_r2_pmp:.3f}\nRMSE = {test_rmse_pmp:.3f} m'
    ax4.text(0.02, 0.98, textstr4, transform=ax4.transAxes, fontsize=9,
             va='top', ha='left',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.75, linewidth=0.5))
    
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.join(CONFIG["outputs_dir"], "comparison_plot.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved to: {output_path}")
    
    plt.show()


# ------------------------
# Main Workflow
# ------------------------

def main():
    base_ws = CONFIG["base_model_ws"]
    test_run_ws = CONFIG["test_run_ws"]
    outputs_dir = CONFIG["outputs_dir"]
    mf6_exe = CONFIG["mf6_exe"]
    obs_hint = CONFIG.get("obs_csv_name_hint")

    Path(test_run_ws).parent.mkdir(parents=True, exist_ok=True)
    Path(outputs_dir).mkdir(parents=True, exist_ok=True)

    print("\n=== Manual Parameter Test Run ===")
    print("\nCurrent parameter settings:")
    for k, v in MANUAL_PARAMS.items():
        if "head" in k or k == "initial_head":
            print(f"  {k}: {v:.3f}")
        else:
            print(f"  {k}: {v:.3e}")

    # 1. Copy model to test directory
    print(f"\nCopying model to: {test_run_ws}")
    copy_model_to_run(base_ws, test_run_ws)

    # 2. Load simulation
    sim = load_simulation(test_run_ws, mf6_exe)
    
    # 3. Load cellids
    print("\nLoading lithology-cell mapping:")
    cellids = load_cellids_by_lithology(sim)

    # 4. Set parameters
    print("\nSetting hydraulic conductivity:")
    set_k_for_zones(sim, MANUAL_PARAMS, cellids)
    
    print("\nSetting boundary conditions:")
    set_boundary_conditions(sim, MANUAL_PARAMS)

    # 5. Run simulation
    print("\nRunning MODFLOW 6...")
    t0 = time.time()
    ok = run_simulation(sim, test_run_ws)
    dt = time.time() - t0

    # 6. Analyze results
    if ok:
        print(f"✅ Run successful! Time elapsed: {dt:.2f} seconds")
        
        obs_csv = find_obs_csv(test_run_ws, hint=obs_hint)
        if obs_csv and os.path.exists(obs_csv):
            print(f"\nReading observation data: {obs_csv}")
            sim_df = make_sim_df(obs_csv)
            
            # Check for failed observation points
            has_failure = False
            for condition in CONFIG["conditions"]:
                if condition in sim_df.columns:
                    simulated = sim_df[condition].values.astype(float)
                    if (simulated < -1e+20).any():
                        print(f"⚠️ {condition} condition has failed observation points")
                        has_failure = True
            
            if not has_failure:
                metrics = compare_to_observations(sim_df)
                
                print("\n=== Model Performance Metrics ===")
                for condition, (rmse, r2) in metrics.items():
                    print(f"{condition}:")
                    print(f"  RMSE: {rmse:.4f} m")
                    print(f"  R²:   {r2:.4f}")
                
                # Display simulated vs observed values
                print("\n=== Simulated vs Observed Values ===")
                obs_path = CONFIG["observed_csv_path"]
                obs_df = pd.read_csv(obs_path, delimiter=CONFIG["observed_csv_delim"])
                if CONFIG["observed_index_col"] in obs_df.columns:
                    obs_df = obs_df.set_index(CONFIG["observed_index_col"])
                obs_df = obs_df.loc[sim_df.index]
                
                for condition in CONFIG["conditions"]:
                    print(f"\n{condition}:")
                    print(f"{'Obs Point':<8} {'Observed(m)':<12} {'Simulated(m)':<12} {'Diff(m)':<10}")
                    print("-" * 45)
                    for idx in sim_df.index:
                        obs_val = obs_df.loc[idx, condition]
                        sim_val = sim_df.loc[idx, condition]
                        diff = sim_val - obs_val
                        print(f"{idx:<8} {obs_val:<12.3f} {sim_val:<12.3f} {diff:<10.3f}")
                
                # Plot comparison
                print("\n=== Plotting Comparison ===")
                model_csv = os.path.join(CONFIG["base_model_obs"], "Steady_cali.ob_gw_out_head.csv")
                if os.path.exists(model_csv):
                    plot_comparison(
                        obs_csv_path=CONFIG["observed_csv_path"],
                        model_csv_path=model_csv,
                        test_csv_path=obs_csv
                    )
                else:
                    print(f"⚠️ Original model output file not found: {model_csv}")
                    print("   Skipping comparison plot")
            else:
                print("\n❌ Model has failed observation points, cannot calculate metrics")
        else:
            print("\n⚠️ Observation data CSV file not found")
    else:
        print(f"❌ Run failed! Time elapsed: {dt:.2f} seconds")
        print("Please check MODFLOW output files and logs.")

    print(f"\nRun results saved to: {test_run_ws}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error occurred: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
