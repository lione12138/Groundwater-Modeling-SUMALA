#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MF6 Single Test Run (Manual Parameter Input - Based on flopy-strict.py)
================================================
Based on flopy-strict.py, but uses manual parameter input instead of LHS sampling.

Key differences:
- Reads lithology distribution from ModMuse .gpt file (not using shapefile)
- Strictly follows ModMuse's internal logic

Dependencies: numpy, pandas, flopy, shapely
"""

from __future__ import annotations
import os
import sys
import shutil
import glob
import time
import re
from pathlib import Path
from typing import Dict, List, Tuple, Any

import numpy as np
import pandas as pd

import flopy
from flopy.mf6 import MFSimulation
from flopy.discretization import StructuredGrid
from flopy.utils import GridIntersect
from shapely.geometry import Polygon

# ------------------------
# User Configuration
# ------------------------
ROOT = Path(__file__).resolve().parent
CONFIG = {
    "base_model_ws": str((ROOT.parent / "Model").resolve()),
    "mf6_exe": r"D:\Program Files\USGS\mf6.6.2_win64\bin\mf6.exe",
    "test_run_ws": str((ROOT.parent / "Out" / "run_flopy_gpt").resolve()),
    "outputs_dir": str((ROOT.parent / "Out").resolve()),
    "obs_csv_name_hint": "Steady_cali.ob_gw_out_head",
    "observed_csv_path": str((ROOT.parent / "steady_state_calibration.csv").resolve()),
    "observed_csv_delim": ";",
    "observed_index_col": "Unnamed: 0",
    "conditions": ["natural", "pumping"],
    "obs_index_labels": ["P1", "Pz2", "Pz3", "Pz4", "Pz5", "Pz6", "Pz7", "Pz8", "Pz9", "Pz10", "Pz11", "Pz12"],
    "modmuse_gpt_path": str((ROOT.parent / "Model" / "Steady_cali.gpt").resolve()),
}

# ====================
# 👇 Manual Parameter Settings
# ====================
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

DEFAULT_PARAMS = [
    "sandy_gravel",
    "limestone",
]

LITHOLOGY_OBJECTS = [
    "clean_gravel",
    "clayey_sand",
    "loamy_sand",
    "sandstone",
]

BOUNDARY_PARAMS = [
    "riv_cond",
    "initial_head",
    "ghb_north_head",
    "ghb_south_head",
]

# ------------------------
# Helper Functions (copied from flopy-strict.py)
# ------------------------

def parse_gpt_screenobjects(gpt_path: str) -> Dict[str, Dict[str, Any]]:
    """Parse ScreenObject definitions from ModMuse .gpt file"""
    result = {}
    
    target_names = {
        "Clean_Gravel": "clean_gravel",
        "Clayey_Sand": "clayey_sand",
        "Loamy_Sand": "loamy_sand",
        "Sand_Stone": "sandstone",
    }
    
    with open(gpt_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i]
        
        found_name = None
        for target_name, normalized_name in target_names.items():
            if f"ScreenObject.Name = '{target_name}'" in line:
                found_name = (target_name, normalized_name)
                break
        
        if not found_name:
            i += 1
            continue
        
        target_name, normalized_name = found_name
        
        item_start = None
        j = i - 1
        while j >= max(0, i - 200):
            if 'item' in lines[j] and lines[j].strip() == 'item':
                has_classtype = False
                for k in range(j, min(j + 15, len(lines))):
                    if "ClassTypeName = 'TScreenObject'" in lines[k]:
                        has_classtype = True
                        break
                if has_classtype:
                    item_start = j
                    break
            j -= 1
        
        if item_start is None:
            i += 1
            continue
        
        item_end = None
        mixture_found = False
        for j in range(i, min(i + 500, len(lines))):
            if "MixtureFormulas.Strings" in lines[j]:
                mixture_found = True
            if mixture_found and lines[j].strip() == 'end':
                item_end = j
                break
        
        if item_end is None:
            i += 1
            continue
        
        block = ''.join(lines[item_start:item_end+1])
        
        # Extract Points
        points = []
        points_match = re.search(r"Points\s*=\s*<(.*?)>", block, re.DOTALL)
        if points_match:
            points_block = points_match.group(1)
            point_items = re.finditer(r"X\s*=\s*([\d.e+-]+)\s+Y\s*=\s*([\d.e+-]+)", points_block)
            for point_match in point_items:
                x = float(point_match.group(1))
                y = float(point_match.group(2))
                points.append((x, y))
        
        # Extract Layer constraint
        higher_elev = None
        lower_elev = None
        higher_match = re.search(r"HigherElevationFormula\s*=\s*'([^']+)'", block)
        lower_match = re.search(r"LowerElevationFormula\s*=\s*'([^']+)'", block)
        if higher_match:
            higher_elev = higher_match.group(1)
        if lower_match:
            lower_elev = lower_match.group(1)
        
        # Extract K value (supports scientific notation)
        k_value = None
        dataset_match = re.search(r"DataSetFormulas\.Strings\s*=\s*\((.*?)\)", block, re.DOTALL)
        if dataset_match:
            formulas_str = dataset_match.group(1)
            # Modified regex to match scientific notation (e.g., 5.37E-5)
            first_line_match = re.search(r"'([\d.]+[Ee]?[+-]?\d*)'", formulas_str)
            if first_line_match:
                k_value = float(first_line_match.group(1))
        
        result[normalized_name] = {
            "polygon": points,
            "layer_constraint": (higher_elev, lower_elev),
            "k_value": k_value,
        }
        
        i = item_end + 1
    
    return result


def load_cellids_from_gpt_objects(
    gwf: flopy.mf6.ModflowGwf,
    screenobjects: Dict[str, Dict[str, Any]]
) -> Dict[str, List[Tuple[int, int, int]]]:
    """Calculate which cells belong to each lithology zone based on ScreenObject polygons from .gpt file"""
    dis = gwf.get_package("DIS")
    mg = StructuredGrid(
        delc=np.array(dis.delc.array),
        delr=np.array(dis.delr.array),
        top=np.array(dis.top.array),
        botm=np.array(dis.botm.array),
        xoff=0.0,
        yoff=0.0,
    )
    
    intersector = GridIntersect(mg, method="structured")
    
    cellids_by_lithology = {}
    
    for lith_name, obj_data in screenobjects.items():
        polygon_coords = obj_data["polygon"]
        if not polygon_coords:
            print(f"⚠️ {lith_name} has no polygon coordinates, skipping")
            continue
        
        poly = Polygon(polygon_coords)
        
        higher_elev, lower_elev = obj_data["layer_constraint"]
        
        if higher_elev == "Model_Top" and lower_elev == "Upper_Aquifer_Bottom":
            target_layer = 0  # Layer 1
        elif higher_elev == "Upper_Aquifer_Bottom" and lower_elev == "Lower_Aquifer_Bottom":
            target_layer = 1  # Layer 2
        else:
            print(f"⚠️ {lith_name} layer_constraint unrecognized: {higher_elev} -> {lower_elev}")
            continue
        
        result = intersector.intersect(poly)
        cellids = []
        
        # GridIntersect returns cellids as (row, col) tuples
        unique_cellids = set()
        for cellid in result["cellids"]:
            row, col = cellid  # cellid is already (row, col) tuple
            unique_cellids.add((target_layer, row, col))
        
        cellids = list(unique_cellids)
        cellids_by_lithology[lith_name] = cellids
        print(f"  - {lith_name}: {len(cellids)} cells (Layer {target_layer + 1})")
    
    return cellids_by_lithology


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


def set_k_for_modmuse_logic(
    sim: MFSimulation, 
    kmap: Dict[str, float],
    cellids_by_lithology: Dict[str, List[Tuple[int, int, int]]]
):
    """Strictly follow ModMuse lithology distribution logic (read from .gpt file)"""
    gwf = sim.get_model(list(sim.model_names)[0])
    npf = gwf.get_package("npf")

    k_arr = np.array(npf.k.array)
    k33_arr = np.array(npf.k33.array)
    
    # Step 1: First set default background values (ModMuse formula)
    if "sandy_gravel" in kmap:
        kx_layer1 = kmap["sandy_gravel"]
        k_arr[0, :, :] = kx_layer1
        k33_arr[0, :, :] = kx_layer1 / 10.0
        print(f"  Set Layer 1 default background sandy_gravel: kx={kx_layer1:.3e}")
    
    if "limestone" in kmap:
        kx_layer2 = kmap["limestone"]
        k_arr[1, :, :] = kx_layer2
        k33_arr[1, :, :] = kx_layer2 / 10.0
        print(f"  Set Layer 2 default background limestone: kx={kx_layer2:.3e}")
    
    # Step 2: Then override specific areas with ScreenObject (K values from kmap parameter)
    for lith_name in LITHOLOGY_OBJECTS:
        if lith_name not in kmap:
            print(f"  ⚠️ {lith_name} not in parameter dict, skipping")
            continue
        
        if lith_name not in cellids_by_lithology:
            print(f"  ⚠️ {lith_name} has no cell info, skipping")
            continue
        
        kx_value = kmap[lith_name]  # Get K value from parameter dict
        kz_value = kx_value / 10.0
        
        cellids = cellids_by_lithology[lith_name]
        
        for lay, row, col in cellids:
            k_arr[lay, row, col] = kx_value
            k33_arr[lay, row, col] = kz_value
        
        print(f"  Override {lith_name}: {len(cellids)} cells, kx={kx_value:.3e}")

    npf.k.set_data(k_arr)
    npf.k33.set_data(k33_arr)
    
    # Don't write here, wait until boundary conditions are set
    # sim.write_simulation()


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
                ghb_data[i]['cond'] = sandy_gravel_kx
                
                if row == 0 and "ghb_north_head" in params:
                    ghb_data[i]['bhead'] = params["ghb_north_head"]
                    layer1_north_count += 1
                elif row == 126 and col <= 91 and "ghb_south_head" in params:
                    ghb_data[i]['bhead'] = params["ghb_south_head"]
                    layer1_south_count += 1
            elif layer == 1:
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
            print(f"  Set Layer 1 North GHB: {layer1_north_count} cells, head={params.get('ghb_north_head', 'N/A'):.3f} m")
        if layer1_south_count > 0:
            print(f"  Set Layer 1 South GHB: {layer1_south_count} cells, head={params.get('ghb_south_head', 'N/A'):.3f} m")
        if layer2_north_count > 0:
            print(f"  Set Layer 2 North GHB: {layer2_north_count} cells, head={params.get('ghb_north_head', 'N/A'):.3f} m")
        if layer2_south_count > 0:
            print(f"  Set Layer 2 South GHB: {layer2_south_count} cells, head={params.get('ghb_south_head', 'N/A'):.3f} m")
    
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
                if not rel_path.startswith("test_run"):
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


def plot_comparison(obs_csv_path: str, test_csv_path: str):
    """Plot comparison: 2 line plots (1 row x 2 columns) - only showing test run results"""
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    
    plt.rcParams['figure.dpi'] = 200
    
    obs_df = pd.read_csv(obs_csv_path, delimiter=CONFIG["observed_csv_delim"])
    if CONFIG["observed_index_col"] in obs_df.columns:
        obs_df = obs_df.set_index(CONFIG["observed_index_col"])
    
    test_df = make_sim_df(test_csv_path)
    
    obs_df = obs_df.loc[test_df.index]
    
    labels = test_df.index.tolist()
    x = np.arange(len(labels))
    
    obs_nat = obs_df['natural'].values
    obs_pmp = obs_df['pumping'].values
    
    test_nat = test_df['natural'].values
    test_pmp = test_df['pumping'].values
    
    # Calculate metrics
    test_resid_nat = test_nat - obs_nat
    test_rmse_nat = np.sqrt(np.mean(test_resid_nat**2))
    test_r2_nat = 1 - (np.sum(test_resid_nat**2) / np.sum((obs_nat - np.mean(obs_nat))**2))
    
    test_resid_pmp = test_pmp - obs_pmp
    test_rmse_pmp = np.sqrt(np.mean(test_resid_pmp**2))
    test_r2_pmp = 1 - (np.sum(test_resid_pmp**2) / np.sum((obs_pmp - np.mean(obs_pmp))**2))
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Subplot 1: Natural
    ax1.plot(x, obs_nat, marker='o', markersize=5, linewidth=1.5, label='Observed', color='#1f77b4')
    ax1.plot(x, test_nat, marker='s', markersize=5, linewidth=1.5, label='Simulated', color='#2ca02c')
    ax1.set_title('Test Run - Natural', fontsize=11, fontweight='bold')
    ax1.set_xlabel('Observation Points', fontsize=10)
    ax1.set_ylabel('Head (m)', fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='best', fontsize=9)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=45, fontsize=8, ha='right')
    ax1.text(0.02, 0.98, f'R² = {test_r2_nat:.3f}\nRMSE = {test_rmse_nat:.3f} m',
             transform=ax1.transAxes, fontsize=9, va='top', ha='left',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.75))
    
    # Subplot 2: Pumping
    ax2.plot(x, obs_pmp, marker='o', markersize=5, linewidth=1.5, label='Observed', color='#1f77b4')
    ax2.plot(x, test_pmp, marker='s', markersize=5, linewidth=1.5, label='Simulated', color='#2ca02c')
    ax2.set_title('Test Run - Pumping', fontsize=11, fontweight='bold')
    ax2.set_xlabel('Observation Points', fontsize=10)
    ax2.set_ylabel('Head (m)', fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best', fontsize=9)
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, rotation=45, fontsize=8, ha='right')
    ax2.text(0.02, 0.98, f'R² = {test_r2_pmp:.3f}\nRMSE = {test_rmse_pmp:.3f} m',
             transform=ax2.transAxes, fontsize=9, va='top', ha='left',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.75))
    
    plt.tight_layout()
    
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
    gpt_path = CONFIG["modmuse_gpt_path"]

    Path(test_run_ws).parent.mkdir(parents=True, exist_ok=True)
    Path(outputs_dir).mkdir(parents=True, exist_ok=True)

    print("\n=== Manual Parameter Test Run (Strict Mode) ===")
    print("Based on flopy-strict.py - reads lithology distribution from .gpt file")
    print("\nCurrent parameter settings:")
    for k, v in MANUAL_PARAMS.items():
        if "head" in k or k == "initial_head":
            print(f"  {k}: {v:.3f}")
        else:
            print(f"  {k}: {v:.3e}")

    # 1. Parse ModMuse .gpt file
    print(f"\n=== Parsing ModMuse .gpt file ===")
    print(f"File path: {gpt_path}")
    screenobjects = parse_gpt_screenobjects(gpt_path)
    print(f"✅ Successfully parsed {len(screenobjects)} lithology objects")

    # 2. Copy model to test directory
    print(f"\nCopying model to: {test_run_ws}")
    copy_model_to_run(base_ws, test_run_ws)

    # 3. Load simulation
    sim = load_simulation(test_run_ws, mf6_exe)
    gwf = sim.get_model(list(sim.model_names)[0])
    
    # 4. Calculate lithology-cell mapping
    print("\n=== Calculating lithology cell distribution ===")
    cellids_by_lithology = load_cellids_from_gpt_objects(gwf, screenobjects)

    # 5. Set parameters
    print("\n=== Setting hydraulic conductivity ===")
    set_k_for_modmuse_logic(sim, MANUAL_PARAMS, cellids_by_lithology)
    
    print("\n=== Setting boundary conditions ===")
    set_boundary_conditions(sim, MANUAL_PARAMS)

    # 6. Run simulation
    print("\n=== Running MODFLOW 6 ===")
    t0 = time.time()
    ok = run_simulation(sim, test_run_ws)
    dt = time.time() - t0

    # 7. Analyze results
    if ok:
        print(f"✅ Run successful! Time elapsed: {dt:.2f} seconds")
        
        # Directly construct CSV file path (don't search to avoid reading old files)
        obs_csv = os.path.join(test_run_ws, "Steady_cali.ob_gw_out_head.csv")
        
        if os.path.exists(obs_csv):
            print(f"\nReading observation data: {obs_csv}")
            print(f"File modification time: {time.ctime(os.path.getmtime(obs_csv))}")
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
                plot_comparison(
                    obs_csv_path=CONFIG["observed_csv_path"],
                    test_csv_path=obs_csv
                )
            else:
                print("\n❌ Model has failed observation points, cannot calculate metrics")
        else:
            print("\n⚠️ Observation data CSV file not found")
    else:
        print(f"❌ Run failed! Time elapsed: {dt:.2f} seconds")

    print(f"\nRun results saved to: {test_run_ws}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error occurred: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
