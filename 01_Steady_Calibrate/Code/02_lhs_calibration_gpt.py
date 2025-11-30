#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MF6 Batch Automation (Strictly Following ModMuse Lithology Distribution Logic - Reading from .gpt File)
========================================================================================================
Key differences from Steady_cali_06.py:
- **Does not use external shapefile**, instead directly parses ScreenObject from ModMuse .gpt file
- **Strictly follows ModMuse internal logic**:
  1. Default formula: If((Layer = 1), sandy_gravel, limestone)
  2. ScreenObject overrides specific regions (Clean_Gravel, Clayey_Sand, Loamy_Sand, Sand_Stone)

This ensures Python/Flopy results are fully consistent with ModMuse GUI results.

Lithology Distribution Logic (Strictly Matching ModMuse):
- Layer 1 background: sandy_gravel
- Layer 2 background: limestone
- Layer 1 specific region overrides: clean_gravel, clayey_sand, loamy_sand (polygons from .gpt file)
- Layer 2 specific region overrides: sandstone (polygons from .gpt file)

Calibration Parameters (10 total):
1. sandy_gravel: Layer 1 background K value
2. limestone: Layer 2 background K value
3. clean_gravel: Layer 1 specific region K value
4. clayey_sand: Layer 1 specific region K value
5. loamy_sand: Layer 1 specific region K value
6. sandstone: Layer 2 specific region K value
7. riv_cond: River boundary conductance
8. initial_head: Initial head
9. ghb_north_head: North boundary head
10. ghb_south_head: South boundary head

Dependencies: numpy, pandas, flopy, scipy, shapely
"""

from __future__ import annotations
import os
import sys
import shutil
import glob
import time
import json
import re
from pathlib import Path
from typing import Dict, List, Tuple, Any

import numpy as np
import pandas as pd

import flopy
from flopy.mf6 import MFSimulation
from flopy.discretization import StructuredGrid
from flopy.utils import GridIntersect
from scipy.optimize import differential_evolution
from shapely.geometry import Polygon

# ------------------------
# User Configuration
# ------------------------
ROOT = Path(__file__).resolve().parent
CONFIG = {
    "base_model_ws": str((ROOT.parent / "Model").resolve()),
    "mf6_exe": r"D:\Program Files\USGS\mf6.6.2_win64\bin\mf6.exe",
    "runs_root": str((ROOT.parent / "Out" / "runs").resolve()),
    "outputs_dir": str((ROOT.parent / "Out").resolve()),
    "n_samples": 5,
    "random_seed": 50,
    "obs_csv_name_hint": "Steady_cali.ob_gw_out_head",
    "observed_csv_path": str((ROOT.parent / "steady_state_calibration.csv").resolve()),
    "observed_csv_delim": ";",
    "observed_index_col": "Unnamed: 0",
    "conditions": ["natural", "pumping"],
    "obs_index_labels": ["P1", "Pz2", "Pz3", "Pz4", "Pz5", "Pz6", "Pz7", "Pz8", "Pz9", "Pz10", "Pz11", "Pz12"],
    # ✅ New: ModMuse .gpt file path for extracting lithology polygons
    "modmuse_gpt_path": str((ROOT.parent / "Model" / "Steady_cali.gpt").resolve()),
}

# ✅ Calibration parameters: 2 default formula parameters + 4 specific region lithology parameters
DEFAULT_PARAMS = [
    "sandy_gravel",   # Layer 1 default (background)
    "limestone",      # Layer 2 default (background)
]

# ✅ Lithology objects defined in ModMuse .gpt file (polygon + K value formula)
LITHOLOGY_OBJECTS = [
    "clean_gravel",   # Layer 1, polygon from .gpt
    "clayey_sand",    # Layer 1, polygon from .gpt
    "loamy_sand",     # Layer 1, polygon from .gpt
    "sandstone",      # Layer 2, polygon from .gpt
]

# ✅ Boundary condition parameters
BOUNDARY_PARAMS = [
    "riv_cond",       # River_Boundary_Whole RIV conductance
    "initial_head",   # Initial head
    "ghb_north_head", # North_Boundary GHB boundary head
    "ghb_south_head", # South_Boundary GHB boundary head
]

# ✅ Parameter ranges (6 lithology parameters + 4 boundary parameters = 10 total)
PARAM_BOUNDS: Dict[str, Tuple[float, float]] = {
    # Default formula parameters (background values)
    "sandy_gravel": (1e-3, 3e-3),      # Layer 1 default
    "limestone": (1e-4, 1e-3),         # Layer 2 default
    # Specific region lithology parameters (override defaults)
    "clean_gravel": (5e-2, 3.5e-1),    # Layer 1 specific region
    "clayey_sand": (1e-5, 5e-4),       # Layer 1 specific region
    "loamy_sand": (1e-5, 1e-3),        # Layer 1 specific region
    "sandstone": (1e-5, 1e-4),         # Layer 2 specific region
    # Boundary condition parameters
    "riv_cond": (8e-5, 2.5e-4),        # River conductance (m2/s)
    "initial_head": (37, 39),          # Initial head (m)
    "ghb_north_head": (37.1, 37.2),    # North boundary head (m)
    "ghb_south_head": (37.5, 38.0),    # South boundary head (m)
}

# ------------------------
# Utility Functions
# ------------------------

def parse_gpt_screenobjects(gpt_path: str) -> Dict[str, Dict[str, Any]]:
    """
    Parse ScreenObject definitions from ModMuse .gpt file
    
    Return format:
    {
        "clean_gravel": {
            "polygon": [(x1, y1), (x2, y2), ...],
            "layer_constraint": ("Model_Top", "Upper_Aquifer_Bottom"),
            "k_value": 0.0001  # or None if not defined
        },
        ...
    }
    """
    result = {}
    
    # Define target lithology object names (case-sensitive)
    target_names = {
        "Clean_Gravel": "clean_gravel",
        "Clayey_Sand": "clayey_sand",
        "Loamy_Sand": "loamy_sand",
        "Sand_Stone": "sandstone",
    }
    
    with open(gpt_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    # Scan lines to find targets
    i = 0
    while i < len(lines):
        line = lines[i]
        
        # Check if line contains target name
        found_name = None
        for target_name, normalized_name in target_names.items():
            if f"ScreenObject.Name = '{target_name}'" in line:
                found_name = (target_name, normalized_name)
                break
        
        if not found_name:
            i += 1
            continue
        
        target_name, normalized_name = found_name
        
        # Search backward to find the start of this item block (look for ClassTypeName)
        item_start = None
        j = i - 1
        while j >= max(0, i - 200):
            if 'item' in lines[j] and lines[j].strip() == 'item':
                # Check if next 15 lines contain ClassTypeName = 'TScreenObject'
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
            print(f"⚠️ {normalized_name}: Cannot find item start")
            i += 1
            continue
        
        # Search forward to find MixtureFormulas.Strings, then find the first 'end' after it
        item_end = None
        mixture_found = False
        for j in range(i, min(i + 500, len(lines))):
            if "MixtureFormulas.Strings" in lines[j]:
                mixture_found = True
            if mixture_found and lines[j].strip() == 'end':
                item_end = j
                break
        
        if item_end is None:
            print(f"⚠️ {normalized_name}: Cannot find item end")
            i += 1
            continue
        
        # Extract this block's content
        block = ''.join(lines[item_start:item_end+1])
        
        # Extract Points (polygon coordinates)
        points = []
        points_match = re.search(r"Points\s*=\s*<(.*?)>", block, re.DOTALL)
        if points_match:
            points_block = points_match.group(1)
            # Each point format: X = ... Y = ...
            point_items = re.finditer(r"X\s*=\s*([\d.e+-]+)\s+Y\s*=\s*([\d.e+-]+)", points_block)
            for point_match in point_items:
                x = float(point_match.group(1))
                y = float(point_match.group(2))
                points.append((x, y))
        
        # Extract Layer constraints
        higher_elev = None
        lower_elev = None
        higher_match = re.search(r"HigherElevationFormula\s*=\s*'([^']+)'", block)
        lower_match = re.search(r"LowerElevationFormula\s*=\s*'([^']+)'", block)
        if higher_match:
            higher_elev = higher_match.group(1)
        if lower_match:
            lower_elev = lower_match.group(1)
        
        # Extract K value (from DataSetFormulas)
        k_value = None
        dataset_match = re.search(r"DataSetFormulas\.Strings\s*=\s*\((.*?)\)", block, re.DOTALL)
        if dataset_match:
            formulas_str = dataset_match.group(1)
            # First line is typically Kx value
            first_line_match = re.search(r"'([\d.e+-]+)'", formulas_str)
            if first_line_match:
                k_value = float(first_line_match.group(1))
        
        result[normalized_name] = {
            "polygon": points,
            "layer_constraint": (higher_elev, lower_elev),
            "k_value": k_value,
        }
        
        # Skip to after item_end and continue
        i = item_end + 1
    
    return result


def load_cellids_from_gpt_objects(
    gwf: flopy.mf6.ModflowGwf,
    screenobjects: Dict[str, Dict[str, Any]]
) -> Dict[str, List[Tuple[int, int, int]]]:
    """
    Calculate which cells belong to each lithology zone based on ScreenObject polygons in .gpt file
    
    Return format:
    {
        "clean_gravel": [(lay, row, col), ...],
        "clayey_sand": [(lay, row, col), ...],
        ...
    }
    """
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
        
        # Determine which layer to apply based on layer_constraint
        higher_elev, lower_elev = obj_data["layer_constraint"]
        
        # Simple check: If Model_Top → Upper_Aquifer_Bottom, then Layer 1
        # If Upper_Aquifer_Bottom → Lower_Aquifer_Bottom, then Layer 2
        if higher_elev == "Model_Top" and lower_elev == "Upper_Aquifer_Bottom":
            target_layer = 0  # Layer 1
        elif higher_elev == "Upper_Aquifer_Bottom" and lower_elev == "Lower_Aquifer_Bottom":
            target_layer = 1  # Layer 2
        else:
            print(f"⚠️ {lith_name} layer_constraint unrecognized: {higher_elev} -> {lower_elev}")
            continue
        
        # Use GridIntersect to find cells intersecting with polygon
        result = intersector.intersect(poly)
        cellids = []
        
        # For structured grid, GridIntersect returns 2D plane intersections
        # cellid is single-layer index (0 to nrow*ncol-1), does not include layer info
        nrow, ncol = mg.nrow, mg.ncol
        
        for item in result["cellids"]:
            cellid, _ = item  # cellid is 2D index: cellid = row * ncol + col
            
            # Convert 2D index to (row, col)
            row = cellid // ncol
            col = cellid % ncol
            
            # Determine which layer to apply based on layer_constraint
            # target_layer has been determined by Higher/LowerElevationFormula
            cellids.append((target_layer, row, col))
        
        cellids_by_lithology[lith_name] = cellids
        print(f"✅ {lith_name}: {len(cellids)} cells (Layer {target_layer + 1})")
    
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


def lhs(n_samples: int, n_dim: int, rng: np.random.Generator) -> np.ndarray:
    """Generate LHS samples in [0,1] interval"""
    cut = np.linspace(0, 1, n_samples + 1)
    u = rng.random((n_samples, n_dim))
    a = cut[:-1][:, None]
    b = cut[1:][:, None]
    pts = a + (b - a) * u
    for j in range(n_dim):
        rng.shuffle(pts[:, j])
    return pts


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
                raise RuntimeError(f"❌ Cannot handle target folder: {e2}") from e2
    
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
    print(f"Loading simulation... sim_ws = {workspace}")
    fix_nam_paths(workspace)
    sim = flopy.mf6.MFSimulation.load(sim_ws=workspace, exe_name=mf6_exe)
    return sim


def set_k_for_modmuse_logic(
    sim: MFSimulation, 
    kmap: Dict[str, float],
    cellids_by_lithology: Dict[str, List[Tuple[int, int, int]]]
):
    """
    ✅ Strictly follow ModMuse lithology distribution logic (read from .gpt file)
    
    ModMuse Logic:
    1. Default formula: If((Layer = 1), sandy_gravel, limestone)
       - Layer 1 background value = sandy_gravel
       - Layer 2 background value = limestone
    
    2. ScreenObject overrides specific regions:
       - Clean_Gravel, Clayey_Sand, Loamy_Sand override specific regions in Layer 1
       - Sand_Stone overrides specific regions in Layer 2
    
    Parameters:
    - kmap: K value dictionary containing all lithology parameters
    - cellids_by_lithology: Cell list parsed from .gpt file
    """
    gwf = sim.get_model(list(sim.model_names)[0])
    npf = gwf.get_package("npf")

    k_arr = np.array(npf.k.array)
    k33_arr = np.array(npf.k33.array)
    
    # Step 1: First set default background values (ModMuse formula)
    if "sandy_gravel" in kmap:
        kx_layer1 = kmap["sandy_gravel"]
        k_arr[0, :, :] = kx_layer1
        k33_arr[0, :, :] = kx_layer1 / 10.0
        print(f"  ✅ Layer 1 default background sandy_gravel: kx={kx_layer1:.3e}")
    
    if "limestone" in kmap:
        kx_layer2 = kmap["limestone"]
        k_arr[1, :, :] = kx_layer2
        k33_arr[1, :, :] = kx_layer2 / 10.0
        print(f"  ✅ Layer 2 default background limestone: kx={kx_layer2:.3e}")
    
    # Step 2: Override specific regions with ScreenObjects
    for lith_name in LITHOLOGY_OBJECTS:
        if lith_name not in kmap:
            continue
        
        if lith_name not in cellids_by_lithology:
            print(f"  ⚠️ {lith_name} has no cell information, skipping")
            continue
        
        kx_value = kmap[lith_name]
        kz_value = kx_value / 10.0
        
        cellids = cellids_by_lithology[lith_name]
        for lay, row, col in cellids:
            k_arr[lay, row, col] = kx_value
            k33_arr[lay, row, col] = kz_value
        
        print(f"  ✅ {lith_name} overrides {len(cellids)} cells: kx={kx_value:.3e}")

    npf.k.set_data(k_arr)
    npf.k33.set_data(k33_arr)
    
    sim.write_simulation()


def set_boundary_conditions(sim: MFSimulation, params: Dict[str, float]):
    """Set boundary condition parameters"""
    gwf = sim.get_model(list(sim.model_names)[0])
    
    # 1. RIV conductance
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
    
    # 2. Initial head
    if "initial_head" in params:
        ic = gwf.get_package("ic")
        if ic is not None:
            strt_arr = np.array(ic.strt.array)
            strt_arr[:, :, :] = params["initial_head"]
            ic.strt.set_data(strt_arr)
            print(f"  Set initial head: {params['initial_head']:.3f} m")
    
    # 3. GHB (calculate conductance using calibrated K values)
    ghb = gwf.get_package("ghb")
    if ghb is not None:
        ghb_data = ghb.stress_period_data.get_data(0)
        sandy_gravel_kx = params.get("sandy_gravel", 0.0005)
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
                # Layer 2: Use limestone kx
                ghb_data[i]['cond'] = limestone_kx
                
                if row == 0 and "ghb_north_head" in params:
                    ghb_data[i]['bhead'] = params["ghb_north_head"]
                    layer2_north_count += 1
                elif row == 126 and col <= 91 and "ghb_south_head" in params:
                    ghb_data[i]['bhead'] = params["ghb_south_head"]
                    layer2_south_count += 1
        
        ghb.stress_period_data.set_data(ghb_data, 0)
        
        if layer1_north_count > 0:
            print(f"  Set Layer 1 North GHB: {layer1_north_count} cells, head={params.get('ghb_north_head', 'N/A'):.3f} m, cond={sandy_gravel_kx:.3e}")
        if layer1_south_count > 0:
            print(f"  Set Layer 1 South GHB: {layer1_south_count} cells, head={params.get('ghb_south_head', 'N/A'):.3f} m, cond={sandy_gravel_kx:.3e}")
        if layer2_north_count > 0:
            print(f"  Set Layer 2 North GHB: {layer2_north_count} cells, head={params.get('ghb_north_head', 'N/A'):.3f} m, cond={limestone_kx:.3e}")
        if layer2_south_count > 0:
            print(f"  Set Layer 2 South GHB: {layer2_south_count} cells, head={params.get('ghb_south_head', 'N/A'):.3f} m, cond={limestone_kx:.3e}")
    
    sim.write_simulation()


def run_simulation(sim: MFSimulation, run_ws: str) -> bool:
    """Run simulation"""
    success, buff = sim.run_simulation(report=True, silent=True)
    if not success:
        print("Run failed, output:\n" + "\n".join(buff))
    return bool(success)


def cleanup_run_folder(run_ws: str, obs_csv_path: str | None = None):
    """Clean up run folder, keep only observation data CSV"""
    if not os.path.exists(run_ws):
        return
    
    obs_csv_name = None
    if obs_csv_path and os.path.exists(obs_csv_path):
        obs_csv_name = os.path.basename(obs_csv_path)
    
    deleted_count = 0
    for item in os.listdir(run_ws):
        item_path = os.path.join(run_ws, item)
        
        if obs_csv_name and item == obs_csv_name:
            continue
        
        if not obs_csv_name and item.endswith('.csv'):
            continue
        
        try:
            if os.path.isfile(item_path):
                os.remove(item_path)
                deleted_count += 1
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)
                deleted_count += 1
        except Exception as e:
            print(f"  ⚠️  Failed to delete {item}: {e}")
    
    if deleted_count > 0:
        print(f"  🧹 Cleanup complete: deleted {deleted_count} file(s)/folder(s)")


def find_obs_csv(run_ws: str, hint: str | List[str] | None = None) -> str | None:
    """Find OBS package exported CSV file"""
    csvs = glob.glob(os.path.join(run_ws, "**", "*.csv"), recursive=True)
    if not csvs:
        return None

    outputs_dir = CONFIG.get("outputs_dir", "")
    observed_csv_path = CONFIG.get("observed_csv_path", "")
    observed_csv_path = os.path.normpath(observed_csv_path) if observed_csv_path else ""

    candidates: List[str] = []
    for p in csvs:
        npth = os.path.normpath(p)
        if outputs_dir:
            abs_p = os.path.abspath(npth)
            abs_out = os.path.abspath(outputs_dir)
            if abs_p.startswith(abs_out + os.sep):
                rel_path = os.path.relpath(abs_p, abs_out)
                if not rel_path.startswith("runs" + os.sep):
                    continue
        if observed_csv_path and os.path.abspath(npth) == os.path.abspath(observed_csv_path):
            continue
        base = os.path.basename(npth).lower()
        if base in {"results_running.csv", "results.csv", "results.xlsx"}:
            continue
        candidates.append(npth)

    if not candidates:
        return None

    def latest(paths: List[str]) -> str:
        return max(paths, key=lambda x: os.path.getmtime(x))

    if hint:
        if isinstance(hint, str):
            hints = [hint]
        else:
            hints = list(hint)
        hints = [h.lower() for h in hints if isinstance(h, str) and h]

        exacts = []
        for p in candidates:
            base = os.path.basename(p).lower()
            for h in hints:
                if base == h or base == f"{h}.csv":
                    exacts.append(p)
        if exacts:
            return latest(exacts)

        matches = [p for p in candidates if any(h in os.path.basename(p).lower() for h in hints)]
        if matches:
            return latest(matches)

    return latest(candidates)


def make_sim_df(csv_path: str) -> pd.DataFrame:
    """Read MF6 OBS CSV, construct DataFrame matching observation data"""
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


# ------------------------
# Main Process
# ------------------------

def main():
    base_ws = CONFIG["base_model_ws"]
    runs_root = CONFIG["runs_root"]
    outputs_dir = CONFIG["outputs_dir"]
    mf6_exe = CONFIG["mf6_exe"]
    n_samples = int(CONFIG["n_samples"]) 
    seed = int(CONFIG["random_seed"]) 
    obs_hint = CONFIG.get("obs_csv_name_hint") or None
    gpt_path = CONFIG["modmuse_gpt_path"]

    Path(runs_root).mkdir(parents=True, exist_ok=True)
    Path(outputs_dir).mkdir(parents=True, exist_ok=True)

    # ✅ Parse lithology objects from ModMuse .gpt file
    print("\n=== Parsing ModMuse .gpt File ===")
    screenobjects = parse_gpt_screenobjects(gpt_path)
    print(f"✅ Successfully parsed {len(screenobjects)} lithology objects:")
    for name, data in screenobjects.items():
        print(f"  - {name}: {len(data['polygon'])} vertices, K={data['k_value']}")
    
    # ✅ Load base model, calculate cell list
    print("\n=== Calculating Lithology Cell Distribution ===")
    sim_temp = load_simulation(base_ws, mf6_exe)
    gwf_temp = sim_temp.get_model(list(sim_temp.model_names)[0])
    cellids_by_lithology = load_cellids_from_gpt_objects(gwf_temp, screenobjects)
    
    rng = np.random.default_rng(seed)

    # ✅ Sample 6 lithology parameters + 4 boundary parameters = 10 parameters
    keys = list(DEFAULT_PARAMS) + list(LITHOLOGY_OBJECTS) + list(BOUNDARY_PARAMS)
    bounds = [PARAM_BOUNDS[k] for k in keys]
    unit = lhs(n_samples, len(keys), rng)
    
    # Handle separately: K values in log space, head values in linear space
    samples = np.empty_like(unit)
    for j, key in enumerate(keys):
        if key in BOUNDARY_PARAMS and "head" in key:
            lo, hi = bounds[j]
            samples[:, j] = lo + unit[:, j] * (hi - lo)
        else:
            lo, hi = bounds[j]
            if lo <= 0 or hi <= 0:
                raise ValueError(f"Log-space sampling requires boundary values > 0, but {key} = ({lo}, {hi})")
            log_lo = np.log10(lo)
            log_hi = np.log10(hi)
            samples[:, j] = 10 ** (log_lo + unit[:, j] * (log_hi - log_lo))
    
    # Round to three significant figures
    def round_to_n_sig_figs(x, n=3):
        if x == 0:
            return 0
        return round(x, -int(np.floor(np.log10(abs(x)))) + (n - 1))
    
    for i in range(samples.shape[0]):
        for j in range(samples.shape[1]):
            samples[i, j] = round_to_n_sig_figs(samples[i, j], 3)
    
    print("\n=== Sampling Statistics (Strict ModMuse Logic) ===")
    print("✅ This version reads lithology distribution from .gpt file, fully matching ModMuse internal logic")
    print("   - Default background: Layer 1 → sandy_gravel, Layer 2 → limestone")
    print("   - ScreenObject overrides: clean_gravel, clayey_sand, loamy_sand, sandstone\n")
    for j, k in enumerate(keys):
        print(f"{k}:")
        print(f"  Range: [{PARAM_BOUNDS[k][0]:.3e}, {PARAM_BOUNDS[k][1]:.3e}]")
        print(f"  Samples: min={samples[:, j].min():.3e}, max={samples[:, j].max():.3e}, median={np.median(samples[:, j]):.3e}")
    print()

    records = []

    for idx in range(n_samples):
        run_id = f"run_{idx+1:04d}"
        run_ws = os.path.join(runs_root, run_id)
        print(f"\n=== Starting {run_id} (ModMuse Formula Mode) ===")

        kmap = {k: float(samples[idx, j]) for j, k in enumerate(keys)}
        print("Parameters:")
        for k, v in kmap.items():
            if "head" in k or k == "initial_head":
                print(f"  {k}: {v:.3f}")
            else:
                print(f"  {k}: {v:.3e}")

        copy_model_to_run(base_ws, run_ws)

        sim = load_simulation(run_ws, mf6_exe)
        # ✅ Use strict ModMuse logic (read lithology distribution from .gpt file)
        set_k_for_modmuse_logic(sim, kmap, cellids_by_lithology)
        set_boundary_conditions(sim, kmap)

        t0 = time.time()
        ok = run_simulation(sim, run_ws)
        dt = time.time() - t0

        rmse, r2 = float("nan"), float("nan")
        sim_obs_rows = 0
        obs_csv = None
        skip_this_run = False
        
        if ok:
            obs_csv = find_obs_csv(run_ws, hint=obs_hint)
            if obs_csv and os.path.exists(obs_csv):
                sim_df = make_sim_df(obs_csv)
                sim_obs_rows = len(sim_df)
                
                for condition in CONFIG["conditions"]:
                    if condition in sim_df.columns:
                        simulated = sim_df[condition].values.astype(float)
                        if (simulated < -1e+20).any():
                            print(f"⚠️ Detected failed observation points in {condition} condition, skipping this run")
                            skip_this_run = True
                            break
                
                if not skip_this_run:
                    metrics = compare_to_observations(sim_df)
                    rmse_nat, r2_nat = metrics.get("natural", (float("nan"), float("nan")))
                    rmse_pmp, r2_pmp = metrics.get("pumping", (float("nan"), float("nan")))
                else:
                    rmse_nat = r2_nat = rmse_pmp = r2_pmp = float("nan")
            else:
                print("OBS CSV file not found")
                skip_this_run = True
                rmse_nat = r2_nat = rmse_pmp = r2_pmp = float("nan")
        else:
            print("Simulation failed")
            skip_this_run = True
            rmse_nat = r2_nat = rmse_pmp = r2_pmp = float("nan")
        
        if skip_this_run:
            print(f"⏭️ Skipping {run_id}")
            continue
        
        cleanup_run_folder(run_ws, obs_csv)

        rec = {
            "run_id": idx + 1,
            "ok": ok,
            "runtime_s": round(dt, 3),
            "obs_rows": sim_obs_rows,
            "RMSE[natural]": rmse_nat,
            "R2[natural]": r2_nat,
            "RMSE[pumping]": rmse_pmp,
            "R2[pumping]": r2_pmp,
        }
        rec.update({f"kx[{k}]": v for k, v in kmap.items()})
        records.append(rec)

        pd.DataFrame(records).to_csv(os.path.join(outputs_dir, "results_running.csv"), index=False)

    results_df = pd.DataFrame(records)
    xlsx_path = os.path.join(outputs_dir, "results.xlsx")
    with pd.ExcelWriter(xlsx_path, engine="xlsxwriter") as writer:
        results_df.to_excel(writer, index=False, sheet_name="runs")
    print(f"\n✅ All completed (ModMuse Formula Mode). Results written to: {xlsx_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("Error occurred:", e)
        sys.exit(1)
