#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MF6 Batch Automation (Modify kx by Lithology) - Minimal Framework
================================================
Lithology distribution logic (strictly matching ModMuse):
- Layer 1 background: sandy_gravel
- Layer 2 background: limestone
- Layer 1 specific zone override: clean_gravel, clayey_sand, loamy_sand (polygon read from .gpt file)
- Layer 2 specific zone override: sandstone (polygon read from .gpt file)

Calibration parameters (10 total):
1. sandy_gravel: Layer 1 background K value
2. limestone: Layer 2 background K value
3. clean_gravel: Layer 1 specific zone K value
4. clayey_sand: Layer 1 specific zone K value
5. loamy_sand: Layer 1 specific zone K value
6. sandstone: Layer 2 specific zone K value
7. riv_cond: River boundary conductance
8. initial_head: Initial head
9. ghb_north_head: North boundary head
10. ghb_south_head: South boundary head

Dependencies: numpy, pandas, flopy, openpyxl/xlsxwriter
Note: No additional libraries needed (includes a lightweight LHS implementation).
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

import numpy as np
import pandas as pd

import flopy
from flopy.mf6 import MFSimulation

# ------------------------
# User Configuration
# ------------------------
# All relative paths are based on the script directory to avoid file-not-found errors when running from different working directories
ROOT = Path(__file__).resolve().parent
CONFIG = {
    # ModMuse exported MF6 base model directory (must contain mfsim.nam and other input files)
    # Assumes it's at the same level as Code/: Steady_cali/ directory contains MF6 input.
    "base_model_ws": str((ROOT.parent / "Model").resolve()),
    # MF6 executable path (can be empty if already in system PATH):
    "mf6_exe": r"D:\Program Files\USGS\mf6.6.2_win64\bin\mf6.exe",
    # Batch run root directory (Out/runs at same level as Code)
    "runs_root": str((ROOT.parent / "Out" / "runs").resolve()),
    # Summary results output (Out/ at same level as Code)
    "outputs_dir": str((ROOT.parent / "Out").resolve()),
    # Sampling settings:
    "n_samples":3000,
    "random_seed": 50,
    # Observation output: OBS package CSV filename hint (for auto-finding in run directory)
    "obs_csv_name_hint": "Steady_cali.ob_gw_out_head",  # ✅ Corrected to actual filename (dot not underscore)

    # ==== Observation comparison settings: adjust according to your observation file ====
    "observed_csv_path": str((ROOT.parent / "steady_state_calibration.csv").resolve()),  # ✅ Corrected path (at project root level)
    "observed_csv_delim": ";",
    "observed_index_col": "Unnamed: 0",
    "conditions": ["natural", "pumping"],
    "obs_index_labels": ["P1", "Pz2", "Pz3", "Pz4", "Pz5", "Pz6", "Pz7", "Pz8", "Pz9", "Pz10", "Pz11", "Pz12"],

    # ==== Lithology shapefile auto-intersection configuration ====
    "lithology_shapes": {
        "path": str((ROOT.parent / "Lithology" / "KX.shp").resolve()),
        "name_field": "OBJECTNAME",          # ✅ Corrected to actual field name
        # Map shapefile attribute values to script internal keys (consistent with LITHOLOGY_TARGETS, recommend lowercase + underscore)
        "value_to_key": {
            "Clean_Gravel": "clean_gravel",   # ✅ Note: case must match shapefile
            "Clayey_Sand": "clayey_sand",
            "Loamy_Sand": "loamy_sand",
            "Sand_Stone": "sandstone"
        },
        # Layer number for each lithology (0-indexed): Layer1=0, Layer2=1
        # Determined by HIGH_Z/LOW_Z fields:
        # - Clayey_Sand, Clean_Gravel, Loamy_Sand: Model_Top -> Upper_Aquifer_Bottom (Layer 0)
        # - Sand_Stone: Upper_Aquifer_Bottom -> Lower_Aquifer_Bottom (Layer 1)
        "layer_by_key": {
            "clean_gravel": 0,
            "clayey_sand": 0,
            "loamy_sand": 0,
            "sandstone": 1
        }
    },
}

from pathlib import Path
base = Path(CONFIG["base_model_ws"])
print("Check mfsim.nam:", (base / "mfsim.nam").exists())
# Find gwf name file
for line in (base / "mfsim.nam").read_text(encoding="utf-8", errors="ignore").splitlines():
    t = line.split()
    if len(t) >= 3 and t[0].upper().startswith("GWF6"):
        gwf_nam = (base / t[2]).resolve()
        print("gwf nam path:", gwf_nam)
        print("Exists?", gwf_nam.exists())
        # List package files in this directory
        print("Nearby package files:", list(gwf_nam.parent.glob("Steady.*")))


# Lithology names to calibrate (keys consistent with the range dictionary below)
LITHOLOGY_TARGETS = [
    "clean_gravel",   # Layer 1 polygon/project
    "clayey_sand",    # Layer 1 polygon/project
    "loamy_sand",     # Layer 1 polygon/project
    "sandstone",      # Layer 2 polygon/project
]

# ✅ New: Two parameters from default formula (not from shapefile, set directly by layer)
DEFAULT_PARAMS = [
    "sandy_gravel",   # Layer 1 default value (corresponds to 0.0005 in original formula)
    "limestone",      # Layer 2 default value (corresponds to 3E-5 in original formula)
]

# ✅ New: Boundary condition parameters (RIV, GHB, initial head)
BOUNDARY_PARAMS = [
    "riv_cond",       # River_Boundary_Whole RIV conductance
    "initial_head",   # Initial head
    "ghb_north_head", # North_Boundary GHB boundary head
    "ghb_south_head", # South_Boundary GHB boundary head
]

# kx value range for each lithology (m/s or your unit; consistent with model units).
# Adjust according to your experience; current values are examples.
# clean_gravel  >  sandy_gravel  >  loamy_sand  >  clayey_sand  >  limestone  >  sandstone
PARAM_BOUNDS: Dict[str, Tuple[float, float]] = {
    "clean_gravel": (5e-3, 1e-1),    
    "clayey_sand": (1e-5, 5e-4),
    "loamy_sand": (1e-5, 1e-3),
    "sandstone": (1e-5, 1e-4),
    # ✅ New: Two default expression parameters
    "sandy_gravel": (1e-3, 3e-3),   # Layer 1 default (verified: >3.5e-4 works better)
    "limestone": (1e-4, 1e-3),        # Layer 2 default
    # ✅ New: Boundary condition parameters
    "riv_cond": (8e-5, 2.5e-4),       # River_Boundary_Whole RIV conductance (m2/s)
    "initial_head": (37, 39),     # ✅ Corrected: narrowed range to avoid conflict with boundary conditions
    "ghb_north_head": (37.1, 37.2),   # ✅ Corrected: expanded lower bound to avoid conflict with initial_head (was 36.8→36.0)
    "ghb_south_head": (37.5, 38.0),   # South_Boundary GHB boundary head (m)
}

# ------------------------
# Functions to be adapted/implemented
# ------------------------

def load_cellids_by_lithology(sim: MFSimulation) -> Dict[str, List[Tuple[int, int, int]]]:
    """Generate (lay,row,col) cellid list for each lithology based on *Shapefile polygon* intersection.

    Requirements:
    - Provide shapefile path, attribute field name, attribute value to script internal key mapping, and layer number for each lithology in CONFIG.
    - Only implements **structured grid (DIS)** (lay,row,col) output; for **DISV**, slight modification needed (assign by (lay, icell2d)).
    """
    from flopy.utils.gridintersect import GridIntersect
    try:
        import fiona  # Lightweight shapefile reader
        from shapely.geometry import shape as shp_shape
    except Exception as e:
        raise RuntimeError("fiona and shapely are required: pip install fiona shapely") from e

    gwf = sim.get_model(list(sim.model_names)[0])
    mg = gwf.modelgrid

    if mg.grid_type.lower() != "structured":
        raise NotImplementedError(
            "Current implementation assumes DIS structured grid (nlay,nrow,ncol). For DISV/unsaturated grid, please notify for modification to (lay,icell2d)."
        )

    shp_cfg = CONFIG.get("lithology_shapes", {})
    shp_path = shp_cfg.get("path")
    name_field = shp_cfg.get("name_field", "NAME")
    value_to_key: Dict[str, str] = shp_cfg.get("value_to_key", {})
    layer_by_key: Dict[str, int] = shp_cfg.get("layer_by_key", {})

    if not shp_path or not os.path.exists(shp_path):
        raise FileNotFoundError(f"Shapefile not found: {shp_path}")

    gi = GridIntersect(mg, method="vertex")

    # Result container
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
            # Map attribute value to internal key, e.g., "clean gravel" -> "clean_gravel"
            key = value_to_key.get(str(raw_val))
            if key is None:
                # If no mapping, but attribute value matches key (already underscored), allow it
                key = str(raw_val)
            if key not in LITHOLOGY_TARGETS:
                # Not relevant to this calibration, skip
                continue
            lay = layer_by_key.get(key)
            if lay is None:
                raise ValueError(f"Missing layer mapping for lithology {key} in layer_by_key[{key}].")

            # Intersect with 2D grid to get (row,col)
            res = gi.intersect(poly)
            # GridIntersect returns record array; in structured grid usually has cellids as (i,j)
            # Defensive handling:
            cell_rcs = []
            try:
                cell_rcs = [tuple(rc) for rc in res['cellids']]
            except Exception:
                # Compatibility with older versions: res is list[CellInfo]
                cell_rcs = [tuple(r.cellid) for r in res]

            # Extend to (lay,row,col)
            out[key].extend([(lay, r, c) for (r, c) in cell_rcs])

    # Print statistics for verification
    for k in LITHOLOGY_TARGETS:
        print(f"  - {k}: {len(out.get(k, []))} cells")

    return out


def compare_to_observations(sim_df: pd.DataFrame) -> Dict[str, Tuple[float, float]]:
    """Calculate RMSE and R² between simulated and observed values
    
    R² formula: R² = 1 - SS_res / SS_tot
    where:
    - SS_res = Σ(observed - simulated)²  (sum of squared residuals)
    - SS_tot = Σ(observed - mean(observed))²  (total sum of squares)
    
    Notes:
    - R² < 0 means model is worse than using mean for prediction
    - R² = 0 means model equals using mean for prediction
    - R² = 1 means perfect fit
    - 🆕 Filters out MODFLOW failed points (value = -1e+30)
    """
    obs_path = CONFIG["observed_csv_path"]
    delim = CONFIG["observed_csv_delim"]
    index_col = CONFIG["observed_index_col"]
    conds = CONFIG["conditions"]

    # Read observations
    obs_df = pd.read_csv(obs_path, delimiter=delim)
    if index_col in obs_df.columns:
        obs_df = obs_df.set_index(index_col)
    
    # Align order (ensure observation and simulation data observation points are in same order)
    obs_df = obs_df.loc[sim_df.index]

    metrics: Dict[str, Tuple[float, float]] = {}
    for condition in conds:
        observed = obs_df[condition].values.astype(float)
        simulated = sim_df[condition].values.astype(float)
        
        # 🆕 Filter out MODFLOW failed points (marked as -1e+30)
        valid_mask = (simulated > -1e+20)  # Filter out extreme negative values
        
        if valid_mask.sum() == 0:
            # All points failed
            metrics[condition] = (float("inf"), float("nan"))
            continue
        
        observed = observed[valid_mask]
        simulated = simulated[valid_mask]
        
        # RMSE: Root mean square error of residuals
        residuals = simulated - observed
        rmse = float(np.sqrt(np.mean(residuals ** 2)))
        
        # R²: Coefficient of determination
        ss_res = float(np.sum(residuals ** 2))  # Residual sum of squares
        ss_tot = float(np.sum((observed - np.mean(observed)) ** 2))  # Total sum of squares
        
        if ss_tot > 1e-10:  # Avoid division by 0
            r2 = float(1.0 - ss_res / ss_tot)
        else:
            r2 = float("nan")
        
        metrics[condition] = (rmse, r2)
    
    return metrics

# ------------------------
# Sampling Tools: Simple LHS (Latin Hypercube Sampling)
# ------------------------

def lhs(n_samples: int, n_dim: int, rng: np.random.Generator) -> np.ndarray:
    """Generate LHS samples in [0,1] interval (n_samples, n_dim).
    Reference: McKay et al. 1979, simplified implementation to avoid external dependencies.
    Fix: Ensure segment endpoints `a`/`b` as (n_samples,1) can broadcast with (n_samples,n_dim).
    """
    cut = np.linspace(0, 1, n_samples + 1)
    u = rng.random((n_samples, n_dim))
    a = cut[:-1][:, None]  # (n_samples,1)
    b = cut[1:][:, None]   # (n_samples,1)
    pts = a + (b - a) * u  # (n_samples,n_dim)
    for j in range(n_dim):
        rng.shuffle(pts[:, j])
    return pts


def scale_unit_lhs(unit_samples: np.ndarray, bounds: List[Tuple[float, float]], log_space: bool = False) -> np.ndarray:
    """Scale [0,1] samples to real parameter space according to bounds, limited to 3 significant figures.
    
    Parameters:
    -----------
    unit_samples : LHS samples in [0,1] interval
    bounds : (lower, upper) bounds for each parameter
    log_space : Whether to sample in log space (default False)
    
    Returns:
    --------
    scaled : Scaled parameter samples (limited to 3 significant figures)
    
    Notes:
    ------
    When log_space=True, first takes log10 of bounds, does LHS in log space, then 10** back to original space.
    This avoids sample compression at high values, better matching the multi-order-of-magnitude distribution of hydrological parameters.
    Finally limits all values to 3 significant figures.
    """
    unit_samples = np.asarray(unit_samples)
    assert unit_samples.shape[1] == len(bounds)
    scaled = np.empty_like(unit_samples)
    
    if log_space:
        # Sample in log space
        for j, (lo, hi) in enumerate(bounds):
            if lo <= 0 or hi <= 0:
                raise ValueError(f"Log-space sampling requires boundary values > 0, but bounds[{j}] = ({lo}, {hi})")
            log_lo = np.log10(lo)
            log_hi = np.log10(hi)
            # Linear interpolation in log space
            scaled[:, j] = 10 ** (log_lo + unit_samples[:, j] * (log_hi - log_lo))
    else:
        # Linear space sampling
        for j, (lo, hi) in enumerate(bounds):
            scaled[:, j] = lo + unit_samples[:, j] * (hi - lo)
    
    # ✅ Round to three significant figures
    def round_to_n_sig_figs(x, n=3):
        """Limit value to n significant figures"""
        if x == 0:
            return 0
        return round(x, -int(np.floor(np.log10(abs(x)))) + (n - 1))
    
    # Apply 3 significant figure limit to each element
    for i in range(scaled.shape[0]):
        for j in range(scaled.shape[1]):
            scaled[i, j] = round_to_n_sig_figs(scaled[i, j], 3)
    
    return scaled

# ------------------------
# MF6 Related Tools
# ------------------------

def copy_model_to_run(base_ws: str, run_ws: str):
    """Enhanced version: Handle OneDrive paths and permission issues"""
    import stat
    
    # If target exists, force delete (handle read-only files)
    if os.path.exists(run_ws):
        def remove_readonly(func, path, _):
            """Clear read-only attribute before deletion"""
            try:
                os.chmod(path, stat.S_IWRITE)
                func(path)
            except Exception:
                pass
        
        try:
            shutil.rmtree(run_ws, onerror=remove_readonly)
        except Exception as e:
            print(f"⚠️ Failed to delete old folder: {e}")
            # Try renaming old folder
            backup = run_ws + "_old_" + str(int(time.time()))
            try:
                os.rename(run_ws, backup)
                print(f"Renamed old folder to: {backup}")
            except Exception as e2:
                raise RuntimeError(f"❌ Cannot process target folder: {e2}\nSuggestions:\n1. Pause OneDrive sync\n2. Move project to local disk (e.g., C:\\Temp)\n3. Manually delete {run_ws}") from e2
    
    # Wait for filesystem sync (OneDrive may delay)
    time.sleep(0.3)
    
    # Copy
    try:
        shutil.copytree(base_ws, run_ws)
    except Exception as e:
        raise RuntimeError(f"❌ Failed to copy model: {e}") from e


def fix_nam_paths(workspace: str):
    """Fix relative paths in .nam files to point to current directory"""
    import re
    
    # Fix mfsim.nam
    mfsim_path = os.path.join(workspace, "mfsim.nam")
    if os.path.exists(mfsim_path):
        with open(mfsim_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        # Change relative paths to current directory
        content = re.sub(r'\.\.[\\/]\.\.[\\/]model[\\/]model\.\w+[\\/]', '', content)
        content = re.sub(r'\.\.[\\/]\.\.[\\/]output[\\/]output\.\w+[\\/]', '', content)
        with open(mfsim_path, 'w', encoding='utf-8') as f:
            f.write(content)
    
    # Fix Steady.nam (or other gwf name files)
    for nam_file in glob.glob(os.path.join(workspace, "*.nam")):
        if os.path.basename(nam_file).lower() == "mfsim.nam":
            continue
        with open(nam_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        # Change relative paths to current directory
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


def set_k_for_zones(sim: MFSimulation, kmap: Dict[str, float], cellids: Dict[str, List[Tuple[int,int,int]]]):
    """
    Update kx by lithology (also ky=kx, kz=kx/10).
    
    ✅ New logic:
    1. First set global default values (by layer):
       - Layer 1 (lay=0): sandy_gravel kx
       - Layer 2 (lay=1): limestone kx
    2. Then override specific cells with shapefile-defined lithology values
    
    This simulates the ModMuse default formula: If((Layer = 1), 0.0005, 3E-5)

    - For MF6 NPF: `k` is used for horizontal hydraulic conductivity (x/y), `k33` for vertical.
    - If your model enables anisotropy (k22/angle etc.), please extend as needed.
    """
    # Get first model (assuming only one gwf)
    gwf = sim.get_model(list(sim.model_names)[0])
    npf = gwf.get_package("npf")

    # Read existing k, k33 as numpy arrays
    k_arr = np.array(npf.k.array)      # shape: (nlay, nrow, ncol)
    k33_arr = np.array(npf.k33.array)  # same
    
    # ✅ Step 1: First set global default values (by layer)
    # Layer 1 (lay=0) → sandy_gravel
    if "sandy_gravel" in kmap:
        kx_layer1 = kmap["sandy_gravel"]
        k_arr[0, :, :] = kx_layer1
        k33_arr[0, :, :] = kx_layer1 / 10.0
        print(f"  Set Layer 1 default (sandy_gravel): kx={kx_layer1:.3e}")
    
    # Layer 2 (lay=1) → limestone
    if "limestone" in kmap:
        kx_layer2 = kmap["limestone"]
        k_arr[1, :, :] = kx_layer2
        k33_arr[1, :, :] = kx_layer2 / 10.0
        print(f"  Set Layer 2 default (limestone): kx={kx_layer2:.3e}")

    # ✅ Step 2: Override specific cells with shapefile-defined lithology values
    for lith, kx in kmap.items():
        if lith in DEFAULT_PARAMS or lith in BOUNDARY_PARAMS:
            # Skip default params and boundary params (handled elsewhere)
            continue
        cells = cellids.get(lith, [])
        if not cells:
            continue
        print(f"  Override {lith}: {len(cells)} cells, kx={kx:.3e}")
        for lay, i, j in cells:
            k_arr[lay, i, j] = kx             # kx=ky
            k33_arr[lay, i, j] = kx / 10.0    # kz = kx/10

    # Write back
    npf.k.set_data(k_arr)
    npf.k33.set_data(k33_arr)
    
    # ✅ Key: Write to file!
    sim.write_simulation()


def set_boundary_conditions(sim: MFSimulation, params: Dict[str, float]):
    """
    Set boundary condition parameters.
    
    Parameters:
    -----------
    sim : MODFLOW 6 simulation object
    params : Parameter dictionary containing:
        - riv_cond: River_Boundary_Whole RIV conductance
        - initial_head: Initial head
        - ghb_north_head: North_Boundary GHB boundary head (row=0, i.e. row 1)
        - ghb_south_head: South_Boundary GHB boundary head (row=126, i.e. row 127, col<=91)
        - sandy_gravel: Used to set GHB conductance
    
    Boundary locations:
    ------------------
    - North boundary: row=0 (row 1)
    - South boundary: row=126 (row 127), col<=91 (col 92 and above is river)
    """
    gwf = sim.get_model(list(sim.model_names)[0])
    
    # 1. Set RIV conductance (River_Boundary_Whole)
    if "riv_cond" in params:
        riv = gwf.get_package("riv")
        if riv is not None:
            riv_data = riv.stress_period_data.get_data(0)  # numpy.recarray
            # ✅ Fix: Use field name access, not index
            # Only modify cells with boundname 'river_boundary_whole'
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
            # Get current head array
            strt_arr = np.array(ic.strt.array)
            # Set initial head for all layers
            strt_arr[:, :, :] = params["initial_head"]
            ic.strt.set_data(strt_arr)
            print(f"  Set initial head: {params['initial_head']:.3f} m")
    
    # 3. Set GHB (North_Boundary and South_Boundary)
    ghb = gwf.get_package("ghb")
    if ghb is not None:
        ghb_data = ghb.stress_period_data.get_data(0)  # numpy.recarray
        sandy_gravel_kx = params.get("sandy_gravel", 0.0005)  # Layer 1 default
        sandstone_kx = params.get("sandstone", 1e-5)  # Layer 2 South boundary
        limestone_kx = params.get("limestone", 3e-5)  # Layer 2 North boundary
        
        # ✅ Fix: Use field name access, cellid is a tuple (layer, row, col)
        # GHB data format: cellid=(layer, row, col), bhead, cond, iface, boundname
        # North boundary: row=0 (row 1)
        # South boundary: row=126 (row 127), col<=91
        # ✅ Layer 1: Use sandy_gravel kx
        # ✅ Layer 2: South boundary uses sandstone kx, North boundary uses limestone kx
        
        layer1_north_count = 0
        layer1_south_count = 0
        layer2_north_count = 0
        layer2_south_count = 0
        
        for i in range(len(ghb_data)):
            cellid = ghb_data[i]['cellid']
            layer, row, col = cellid[0], cellid[1], cellid[2]
            
            if layer == 0:
                # Layer 1 (layer=0): Use sandy_gravel kx
                ghb_data[i]['cond'] = sandy_gravel_kx
                
                # Set boundary head based on exact position
                if row == 0 and "ghb_north_head" in params:
                    # North boundary (row 1, 0-indexed = 0)
                    ghb_data[i]['bhead'] = params["ghb_north_head"]
                    layer1_north_count += 1
                elif row == 126 and col <= 91 and "ghb_south_head" in params:
                    # South boundary (row 127, 0-indexed = 126, col <= 92, 0-indexed <= 91)
                    ghb_data[i]['bhead'] = params["ghb_south_head"]
                    layer1_south_count += 1
            elif layer == 1:
                # Layer 2 (layer=1): South boundary uses sandstone, North boundary uses limestone
                if row == 0:
                    # North boundary: Use limestone kx
                    ghb_data[i]['cond'] = limestone_kx
                    if "ghb_north_head" in params:
                        ghb_data[i]['bhead'] = params["ghb_north_head"]
                    layer2_north_count += 1
                elif row == 126 and col <= 91:
                    # South boundary: Use sandstone kx
                    ghb_data[i]['cond'] = sandstone_kx
                    if "ghb_south_head" in params:
                        ghb_data[i]['bhead'] = params["ghb_south_head"]
                    layer2_south_count += 1
        
        ghb.stress_period_data.set_data(ghb_data, 0)
        
        # Print setting information
        if layer1_north_count > 0:
            print(f"  Set Layer 1 North GHB (row=0): {layer1_north_count} cells, head={params.get('ghb_north_head', 'N/A'):.3f} m, cond={sandy_gravel_kx:.3e}")
        if layer1_south_count > 0:
            print(f"  Set Layer 1 South GHB (row=126, col<=91): {layer1_south_count} cells, head={params.get('ghb_south_head', 'N/A'):.3f} m, cond={sandy_gravel_kx:.3e}")
        if layer2_north_count > 0:
            print(f"  Set Layer 2 North GHB (row=0): {layer2_north_count} cells, head={params.get('ghb_north_head', 'N/A'):.3f} m, cond={limestone_kx:.3e}")
        if layer2_south_count > 0:
            print(f"  Set Layer 2 South GHB (row=126, col<=91): {layer2_south_count} cells, head={params.get('ghb_south_head', 'N/A'):.3f} m, cond={sandstone_kx:.3e}")
    
    # Write to file
    sim.write_simulation()


def relax_solver_settings(sim: MFSimulation):
    """
    Relax solver settings to improve convergence success rate
    
    For parameter calibration, we prioritize model running over high precision
    """
    ims = sim.get_solution_groups()[0][1]  # Get IMS solver
    
    # Relax convergence criteria
    ims.outer_maximum = 50             # Increase outer iterations (was 25→50)
    ims.outer_dvclose = 1e-2           # Relax outer convergence (was 1e-3→1e-2)
    ims.inner_maximum = 100            # Increase inner iterations (was 50→100)
    ims.inner_dvclose = 1e-3           # Relax inner convergence (was 1e-4→1e-3)
    ims.rcloserecord = [0.1, "strict"]  # Keep residual criterion


def run_simulation(sim: MFSimulation, run_ws: str) -> bool:
    success, buff = sim.run_simulation(report=True, silent=True)
    if not success:
        # Print output to help debugging
        print("Run failed, output:\n" + "\n".join(buff))
    return bool(success)


def cleanup_run_folder(run_ws: str, obs_csv_path: str | None = None):
    """
    Clean up run folder, keep only observation data CSV file.
    
    Parameters:
        run_ws: Run folder path
        obs_csv_path: Full path to observation data CSV (if known)
    """
    if not os.path.exists(run_ws):
        return
    
    # If obs_csv_path not specified, try to find it
    obs_csv_name = None
    if obs_csv_path and os.path.exists(obs_csv_path):
        obs_csv_name = os.path.basename(obs_csv_path)
    
    # Iterate through all files in the folder
    deleted_count = 0
    for item in os.listdir(run_ws):
        item_path = os.path.join(run_ws, item)
        
        # If it's the observation data CSV, keep it
        if obs_csv_name and item == obs_csv_name:
            continue
        
        # If no explicit obs_csv_name, keep all .csv files (just in case)
        if not obs_csv_name and item.endswith('.csv'):
            continue
        
        # Delete all other files and folders
        try:
            if os.path.isfile(item_path):
                os.remove(item_path)
                deleted_count += 1
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)
                deleted_count += 1
        except Exception as e:
            print(f"  ⚠️  Cannot delete {item}: {e}")
    
    if deleted_count > 0:
        print(f"  🧹 Cleanup complete: deleted {deleted_count} files/folders, kept observation data CSV")


def find_obs_csv(run_ws: str, hint: str | List[str] | None = None) -> str | None:
    """Robustly find the OBS package exported CSV file in *single run directory*.
    Strategy:
    1) Only recursively search for .csv in run_ws; exclude outputs_dir (summary results), exclude observation comparison CSV itself.
    2) If hint provided (can be string or list of strings), first try **exact filename match** (with or without .csv), otherwise fallback to **contains keyword** match;
       When multiple hits, take the file with **latest modification time**.
    3) If hint doesn't match, select the **latest modified** CSV from candidates.
    """
    # 1) All CSVs
    csvs = glob.glob(os.path.join(run_ws, "**", "*.csv"), recursive=True)
    if not csvs:
        return None

    # Exclude outputs_dir and observation comparison CSV
    outputs_dir = CONFIG.get("outputs_dir", "")
    observed_csv_path = CONFIG.get("observed_csv_path", "")
    observed_csv_path = os.path.normpath(observed_csv_path) if observed_csv_path else ""

    candidates: List[str] = []
    for p in csvs:
        npth = os.path.normpath(p)
        # Fix: Only exclude files directly under outputs_dir, allow files under outputs_dir/runs/...
        # Don't exclude files under run directory (run_ws)
        if outputs_dir:
            abs_p = os.path.abspath(npth)
            abs_out = os.path.abspath(outputs_dir)
            # Exclude files directly in outputs_dir (not in runs subfolder)
            # If file is under outputs_dir/runs/..., don't exclude
            if abs_p.startswith(abs_out + os.sep):
                # Check if it's in runs subdirectory
                rel_path = os.path.relpath(abs_p, abs_out)
                if not rel_path.startswith("runs" + os.sep):
                    continue
        if observed_csv_path and os.path.abspath(npth) == os.path.abspath(observed_csv_path):
            continue  # Exclude observation comparison file
        base = os.path.basename(npth).lower()
        if base in {"results_running.csv", "results.csv", "results.xlsx"}:
            continue
        candidates.append(npth)

    if not candidates:
        return None

    def latest(paths: List[str]) -> str:
        return max(paths, key=lambda x: os.path.getmtime(x))

    # 2) Handle hint matching (priority: exact filename, then contains keyword)
    if hint:
        if isinstance(hint, str):
            hints = [hint]
        else:
            hints = list(hint)
        hints = [h.lower() for h in hints if isinstance(h, str) and h]

        # 2a) Exact filename match (with or without .csv)
        exacts = []
        for p in candidates:
            base = os.path.basename(p).lower()
            for h in hints:
                if base == h or base == f"{h}.csv":
                    exacts.append(p)
        if exacts:
            return latest(exacts)

        # 2b) Keyword contains
        matches = [p for p in candidates if any(h in os.path.basename(p).lower() for h in hints)]
        if matches:
            return latest(matches)

    # 3) Fallback: take the latest modified CSV
    return latest(candidates)

def make_sim_df(csv_path: str) -> pd.DataFrame:
    """Read MF6 OBS CSV (text), construct DataFrame matching observation data:
    - Rows: Observation points (corresponding to CONFIG['obs_index_labels'])
    - Columns: Conditions ['natural', 'pumping'] (corresponding to time=1, time=2)
    """
    from flopy.utils.observationfile import Mf6Obs
    conds = CONFIG["conditions"]
    labels = CONFIG["obs_index_labels"]

    # 1. Use Mf6Obs to read (will automatically handle duplicate column names as _1, _2, etc.)
    sim_out = Mf6Obs(csv_path, isBinary=False)
    df = sim_out.get_dataframe()  # rows=time steps, columns=observation points
    
    # 2. Transpose: rows=observation points, columns=time steps (dates)
    df_t = df.transpose()
    
    # 3. Delete first row (totim row)
    df_t = df_t.iloc[1:]
    
    # 4. Remove duplicate rows (if any)
    # Check for duplicate rows caused by duplicate column names
    df_t = df_t[~df_t.index.duplicated(keep='first')]
    
    # 5. Rename columns to condition names
    if len(df_t.columns) >= len(conds):
        df_t.columns = conds[:len(df_t.columns)]
    
    # 6. Ensure observation point count matches
    if len(df_t) != len(labels):
        print(f"⚠️ Warning: OBS data has {len(df_t)} observation points, but config requires {len(labels)}")
        print(f"   OBS points: {list(df_t.index)}")
        print(f"   Config labels: {labels}")
    
    # 7. Set index to standard labels
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

    Path(runs_root).mkdir(parents=True, exist_ok=True)
    Path(outputs_dir).mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(seed)

    # 1) Sample parameters (perform LHS sampling in log space)
    # ✅ Include shapefile lithology + default formula parameters + boundary condition parameters
    keys = list(LITHOLOGY_TARGETS) + list(DEFAULT_PARAMS) + list(BOUNDARY_PARAMS)
    bounds = [PARAM_BOUNDS[k] for k in keys]
    unit = lhs(n_samples, len(keys), rng)
    
    # ✅ Handle separately: K values in log space, head values in linear space
    samples = np.empty_like(unit)
    for j, key in enumerate(keys):
        if key in BOUNDARY_PARAMS and "head" in key:
            # Head parameters in linear space
            lo, hi = bounds[j]
            samples[:, j] = lo + unit[:, j] * (hi - lo)
        else:
            # K values and conductance in log space
            lo, hi = bounds[j]
            if lo <= 0 or hi <= 0:
                raise ValueError(f"Log-space sampling requires boundary values > 0, but {key} = ({lo}, {hi})")
            log_lo = np.log10(lo)
            log_hi = np.log10(hi)
            samples[:, j] = 10 ** (log_lo + unit[:, j] * (log_hi - log_lo))
    
    # ✅ Round to three significant figures
    def round_to_n_sig_figs(x, n=3):
        """Round value to n significant figures"""
        if x == 0:
            return 0
        return round(x, -int(np.floor(np.log10(abs(x)))) + (n - 1))
    
    for i in range(samples.shape[0]):
        for j in range(samples.shape[1]):
            samples[i, j] = round_to_n_sig_figs(samples[i, j], 3)
    
    print("\n=== Sampling Statistics (Log/Linear Space LHS + Three Significant Figures) ===")
    for j, k in enumerate(keys):
        print(f"{k}:")
        print(f"  Range: [{PARAM_BOUNDS[k][0]:.3e}, {PARAM_BOUNDS[k][1]:.3e}]")
        print(f"  Samples: min={samples[:, j].min():.3e}, max={samples[:, j].max():.3e}, median={np.median(samples[:, j]):.3e}")
    print()

    # 2) Results table
    records = []

    # 3) Pre-load simulation once to parse grid and cellids (for potential info; actual runs use their own copies)
    tmp_sim = load_simulation(base_ws, mf6_exe)
    cellids = load_cellids_by_lithology(tmp_sim)

    for idx in range(n_samples):
        run_id = f"run_{idx+1:04d}"
        run_ws = os.path.join(runs_root, run_id)
        print(f"\n=== Starting {run_id} ===")

        # Assemble parameter mapping
        kmap = {k: float(samples[idx, j]) for j, k in enumerate(keys)}
        print("Parameters:")
        for k, v in kmap.items():
            if "head" in k or k == "initial_head":
                print(f"  {k}: {v:.3f}")
            else:
                print(f"  {k}: {v:.3e}")

        # Copy model to run directory
        copy_model_to_run(base_ws, run_ws)

        # Load and write parameters
        sim = load_simulation(run_ws, mf6_exe)
        set_k_for_zones(sim, kmap, cellids)
        # ✅ Set boundary condition parameters
        set_boundary_conditions(sim, kmap)
        # ✅ Relax solver settings to improve convergence success rate
        # relax_solver_settings(sim)  # ⚠️ Temporarily disabled, solver modifications have issues

        # Run
        t0 = time.time()
        ok = run_simulation(sim, run_ws)
        dt = time.time() - t0

        rmse, r2 = float("nan"), float("nan")
        sim_obs_rows = 0
        obs_csv = None  # 🆕 Initialize
        skip_this_run = False  # 🆕 Failure flag
        
        if ok:
            # Read OBS → construct sim_df (consistent with your calculation method)
            obs_csv = find_obs_csv(run_ws, hint=obs_hint)
            if obs_csv and os.path.exists(obs_csv):
                sim_df = make_sim_df(obs_csv)
                sim_obs_rows = len(sim_df)
                
                # 🆕 Check for failed observation points (MODFLOW failure marker is -1e+30)
                for condition in CONFIG["conditions"]:
                    if condition in sim_df.columns:
                        simulated = sim_df[condition].values.astype(float)
                        if (simulated < -1e+20).any():
                            print(f"⚠️ Detected failed observation points in {condition} condition, skipping this run")
                            skip_this_run = True
                            break
                
                if not skip_this_run:
                    # Calculate metrics for each condition
                    metrics = compare_to_observations(sim_df)
                    # Extract (if exists)
                    rmse_nat, r2_nat = metrics.get("natural", (float("nan"), float("nan")))
                    rmse_pmp, r2_pmp = metrics.get("pumping", (float("nan"), float("nan")))
                else:
                    rmse_nat = r2_nat = rmse_pmp = r2_pmp = float("nan")
            else:
                print("OBS CSV file not found, please check OBS package output settings or CONFIG['obs_csv_name_hint'].")
                skip_this_run = True
                rmse_nat = r2_nat = rmse_pmp = r2_pmp = float("nan")
        else:
            print("Simulation failed, recording NaN metrics.")
            skip_this_run = True
            rmse_nat = r2_nat = rmse_pmp = r2_pmp = float("nan")
        
        # 🆕 If failure detected, skip saving and cleanup, proceed to next run
        if skip_this_run:
            print(f"⏭️ Skipping {run_id}, not saving results")
            continue
        
        # 🆕 Clean up run folder (keep only observation data CSV)
        cleanup_run_folder(run_ws, obs_csv)

        # Record
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

        # Save parameters and metrics to CSV for each run (for monitoring progress)
        pd.DataFrame(records).to_csv(os.path.join(outputs_dir, "results_running.csv"), index=False)

    # Output Excel after all runs complete
    results_df = pd.DataFrame(records)
    xlsx_path = os.path.join(outputs_dir, "results.xlsx")
    with pd.ExcelWriter(xlsx_path, engine="xlsxwriter") as writer:
        results_df.to_excel(writer, index=False, sheet_name="runs")
    print(f"\nAll completed. Results written to: {xlsx_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("Error occurred:", e)
        sys.exit(1)
