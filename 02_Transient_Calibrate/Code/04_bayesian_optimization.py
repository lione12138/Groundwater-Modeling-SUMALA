#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Bayesian Optimization for Transient Model Calibration
=====================================================
Use existing calibration data to train a surrogate model (Gaussian Process),
then use Bayesian Optimization to find optimal parameters efficiently.

Advantages:
- Learns from existing 1100+ runs
- Intelligently explores parameter space
- Balances exploration vs exploitation
- Much faster than random/LHS sampling
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
import time
import shutil
import glob
import os
from datetime import datetime

# Bayesian Optimization
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, Matern
from sklearn.preprocessing import StandardScaler
from scipy.optimize import minimize
from scipy.stats import norm

import matplotlib.pyplot as plt
import seaborn as sns

# MODFLOW
import flopy
from flopy.mf6 import MFSimulation

# Spatial libraries
try:
    import fiona
    from shapely.geometry import shape, Point
    from flopy.utils import GridIntersect
    HAS_SPATIAL_LIBS = True
except ImportError:
    HAS_SPATIAL_LIBS = False
    print("Warning: fiona and/or shapely not available.")

# Configuration
ROOT = Path(__file__).resolve().parent
CONFIG = {
    "existing_results": str((ROOT.parent / "Output" / "lhs_runs" / "calibration_results.csv").resolve()),
    "output_dir": str((ROOT.parent / "Bayesian_Optimization").resolve()),
    "n_iterations": 50,  # Number of Bayesian Optimization iterations
    "n_restarts": 10,    # Random restarts for acquisition function optimization
    
    # MODFLOW Configuration (same as 01_lhs_calibration.py)
    "base_model_ws": str((ROOT.parent / "Model").resolve()),
    "mf6_exe": r"D:\Program Files\USGS\mf6.6.2_win64\bin\mf6.exe",
    "runs_root": str((ROOT.parent / "Bayesian_Optimization" / "runs").resolve()),
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

# Parameter bounds (same as LHS_CONFIG)
PARAM_BOUNDS = {
    # Sy bounds (linear scale)
    "sy_clean_gravel": [0.20, 0.35],    # High permeability → high Sy
    "sy_clayey_sand": [0.10, 0.25],     # Clay reduces Sy
    "sy_loamy_sand": [0.15, 0.30],      # Medium Sy
    "sy_sandy_gravel": [0.1, 0.35],    # High Sy
    "sy_sandstone": [0.001, 0.001],     # Confined layer, Sy irrelevant (fixed)
    "sy_limestone": [0.001, 0.001],     # Confined layer, Sy irrelevant (fixed)

    # Ss bounds (log scale, but we'll work in linear space)
    "ss_clean_gravel": [-5.5, -4.5],    # 10^-5.5 to 10^-4.5 (3.16e-6 to 3.16e-5)
    "ss_clayey_sand": [-5.3, -4.7],     # Higher due to clay compressibility
    "ss_loamy_sand": [-5.4, -4.6],      # Medium range
    "ss_sandy_gravel": [-5.5, -4.5],    # Similar to clean gravel
    "ss_sandstone": [-5.0, -4.0],       # Layer 2: critical parameter (1e-5 to 1e-4)
    "ss_limestone": [-6.0, -5.0],       # Lower compressibility (1e-6 to 1e-5)
}

TARGET_METRIC = "Overall_R2"  # Maximize this

LITHOLOGY_TARGETS = [
    "clean_gravel",
    "clayey_sand",
    "loamy_sand",
    "sandstone",
    "sandy_gravel",
    "limestone"
]


# ============================================================================
# MODFLOW Helper Functions (from 01_lhs_calibration.py)
# ============================================================================

def load_cellids_by_lithology(sim: MFSimulation, shapefile_path: str,
                               name_field: str, value_to_key: dict,
                               layer_by_key: dict):
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


def apply_storage_parameters(sim: MFSimulation, params: dict, cellids_by_lith: dict):
    """Apply Sy and Ss parameters to storage package."""
    gwf = sim.get_model()
    sto = gwf.sto
    dis = gwf.dis
    
    nlay, nrow, ncol = dis.nlay.array, dis.nrow.array, dis.ncol.array
    
    sy_arrays = {}
    ss_arrays = {}
    
    for lay in range(nlay):
        sy_arrays[lay] = np.full((nrow, ncol), 0.2)
        ss_arrays[lay] = np.full((nrow, ncol), 1e-5)
    
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
    
    sto.sy.set_data(sy_arrays)
    sto.ss.set_data(ss_arrays)
    sim.write_simulation()


def run_model(sim: MFSimulation, mf6_exe: str, verbose: bool = False):
    """Run MODFLOW6 simulation."""
    success, _ = sim.run_simulation(silent=not verbose)
    return success


def extract_drawdown(sim: MFSimulation, obs_csv_hint: str):
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
    
    obs_mapping = CONFIG.get("obs_name_mapping", {})
    
    for mf_name, user_name in obs_mapping.items():
        if mf_name in sim_df.columns:
            initial_head = sim_df[mf_name].iloc[0]
            drawdown_df[user_name] = initial_head - sim_df[mf_name]
    
    return drawdown_df


def compare_to_observations(sim_df: pd.DataFrame, obs_csv_path: str):
    """Compare simulation to observed data and calculate metrics."""
    obs_df = pd.read_csv(obs_csv_path)
    
    if obs_df.columns[0].count(';') > 0:
        obs_df = obs_df.iloc[:, 1:]
    
    if 'time (s)' not in obs_df.columns and obs_df.columns[0] != 'time (s)':
        obs_df.rename(columns={obs_df.columns[0]: 'time (s)'}, inplace=True)
    
    for col in obs_df.columns:
        obs_df[col] = pd.to_numeric(obs_df[col], errors='coerce')
    
    obs_df = obs_df.dropna()
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
        
        metrics[f"{obs_point}_RMSE"] = rmse
        metrics[f"{obs_point}_R2"] = r2
    
    rmse_values = [v for k, v in metrics.items() if 'RMSE' in k]
    r2_values = [v for k, v in metrics.items() if 'R2' in k]
    
    if rmse_values:
        metrics['Overall_RMSE'] = np.mean(rmse_values)
        metrics['Overall_R2'] = np.mean(r2_values)
    
    return metrics


def run_modflow_with_params(params: dict, iteration: int, cellids_by_lith: dict, debug: bool = False):
    """
    Run MODFLOW with given parameters and return performance metrics.
    Returns (success: bool, metric_value: float, all_metrics: dict, error_msg: str)
    """
    iter_ws = Path(CONFIG["runs_root"]) / f"bo_iter_{iteration:03d}"
    error_msg = ""
    
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
            pass
        time.sleep(0.1)
    
    try:
        shutil.copytree(CONFIG["base_model_ws"], iter_ws)
    except Exception as e:
        error_msg = f"Copy failed: {e}"
        if debug:
            print(f"\n  DEBUG: {error_msg}")
        return False, np.nan, {}, error_msg
    
    try:
        # Load simulation
        sim = MFSimulation.load(
            sim_ws=str(iter_ws),
            exe_name=CONFIG["mf6_exe"],
            verbosity_level=0
        )
        
        # Apply parameters
        apply_storage_parameters(sim, params, cellids_by_lith)
        
        # Run model
        success = run_model(sim, CONFIG["mf6_exe"], verbose=debug)
        
        if not success:
            error_msg = "MODFLOW run failed"
            if debug:
                # Check MODFLOW list file for errors
                list_files = list(iter_ws.glob("*.lst"))
                if list_files:
                    with open(list_files[0], 'r') as f:
                        lines = f.readlines()
                        for i, line in enumerate(lines[-50:]):  # Last 50 lines
                            if 'ERROR' in line.upper() or 'FAIL' in line.upper():
                                print(f"\n  DEBUG MODFLOW: {line.strip()}")
            return False, np.nan, {}, error_msg
        
        # Extract and compare results
        try:
            sim_drawdown = extract_drawdown(sim, CONFIG["obs_csv_name_hint"])
        except Exception as e:
            error_msg = f"Extract failed: {e}"
            if debug:
                print(f"\n  DEBUG: {error_msg}")
            return False, np.nan, {}, error_msg
        
        if len(sim_drawdown.columns) <= 1:
            error_msg = "No observation points in output"
            if debug:
                print(f"\n  DEBUG: {error_msg}")
            return False, np.nan, {}, error_msg
        
        try:
            metrics = compare_to_observations(sim_drawdown, CONFIG["observed_csv_path"])
        except Exception as e:
            error_msg = f"Metrics calculation failed: {e}"
            if debug:
                print(f"\n  DEBUG: {error_msg}")
            return False, np.nan, {}, error_msg
        
        if not metrics or TARGET_METRIC not in metrics:
            error_msg = f"Target metric '{TARGET_METRIC}' not found"
            if debug:
                print(f"\n  DEBUG: {error_msg}, available: {list(metrics.keys())}")
            return False, np.nan, metrics, error_msg
        
        # Cleanup
        try:
            for item in iter_ws.iterdir():
                if item.is_file() and not item.name.endswith('.csv'):
                    try:
                        os.chmod(item, 0o777)
                        item.unlink()
                    except Exception:
                        pass
                elif item.is_dir():
                    try:
                        shutil.rmtree(item, ignore_errors=True)
                    except Exception:
                        pass
        except Exception:
            pass
        
        return True, metrics[TARGET_METRIC], metrics, ""
        
    except Exception as e:
        error_msg = f"Unexpected error: {e}"
        if debug:
            import traceback
            print(f"\n  DEBUG TRACEBACK:\n{traceback.format_exc()}")
        return False, np.nan, {}, error_msg


# ============================================================================
# Bayesian Optimization Class
# ============================================================================

class BayesianOptimizer:
    """Bayesian Optimization with Gaussian Process surrogate model."""
    
    def __init__(self, param_bounds, target_metric="Overall_R2"):
        self.param_bounds = param_bounds
        self.target_metric = target_metric
        self.param_names = list(param_bounds.keys())
        
        # Filter out fixed parameters
        self.variable_params = {k: v for k, v in param_bounds.items() 
                                if v[0] != v[1]}
        self.variable_param_names = list(self.variable_params.keys())
        self.n_params = len(self.variable_param_names)
        
        # Scalers for normalization
        self.X_scaler = StandardScaler()
        self.y_scaler = StandardScaler()
        
        # Gaussian Process model
        kernel = ConstantKernel(1.0) * Matern(length_scale=1.0, nu=2.5)
        self.gp = GaussianProcessRegressor(
            kernel=kernel,
            alpha=1e-6,
            normalize_y=False,  # We do our own normalization
            n_restarts_optimizer=5,
            random_state=42
        )
        
        # Training data
        self.X_train = None
        self.y_train = None
        self.best_y = -np.inf
        self.best_x = None
        
    def load_existing_data(self, csv_path):
        """Load and prepare existing calibration results."""
        df = pd.read_csv(csv_path)
        df = df[df['success'] == True].copy()
        
        # Extract variable parameters
        X = df[self.variable_param_names].values
        y = df[self.target_metric].values
        
        # Remove NaN values
        valid_idx = ~np.isnan(y)
        X = X[valid_idx]
        y = y[valid_idx]
        
        print(f"Loaded {len(X)} successful runs with valid {self.target_metric}")
        print(f"Current best {self.target_metric}: {y.max():.6f}")
        print(f"Variable parameters: {self.n_params}")
        
        return X, y
    
    def fit(self, X, y):
        """Fit Gaussian Process on existing data."""
        # Normalize features and target
        self.X_train = self.X_scaler.fit_transform(X)
        self.y_train = self.y_scaler.fit_transform(y.reshape(-1, 1)).ravel()
        
        # Fit GP
        print("\nFitting Gaussian Process...")
        self.gp.fit(self.X_train, self.y_train)
        print(f"✓ GP fitted with {len(X)} samples")
        
        # Track best observed value
        self.best_y = y.max()
        self.best_x = X[y.argmax()]
        
        # Calculate training R²
        y_pred = self.predict(X)
        r2 = 1 - np.sum((y - y_pred)**2) / np.sum((y - y.mean())**2)
        print(f"✓ GP Training R²: {r2:.4f}")
        
    def predict(self, X, return_std=False):
        """Predict using GP (returns in original scale)."""
        X_scaled = self.X_scaler.transform(X)
        
        if return_std:
            y_pred, y_std = self.gp.predict(X_scaled, return_std=True)
            # Transform back to original scale
            y_pred = self.y_scaler.inverse_transform(y_pred.reshape(-1, 1)).ravel()
            # Approximate std in original scale
            y_std = y_std * self.y_scaler.scale_[0]
            return y_pred, y_std
        else:
            y_pred = self.gp.predict(X_scaled)
            y_pred = self.y_scaler.inverse_transform(y_pred.reshape(-1, 1)).ravel()
            return y_pred
    
    def acquisition_function(self, X, xi=0.01):
        """Expected Improvement (EI) acquisition function."""
        mu, sigma = self.predict(X, return_std=True)
        
        # Avoid division by zero
        sigma = np.maximum(sigma, 1e-9)
        
        # Expected Improvement
        improvement = mu - self.best_y - xi
        Z = improvement / sigma
        ei = improvement * norm.cdf(Z) + sigma * norm.pdf(Z)
        
        return ei
    
    def propose_location(self, n_restarts=10, xi=0.01):
        """Propose next sampling location by maximizing acquisition function."""
        best_ei = -np.inf
        best_x = None
        
        # Get bounds for variable parameters
        bounds = np.array([self.variable_params[name] for name in self.variable_param_names])
        
        # Multiple random restarts
        for restart in range(n_restarts):
            # Random starting point
            x0 = np.random.uniform(bounds[:, 0], bounds[:, 1])
            
            # Maximize acquisition function (minimize negative)
            def neg_acquisition(x):
                return -self.acquisition_function(x.reshape(1, -1), xi=xi)[0]
            
            result = minimize(
                neg_acquisition,
                x0,
                bounds=bounds,
                method='L-BFGS-B'
            )
            
            if -result.fun > best_ei:
                best_ei = -result.fun
                best_x = result.x
        
        # Round to 3 significant figures
        best_x = np.array([float(f'{x:.3g}') for x in best_x])
        
        return best_x, best_ei
    
    def update(self, X_new, y_new):
        """Update GP with new observations."""
        # Add to training data
        X_new_scaled = self.X_scaler.transform(X_new.reshape(1, -1))
        y_new_scaled = self.y_scaler.transform(np.array([[y_new]])).ravel()
        
        self.X_train = np.vstack([self.X_train, X_new_scaled])
        self.y_train = np.append(self.y_train, y_new_scaled)
        
        # Update best observed
        if y_new > self.best_y:
            self.best_y = y_new
            self.best_x = X_new
        
        # Refit GP
        self.gp.fit(self.X_train, self.y_train)
    
    def params_to_dict(self, x_array):
        """Convert parameter array to dictionary (including fixed params)."""
        params = {}
        var_idx = 0
        for name in self.param_names:
            if name in self.variable_param_names:
                value = float(x_array[var_idx])
                
                # Convert log-scale Ss values to linear scale
                if name.startswith('ss_'):
                    value = 10 ** value  # Convert from log to linear
                
                params[name] = value
                var_idx += 1
            else:
                # Fixed parameter
                fixed_value = self.param_bounds[name][0]
                # Also convert fixed Ss values from log to linear if needed
                if name.startswith('ss_') and fixed_value < 0:
                    fixed_value = 10 ** fixed_value
                params[name] = fixed_value
        return params


def plot_optimization_progress(results, output_path, initial_best):
    """Plot Bayesian Optimization progress."""
    fig, axes = plt.subplots(3, 1, figsize=(12, 14))
    
    iterations = [r['iteration'] for r in results]
    ei_values = [r['expected_improvement'] for r in results]
    predicted_values = [r['predicted_metric'] for r in results]
    actual_values = [r.get('actual_metric', np.nan) for r in results]
    
    # Cumulative best
    cumulative_best = []
    best_so_far = initial_best
    for actual in actual_values:
        if not np.isnan(actual):
            best_so_far = max(best_so_far, actual)
        cumulative_best.append(best_so_far)
    
    # Plot 1: Expected Improvement
    ax1 = axes[0]
    ax1.plot(iterations, ei_values, 'b-', linewidth=2, label='Expected Improvement')
    ax1.fill_between(iterations, 0, ei_values, alpha=0.3)
    ax1.set_xlabel('Iteration', fontsize=12)
    ax1.set_ylabel('Expected Improvement', fontsize=12)
    ax1.set_title('Acquisition Function - Expected Improvement', fontsize=14, fontweight='bold')
    ax1.grid(alpha=0.3)
    ax1.legend()
    
    # Plot 2: Predicted vs Actual
    ax2 = axes[1]
    ax2.plot(iterations, predicted_values, 'g--', linewidth=2, marker='o', markersize=4,
             label='GP Predicted', alpha=0.7)
    valid_idx = [i for i, v in enumerate(actual_values) if not np.isnan(v)]
    ax2.plot([iterations[i] for i in valid_idx], [actual_values[i] for i in valid_idx],
             'r-', linewidth=2, marker='s', markersize=6, label='Actual (MODFLOW)')
    ax2.axhline(y=initial_best, color='gray', linestyle=':', linewidth=2,
                label=f'Initial Best: {initial_best:.4f}')
    ax2.set_xlabel('Iteration', fontsize=12)
    ax2.set_ylabel('Overall R²', fontsize=12)
    ax2.set_title('Predicted vs Actual Performance', fontsize=14, fontweight='bold')
    ax2.grid(alpha=0.3)
    ax2.legend()
    
    # Plot 3: Cumulative Best
    ax3 = axes[2]
    ax3.plot(iterations, cumulative_best, 'purple', linewidth=3, marker='D', markersize=5,
             label='Best R² Found')
    ax3.axhline(y=initial_best, color='gray', linestyle=':', linewidth=2,
                label=f'Initial Best: {initial_best:.4f}')
    ax3.fill_between(iterations, initial_best, cumulative_best, alpha=0.3, color='purple')
    ax3.set_xlabel('Iteration', fontsize=12)
    ax3.set_ylabel('Best Overall R²', fontsize=12)
    ax3.set_title('Cumulative Best Performance', fontsize=14, fontweight='bold')
    ax3.grid(alpha=0.3)
    ax3.legend()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_parameter_evolution(results, output_path):
    """Plot how parameters evolve during optimization."""
    # Extract variable parameters
    iterations = [r['iteration'] for r in results]
    param_names = [k for k in results[0]['parameters'].keys() 
                   if k not in ['sy_sandstone', 'sy_limestone']]  # Skip fixed params
    
    n_params = len(param_names)
    n_cols = 3
    n_rows = (n_params + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4*n_rows))
    axes = axes.flatten() if n_params > 1 else [axes]
    
    for idx, param in enumerate(param_names):
        values = [r['parameters'][param] for r in results]
        
        ax = axes[idx]
        ax.plot(iterations, values, 'o-', linewidth=2, markersize=4)
        ax.set_xlabel('Iteration', fontsize=10)
        ax.set_ylabel('Parameter Value', fontsize=10)
        ax.set_title(param.replace('_', ' ').title(), fontsize=11, fontweight='bold')
        ax.grid(alpha=0.3)
    
    # Hide unused subplots
    for idx in range(n_params, len(axes)):
        axes[idx].axis('off')
    
    plt.suptitle('Parameter Evolution During Bayesian Optimization', 
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def main():
    """Main Bayesian Optimization execution."""
    print("\n" + "="*80)
    print("BAYESIAN OPTIMIZATION FOR TRANSIENT MODEL CALIBRATION")
    print("="*80 + "\n")
    
    start_time = time.time()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Create output directory
    output_dir = Path(CONFIG["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create runs directory
    runs_root = Path(CONFIG["runs_root"])
    runs_root.mkdir(parents=True, exist_ok=True)
    
    # Initialize optimizer
    print("Initializing Bayesian Optimizer...")
    optimizer = BayesianOptimizer(PARAM_BOUNDS, TARGET_METRIC)
    print()
    
    # Load existing data
    print("Loading existing calibration results...")
    X_existing, y_existing = optimizer.load_existing_data(CONFIG["existing_results"])
    print()
    
    # Fit GP on existing data
    optimizer.fit(X_existing, y_existing)
    print()
    
    # Load lithology mapping
    print("Loading lithology mapping...")
    temp_ws = runs_root / "temp_mapping"
    if temp_ws.exists():
        import stat
        def remove_readonly(func, path, _):
            try:
                os.chmod(path, stat.S_IWRITE)
                func(path)
            except Exception:
                pass
        try:
            shutil.rmtree(temp_ws, onerror=remove_readonly)
        except Exception:
            pass
        time.sleep(0.1)
    
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
        print(f"✓ Loaded lithology mapping for {len(cellids_by_lith)} lithologies")
    
    # Clean up temp directory
    import stat
    def remove_readonly(func, path, _):
        try:
            os.chmod(path, stat.S_IWRITE)
            func(path)
        except Exception:
            pass
    try:
        shutil.rmtree(temp_ws, onerror=remove_readonly)
    except Exception:
        pass
    print()
    
    # Bayesian Optimization iterations
    print("="*80)
    print(f"Starting Bayesian Optimization ({CONFIG['n_iterations']} iterations)")
    print("="*80 + "\n")
    
    results = []
    iteration_counter = 0
    failed_count = 0
    
    for i in range(CONFIG['n_iterations']):
        iter_start = time.time()
        print(f"[{i+1}/{CONFIG['n_iterations']}] ", end='', flush=True)
        
        # Propose next point with exploration parameter
        # Start with more exploration (higher xi), gradually decrease
        xi = 0.01 + 0.05 * (1 - i / CONFIG['n_iterations'])
        x_next, ei = optimizer.propose_location(n_restarts=CONFIG['n_restarts'], xi=xi)
        
        # Predict performance
        y_pred, y_std = optimizer.predict(x_next.reshape(1, -1), return_std=True)
        y_pred = y_pred[0]
        y_std = y_std[0]
        
        print(f"EI={ei:.6f}, Pred R²={y_pred:.4f}±{y_std:.4f} ", end='', flush=True)
        
        # Convert to full parameter dict
        params = optimizer.params_to_dict(x_next)
        
        # Run MODFLOW - enable debug for first 3 failures
        debug = (failed_count < 3)
        success, actual_metric, all_metrics, error_msg = run_modflow_with_params(
            params, iteration_counter, cellids_by_lith, debug=debug
        )
        iteration_counter += 1
        
        if success:
            # Update GP with actual result
            optimizer.update(x_next, actual_metric)
            
            iter_time = time.time() - iter_start
            improvement = actual_metric - optimizer.best_y if actual_metric > optimizer.best_y else 0
            
            print(f"→ Actual R²={actual_metric:.4f} ", end='')
            if improvement > 0:
                print(f"✓ NEW BEST! (+{improvement:.4f}) [{iter_time:.1f}s]")
            else:
                print(f"[{iter_time:.1f}s]")
            
            # Store results
            result = {
                'iteration': i,
                'expected_improvement': float(ei),
                'predicted_metric': float(y_pred),
                'predicted_std': float(y_std),
                'actual_metric': float(actual_metric),
                'is_new_best': improvement > 0,
                'error': '',
                'parameters': params,
                **all_metrics
            }
        else:
            failed_count += 1
            print(f"✗ FAILED: {error_msg} [{time.time() - iter_start:.1f}s]")
            result = {
                'iteration': i,
                'expected_improvement': float(ei),
                'predicted_metric': float(y_pred),
                'predicted_std': float(y_std),
                'actual_metric': np.nan,
                'is_new_best': False,
                'error': error_msg,
                'parameters': params
            }
            
            # Stop if too many consecutive failures
            if failed_count >= 10:
                print(f"\n\n⚠️ Stopping after {failed_count} consecutive failures.")
                print("Please check the model configuration and parameters.")
                break
        
        results.append(result)
        
        # Save intermediate results every 5 iterations
        if (i + 1) % 5 == 0:
            temp_df = pd.DataFrame(results)
            temp_df.to_csv(output_dir / f"bo_results_iter{i+1}.csv", index=False)
    
    print()
    
    # Find best actual result
    successful_results = [r for r in results if not np.isnan(r.get('actual_metric', np.nan))]
    
    if successful_results:
        best_result = max(successful_results, key=lambda x: x['actual_metric'])
        best_from_bo = best_result['actual_metric']
    else:
        best_result = max(results, key=lambda x: x['predicted_metric'])
        best_from_bo = best_result['predicted_metric']
    
    print("="*80)
    print("OPTIMIZATION RESULTS")
    print("="*80)
    print(f"\nInitial best (from existing data):")
    print(f"  {TARGET_METRIC}: {y_existing.max():.6f}")
    print(f"\nFinal best (after Bayesian Optimization):")
    print(f"  {TARGET_METRIC}: {optimizer.best_y:.6f}")
    
    improvement = optimizer.best_y - y_existing.max()
    print(f"\n{'✓' if improvement > 0 else '→'} Improvement: {improvement:+.6f} ({improvement/y_existing.max()*100:+.2f}%)")
    
    print(f"\nSuccessful runs: {len(successful_results)}/{len(results)}")
    
    print(f"\n--- Optimal Parameters (Predicted) ---")
    print("\nSpecific Yield (Sy):")
    for key, val in best_result['parameters'].items():
        if key.startswith('sy_'):
            lith = key.replace('sy_', '')
            print(f"  {lith:20s}: {val:.4f}")
    
    print("\nSpecific Storage (Ss):")
    for key, val in best_result['parameters'].items():
        if key.startswith('ss_'):
            lith = key.replace('ss_', '')
            print(f"  {lith:20s}: {val:.2e} (1/m)")
    
    # Save results
    results_df = pd.DataFrame(results)
    
    results_file = output_dir / "bayesian_optimization_results.csv"
    results_df.to_csv(results_file, index=False)
    print(f"\n✓ Results saved to: {results_file}")
    
    # Save best parameters
    best_params = {lith: {} for lith in LITHOLOGY_TARGETS}
    for lith in LITHOLOGY_TARGETS:
        sy_key = f"sy_{lith}"
        ss_key = f"ss_{lith}"
        if sy_key in best_result['parameters']:
            best_params[lith]['Sy'] = float(best_result['parameters'][sy_key])
        if ss_key in best_result['parameters']:
            best_params[lith]['Ss'] = float(best_result['parameters'][ss_key])
    
    best_params_file = output_dir / "best_parameters_bayesian.json"
    with open(best_params_file, 'w') as f:
        json.dump({
            'timestamp': timestamp,
            'method': 'Bayesian Optimization',
            'n_existing_samples': len(X_existing),
            'n_bo_iterations': CONFIG['n_iterations'],
            'n_successful': len(successful_results),
            'target_metric': TARGET_METRIC,
            'initial_best': float(y_existing.max()),
            'final_best': float(optimizer.best_y),
            'improvement': float(improvement),
            'improvement_percent': float(improvement/y_existing.max()*100),
            'best_actual_metric': float(best_result.get('actual_metric', best_result['predicted_metric'])),
            'parameters': best_params
        }, f, indent=2)
    print(f"✓ Best parameters saved to: {best_params_file}")
    
    # Generate plots
    print("\nGenerating visualization plots...")
    plot_optimization_progress(results, output_dir / "optimization_progress.png", y_existing.max())
    print("✓ Saved: optimization_progress.png")
    
    plot_parameter_evolution(results, output_dir / "parameter_evolution.png")
    print("✓ Saved: parameter_evolution.png")
    
    # Summary
    elapsed = time.time() - start_time
    print(f"\n{'='*80}")
    print(f"✓ Bayesian Optimization completed in {elapsed:.1f} seconds")
    print(f"✓ All results saved to: {output_dir}")
    print(f"{'='*80}\n")
    
    print("\n" + "="*80)
    print("NEXT STEPS")
    print("="*80)
    print("\n1. Review optimization progress in 'optimization_progress.png'")
    print("2. Check best parameters in 'best_parameters_bayesian.json'")
    print("3. If needed, run more iterations by re-executing this script")
    print("   (it will load all previous results and continue optimizing)")
    print("4. Use the best parameters in your final transient model")
    print()


if __name__ == "__main__":
    main()
