# -*- coding: utf-8 -*-
"""
Comparison of Different Interpolation Methods for Transient Model Calibration

Interpolation Methods:
- Nearest: Nearest neighbor (no interpolation)
- Linear: Linear interpolation
- Cubic: Cubic spline interpolation
- PCHIP: Piecewise Cubic Hermite Interpolating Polynomial (shape-preserving)
- Akima: Akima spline interpolation
"""

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d, PchipInterpolator, Akima1DInterpolator
import os
from pathlib import Path

plt.rcParams['figure.dpi'] = 200
plt.rcParams['font.size'] = 9

# ========================================================================
# Path Configuration
# ========================================================================
# Get the script directory and project root
SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = SCRIPT_DIR.parent.parent  # Model for Exam folder

# Relative paths
directory_sim = PROJECT_ROOT / "02_Transient_Calibrate" / "Output" / "run_test"
# Alternative paths:
# directory_sim = PROJECT_ROOT / "03_Dewatering" / "Minimum_Time"
# directory_sim = PROJECT_ROOT / "02_Transient_Calibrate" / "Remodel"
directory_obs = PROJECT_ROOT / "02_Transient_Calibrate"
model = 'Transient_cali'
obs_csv = 'calibration_transient.csv'

# ========================================================================
# Read Data
# ========================================================================

sim_file = os.path.join(directory_sim, f'{model}.ob_gw_out_head.csv')
sim_df = pd.read_csv(sim_file)

obs_file = os.path.join(directory_obs, obs_csv)
obs_df = pd.read_csv(obs_file, delimiter=';')
if obs_df.shape[1] > 3:
    obs_df = obs_df.iloc[:, :3]
obs_df.columns = [col.split(',')[0].strip() for col in obs_df.columns]
for col in obs_df.columns:
    if obs_df[col].dtype == 'object':
        obs_df[col] = obs_df[col].astype(str).str.split(',').str[0]
        obs_df[col] = pd.to_numeric(obs_df[col], errors='coerce')

# ========================================================================
# Convert to Drawdown
# ========================================================================

point_mapping = {
    'P1': 'HD__1_63_61',
    'Pz12': 'HD__2_63_61'
}

for obs_point, sim_col in point_mapping.items():
    initial_head = sim_df[sim_col].iloc[0]
    sim_df[sim_col + '_drawdown'] = initial_head - sim_df[sim_col]

obs_times = obs_df['time (s)'].values

# ========================================================================
# Define Different Interpolation Methods
# ========================================================================

def calculate_metrics(observed, simulated):
    """Calculate RMSE and R²"""
    residuals = simulated - observed
    rmse = np.sqrt(np.mean(residuals**2))
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((observed - np.mean(observed))**2)
    r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    return rmse, r2

interpolation_methods = {
    'Nearest': None,  # Special marker for nearest neighbor
    'Linear': 'linear',
    'Cubic': 'cubic',
    'PCHIP': 'pchip',
    'Akima': 'akima'
}

# ========================================================================
# Create plots for each observation point
# ========================================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle('Comparison of Interpolation Methods - P1 and Pz12', fontsize=14, fontweight='bold')

for idx, (obs_point, sim_col) in enumerate(point_mapping.items()):
    sim_col_dd = sim_col + '_drawdown'
    
    # Extract data
    sim_time = sim_df['time'].values
    sim_drawdown = sim_df[sim_col_dd].values
    obs_drawdown = obs_df[obs_point].values
    
    # Store results from different methods
    results = {}
    
    for method_name, method_type in interpolation_methods.items():
        if method_type is None:
            # Nearest neighbor: find closest simulation time point
            interp_values = np.array([
                sim_drawdown[np.argmin(np.abs(sim_time - t))] 
                for t in obs_times
            ])
        elif method_type == 'pchip':
            interpolator = PchipInterpolator(sim_time, sim_drawdown)
            interp_values = interpolator(obs_times)
        elif method_type == 'akima':
            interpolator = Akima1DInterpolator(sim_time, sim_drawdown)
            interp_values = interpolator(obs_times)
        else:
            interpolator = interp1d(sim_time, sim_drawdown, kind=method_type, 
                                   bounds_error=False, fill_value='extrapolate')
            interp_values = interpolator(obs_times)
        
        rmse, r2 = calculate_metrics(obs_drawdown, interp_values)
        results[method_name] = {
            'values': interp_values,
            'rmse': rmse,
            'r2': r2
        }
    
    # Plot time series
    ax = axes[idx]
    ax.plot(obs_times, obs_drawdown, 'o', color='red', 
                 markersize=6, label='Observed', zorder=10)
    ax.plot(sim_time, sim_drawdown, 's', color='gray', 
                 markersize=3, alpha=0.5, label='Simulated timesteps')
    
    colors = ['blue', 'green', 'orange', 'purple', 'brown']
    linestyles = ['-', '--', '-.', ':', '-']
    
    for i, (method_name, result) in enumerate(results.items()):
        ax.plot(obs_times, result['values'], 
                    color=colors[i], linestyle=linestyles[i], 
                    linewidth=1.5, alpha=0.7,
                    label=f"{method_name} (R²={result['r2']:.3f})")
    
    ax.set_xlabel('Time (seconds)', fontsize=11)
    ax.set_ylabel('Drawdown (meters)', fontsize=11)
    ax.set_title(f'{obs_point} - Time Series Comparison', fontsize=12, fontweight='bold')
    ax.set_xscale('linear')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc='best')
    
    # Print results
    print(f"\n{'='*70}")
    print(f"{obs_point} - Evaluation Metrics for Different Interpolation Methods:")
    print(f"{'='*70}")
    print(f"{'Method':<20} {'RMSE (m)':<15} {'R²':<15}")
    print("-"*70)
    for method_name, result in results.items():
        print(f"{method_name:<20} {result['rmse']:<15.4f} {result['r2']:<15.4f}")

plt.tight_layout()
output_file = os.path.join(directory_obs, 'interpolation_comparison.png')
plt.savefig(output_file, dpi=200, bbox_inches='tight')
print(f"\n✅ Figure saved: {output_file}")
plt.show()

print("\n" + "="*70)
print("Conclusions:")
print("="*70)
print("1. If all methods have similar R², choose 'Nearest' to avoid interpolation uncertainty")
print("2. If Linear interpolation has better R², data changes smoothly")
print("3. If PCHIP or Akima have better R², data has nonlinear characteristics")
print("4. Note: Systematic bias in simulated values indicates model parameter issues, not interpolation")
print("="*70)
