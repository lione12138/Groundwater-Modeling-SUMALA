"""
Final Comparison Script for Steady-State Calibration Results
Compares three model runs: run_flopy_shp (FloPy), run_modelmuse (ModelMuse), run_flopy_gpt (FloPy Strict)
Generates publication-quality figure for report
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.metrics import r2_score, mean_squared_error

# ============================================================================
# Configuration
# ============================================================================
# Get the script directory and project root
SCRIPT_DIR = Path(__file__).parent.resolve()
BASE_DIR = SCRIPT_DIR.parent.parent  # Model for Exam folder
OUT_DIR = BASE_DIR / "01_Steady_Calibrate" / "Out"

# Three model runs to compare
RUNS = {
    'FloPy (Shapefile)': OUT_DIR / "run_flopy_shp" / "Steady_cali.ob_gw_out_head.csv",
    'ModelMuse GUI': OUT_DIR / "run_modelmuse" / "Steady_cali.ob_gw_out_head.csv",
    'FloPy (.Gpt)': OUT_DIR / "run_flopy_gpt" / "Steady_cali.ob_gw_out_head.csv",
}

# Observation data file
OBS_FILE = BASE_DIR / "Instructions" / "steady_state_calibration.csv"

# Output figure path
OUTPUT_FIG = OUT_DIR / "fig_method_comparison.png"

# Observation point names (matching CSV column order)
OBS_NAMES = ['P1', 'Pz2', 'Pz3', 'Pz4', 'Pz5', 'Pz6', 'Pz7', 'Pz8', 'Pz9', 'Pz10', 'Pz11', 'Pz12']

# ============================================================================
# Load Observation Data
# ============================================================================
def load_observations(obs_file):
    """Load observed heads from CSV file"""
    df = pd.read_csv(obs_file, sep=';', index_col=0)
    natural = df['natural'].values
    pumping = df['pumping'].values
    return natural, pumping

# ============================================================================
# Load Simulation Results
# ============================================================================
def load_simulation(sim_file):
    """Load simulated heads from MODFLOW OBS output"""
    df = pd.read_csv(sim_file)
    # Row 0 = natural (time=1), Row 1 = pumping (time=2)
    # Columns 1-12 are the 12 observation points (Layer 1), column 13 is Layer 2
    natural = df.iloc[0, 1:13].values.astype(float)
    pumping = df.iloc[1, 1:13].values.astype(float)
    return natural, pumping

# ============================================================================
# Calculate Metrics
# ============================================================================
def calculate_metrics(obs, sim):
    """Calculate R² and RMSE"""
    r2 = r2_score(obs, sim)
    rmse = np.sqrt(mean_squared_error(obs, sim))
    return r2, rmse

# ============================================================================
# Main Plotting Function
# ============================================================================
def create_comparison_figure():
    # Load observations
    obs_natural, obs_pumping = load_observations(OBS_FILE)
    
    # Load all simulations and calculate metrics
    results = {}
    for name, filepath in RUNS.items():
        sim_natural, sim_pumping = load_simulation(filepath)
        r2_nat, rmse_nat = calculate_metrics(obs_natural, sim_natural)
        r2_pump, rmse_pump = calculate_metrics(obs_pumping, sim_pumping)
        results[name] = {
            'sim_natural': sim_natural,
            'sim_pumping': sim_pumping,
            'r2_natural': r2_nat,
            'rmse_natural': rmse_nat,
            'r2_pumping': r2_pump,
            'rmse_pumping': rmse_pump,
        }
    
    # Print metrics to console
    print("\n" + "="*70)
    print("STEADY-STATE CALIBRATION COMPARISON")
    print("="*70)
    for name, data in results.items():
        print(f"\n{name}:")
        print(f"  Natural:  R² = {data['r2_natural']:.4f}, RMSE = {data['rmse_natural']:.4f} m")
        print(f"  Pumping:  R² = {data['r2_pumping']:.4f}, RMSE = {data['rmse_pumping']:.4f} m")
    
    # Calculate differences between runs
    print("\n" + "-"*70)
    print("DIFFERENCES (relative to ModelMuse GUI):")
    print("-"*70)
    gui_data = results['ModelMuse GUI']
    for name, data in results.items():
        if name != 'ModelMuse GUI':
            diff_r2_nat = (data['r2_natural'] - gui_data['r2_natural']) / gui_data['r2_natural'] * 100
            diff_rmse_nat = (data['rmse_natural'] - gui_data['rmse_natural']) / gui_data['rmse_natural'] * 100
            diff_r2_pump = (data['r2_pumping'] - gui_data['r2_pumping']) / gui_data['r2_pumping'] * 100
            diff_rmse_pump = (data['rmse_pumping'] - gui_data['rmse_pumping']) / gui_data['rmse_pumping'] * 100
            print(f"\n{name} vs GUI:")
            print(f"  Natural:  ΔR² = {diff_r2_nat:+.2f}%, ΔRMSE = {diff_rmse_nat:+.2f}%")
            print(f"  Pumping:  ΔR² = {diff_r2_pump:+.2f}%, ΔRMSE = {diff_rmse_pump:+.2f}%")
    
    # ========================================================================
    # Create Figure
    # ========================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Colors and markers for each run
    colors = {'FloPy (Shapefile)': '#2ecc71', 'ModelMuse GUI': '#3498db', 'FloPy (.Gpt)': '#e74c3c'}
    markers = {'FloPy (Shapefile)': 'o', 'ModelMuse GUI': 's', 'FloPy (.Gpt)': '^'}
    
    # Get axis limits
    all_values = np.concatenate([obs_natural, obs_pumping] + 
                                 [results[n]['sim_natural'] for n in RUNS] + 
                                 [results[n]['sim_pumping'] for n in RUNS])
    vmin, vmax = all_values.min() - 0.1, all_values.max() + 0.1
    
    # ========================================================================
    # Left Panel: Natural Condition
    # ========================================================================
    ax1 = axes[0]
    ax1.plot([vmin, vmax], [vmin, vmax], 'k--', lw=1.5, label='1:1 Line', alpha=0.7)
    
    for name, data in results.items():
        ax1.scatter(obs_natural, data['sim_natural'], 
                   c=colors[name], marker=markers[name], s=100, 
                   edgecolors='white', linewidths=0.5,
                   label=f"{name}", alpha=0.85, zorder=5)
    
    ax1.set_xlabel('Observed Head (m)', fontsize=12)
    ax1.set_ylabel('Simulated Head (m)', fontsize=12)
    ax1.set_title('(a) Natural Condition (No Pumping)', fontsize=13, fontweight='bold')
    ax1.set_xlim(vmin, vmax)
    ax1.set_ylim(vmin, vmax)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper left', fontsize=9)
    
    # Add metrics text box for Natural
    gui_nat = results['ModelMuse GUI']
    shp_nat = results['FloPy (Shapefile)']
    gpt_nat = results['FloPy (.Gpt)']
    
    # Calculate percentage changes relative to GUI
    shp_r2_nat_pct = (shp_nat['r2_natural'] - gui_nat['r2_natural']) / gui_nat['r2_natural'] * 100
    shp_rmse_nat_pct = (shp_nat['rmse_natural'] - gui_nat['rmse_natural']) / gui_nat['rmse_natural'] * 100
    gpt_r2_nat_pct = (gpt_nat['r2_natural'] - gui_nat['r2_natural']) / gui_nat['r2_natural'] * 100
    gpt_rmse_nat_pct = (gpt_nat['rmse_natural'] - gui_nat['rmse_natural']) / gui_nat['rmse_natural'] * 100
    
    textstr_nat = '\n'.join([
        f"{'Method':<18} {'R²':>12} {'RMSE (m)':>14}",
        '-'*46,
        f"{'ModelMuse GUI':<18} {gui_nat['r2_natural']:>12.4f} {gui_nat['rmse_natural']:>14.4f}",
        f"{'FloPy (Shapefile)':<18} {shp_nat['r2_natural']:.4f}({shp_r2_nat_pct:+.1f}%) {shp_nat['rmse_natural']:.4f}({shp_rmse_nat_pct:+.1f}%)",
        f"{'FloPy (.Gpt)':<18} {gpt_nat['r2_natural']:.4f}({gpt_r2_nat_pct:+.1f}%) {gpt_nat['rmse_natural']:.4f}({gpt_rmse_nat_pct:+.1f}%)",
    ])
    props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray')
    ax1.text(0.98, 0.02, textstr_nat, transform=ax1.transAxes, fontsize=8,
             verticalalignment='bottom', horizontalalignment='right',
             bbox=props, family='monospace')
    
    # ========================================================================
    # Right Panel: Pumping Condition
    # ========================================================================
    ax2 = axes[1]
    ax2.plot([vmin, vmax], [vmin, vmax], 'k--', lw=1.5, label='1:1 Line', alpha=0.7)
    
    for name, data in results.items():
        ax2.scatter(obs_pumping, data['sim_pumping'], 
                   c=colors[name], marker=markers[name], s=100, 
                   edgecolors='white', linewidths=0.5,
                   label=f"{name}", alpha=0.85, zorder=5)
    
    ax2.set_xlabel('Observed Head (m)', fontsize=12)
    ax2.set_ylabel('Simulated Head (m)', fontsize=12)
    ax2.set_title('(b) Pumping Condition (P1 = 100 m³/h)', fontsize=13, fontweight='bold')
    ax2.set_xlim(vmin, vmax)
    ax2.set_ylim(vmin, vmax)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper left', fontsize=9)
    
    # Add metrics text box for Pumping
    gui_pump = results['ModelMuse GUI']
    shp_pump = results['FloPy (Shapefile)']
    gpt_pump = results['FloPy (.Gpt)']
    
    # Calculate percentage changes relative to GUI
    shp_r2_pump_pct = (shp_pump['r2_pumping'] - gui_pump['r2_pumping']) / gui_pump['r2_pumping'] * 100
    shp_rmse_pump_pct = (shp_pump['rmse_pumping'] - gui_pump['rmse_pumping']) / gui_pump['rmse_pumping'] * 100
    gpt_r2_pump_pct = (gpt_pump['r2_pumping'] - gui_pump['r2_pumping']) / gui_pump['r2_pumping'] * 100
    gpt_rmse_pump_pct = (gpt_pump['rmse_pumping'] - gui_pump['rmse_pumping']) / gui_pump['rmse_pumping'] * 100
    
    textstr_pump = '\n'.join([
        f"{'Method':<18} {'R²':>12} {'RMSE (m)':>14}",
        '-'*46,
        f"{'ModelMuse GUI':<18} {gui_pump['r2_pumping']:>12.4f} {gui_pump['rmse_pumping']:>14.4f}",
        f"{'FloPy (Shapefile)':<18} {shp_pump['r2_pumping']:.4f}({shp_r2_pump_pct:+.1f}%) {shp_pump['rmse_pumping']:.4f}({shp_rmse_pump_pct:+.1f}%)",
        f"{'FloPy (.Gpt)':<18} {gpt_pump['r2_pumping']:.4f}({gpt_r2_pump_pct:+.1f}%) {gpt_pump['rmse_pumping']:.4f}({gpt_rmse_pump_pct:+.1f}%)",
    ])
    ax2.text(0.98, 0.02, textstr_pump, transform=ax2.transAxes, fontsize=8,
             verticalalignment='bottom', horizontalalignment='right',
             bbox=props, family='monospace')
    
    # ========================================================================
    # Final Adjustments
    # ========================================================================
    plt.tight_layout()
    
    # Save figure
    plt.savefig(OUTPUT_FIG, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\n✓ Figure saved to: {OUTPUT_FIG}")
    
    plt.show()
    
    return results

# ============================================================================
# Run
# ============================================================================
if __name__ == "__main__":
    results = create_comparison_figure()
