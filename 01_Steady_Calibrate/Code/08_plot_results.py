# -*- coding: utf-8 -*-
"""
Script to plot the simulated heads of a modflow 6 model (steady state)

Things to be adjusted by yourself:
    directory_S
    direcory_O
    model
    csv
    
    condition
    error
"""

import matplotlib
matplotlib.use('TkAgg')  # Use TkAgg backend to avoid Qt font warnings
import matplotlib.pyplot as plt
import flopy
import math
import pandas as pd
import numpy as np
import os
from pathlib import Path

plt.rcParams['figure.dpi'] = 200

# Get the script directory and project root
SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = SCRIPT_DIR.parent.parent  # Model for Exam folder

# path to the directory where your results are (relative paths)
directory_S = PROJECT_ROOT / "01_Steady_Calibrate" / "Out" / "run_modelmuse"
# Alternative paths:
# directory_S = PROJECT_ROOT / "02_Transient_Calibrate" / "Output" / "Steady"
# directory_S = PROJECT_ROOT / "01_Steady_Calibrate" / "Final Model"

# path to the directory where the observed heads are
directory_O = PROJECT_ROOT / "Instructions"

# name of your model
model = 'Steady_cali'

# name of observed heads file (the one from Ufora)
csv = 'steady_state_calibration'

# open the simulation results
sim_out = flopy.utils.observationfile.Mf6Obs(os.path.join(directory_S, '{}.ob_gw_out_head.csv'.format(model)), isBinary=False)

# customize dataframe
sim_df = sim_out.get_dataframe().transpose()
sim_df.columns = ['natural', 'pumping']
sim_df = sim_df.iloc[1:]

index_labels = ["P1", "Pz2", "Pz3", "Pz4", "Pz5", "Pz6", "Pz7", "Pz8", "Pz9", "Pz10", "Pz11", "Pz12"] #CHANGE THIS SO IT MATCHES THE ORDER OF YOUR SIMULATED FILE

sim_df.index=index_labels

#open the file with the observed heads
obs_df = pd.read_csv(os.path.join(directory_O, '{}.csv'.format(csv)),delimiter=';')
obs_df = obs_df.set_index('Unnamed: 0')

#calculate residual between simulated and observed heads
sim_df['natural_R'] = sim_df['natural'] - obs_df['natural']
sim_df['pumping_R'] = sim_df['pumping'] - obs_df['pumping']
#%%
# Select the condition to analyze: pumping/natural
condition = 'pumping'

# Calculate various evaluation metrics
obs_df = obs_df.loc[sim_df.index]
observed = obs_df[condition].values
simulated = sim_df[condition].values
residuals = sim_df['{}_R'.format(condition)].values

# 1. RMSE (Root Mean Square Error)
rmse = np.sqrt(np.mean(residuals**2))

# 2. MAE (Mean Absolute Error)
mae = np.mean(np.abs(residuals))

# 3. R² (Coefficient of Determination)
ss_res = np.sum(residuals**2)  # Residual sum of squares
ss_tot = np.sum((observed - np.mean(observed))**2)  # Total sum of squares
r2 = 1 - (ss_res / ss_tot)

# 4. Max residual
max_residual = np.max(np.abs(residuals))

# Print evaluation metrics
print("\n" + "="*60)
print(f"Model Evaluation Metrics - {condition.upper()} Condition")
print("="*60)
print(f"RMSE (Root Mean Square Error):    {rmse:.4f} m")
print(f"MAE (Mean Absolute Error):         {mae:.4f} m")
print(f"R² (Coefficient of Determination): {r2:.4f}")
print(f"Max Residual:                      {max_residual:.4f} m")
print(f"Number of Observation Points:      {len(observed)}")
print("="*60)

# ===========================================================
# Combined output: 4 subplots
# Top-left: pumping parity plot (RMSE + R²)
# Top-right: natural parity plot (RMSE + R²)
# Bottom-left: natural line plot
# Bottom-right: pumping line plot
# ===========================================================

# Ensure same station order as sim_df
labels = sim_df.index.tolist()
x = np.arange(len(labels))

# Data for line plots
obs_nat = obs_df['natural'].values
sim_nat = sim_df['natural'].values
obs_pmp = obs_df['pumping'].values
sim_pmp = sim_df['pumping'].values

# Pumping scatter data
resid_pmp = sim_df['pumping_R'].values
min_val_pmp = min(obs_pmp.min(), sim_pmp.min()) - 0.3
max_val_pmp = max(obs_pmp.max(), sim_pmp.max()) + 0.3

# Natural scatter data
resid_nat = sim_df['natural_R'].values
min_val_nat = min(obs_nat.min(), sim_nat.min()) - 0.3
max_val_nat = max(obs_nat.max(), sim_nat.max()) + 0.3

# Calculate pumping evaluation metrics
rmse_pmp = np.sqrt(np.mean(resid_pmp**2))
mae_pmp = np.mean(np.abs(resid_pmp))
ss_res_pmp = np.sum(resid_pmp**2)
ss_tot_pmp = np.sum((obs_pmp - np.mean(obs_pmp))**2)
r2_pmp = 1 - (ss_res_pmp / ss_tot_pmp)

# Calculate natural evaluation metrics
rmse_nat = np.sqrt(np.mean(resid_nat**2))
mae_nat = np.mean(np.abs(resid_nat))
ss_res_nat = np.sum(resid_nat**2)
ss_tot_nat = np.sum((obs_nat - np.mean(obs_nat))**2)
r2_nat = 1 - (ss_res_nat / ss_tot_nat)

# Print evaluation metrics for NATURAL condition
print("\n" + "="*60)
print(f"Model Evaluation Metrics - NATURAL Condition")
print("="*60)
print(f"RMSE (Root Mean Square Error):    {rmse_nat:.4f} m")
print(f"MAE (Mean Absolute Error):         {mae_nat:.4f} m")
print(f"R² (Coefficient of Determination): {r2_nat:.4f}")
print(f"Max Residual:                      {np.max(np.abs(resid_nat)):.4f} m")
print(f"Number of Observation Points:      {len(obs_nat)}")
print("="*60)

print("\n" + "="*60)
print(f"Model Evaluation Metrics - PUMPING Condition")
print("="*60)
print(f"RMSE (Root Mean Square Error):    {rmse_pmp:.4f} m")
print(f"MAE (Mean Absolute Error):         {mae_pmp:.4f} m")
print(f"R² (Coefficient of Determination): {r2_pmp:.4f}")
print(f"Max Residual:                      {np.max(np.abs(resid_pmp)):.4f} m")
print(f"Number of Observation Points:      {len(obs_pmp)}")
print("="*60)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))

# ---------- Subplot 1: natural line plot ----------
ax1.plot(x, obs_nat, marker='o', markersize=5, linewidth=1.5, label='Observed', color='#1f77b4')
ax1.plot(x, sim_nat, marker='s', markersize=5, linewidth=1.5, label='Simulated', color='#ff7f0e')
ax1.set_title('Heads Comparison - Natural', fontsize=11, fontweight='bold')
ax1.set_xlabel('Observation Points', fontsize=10)
ax1.set_ylabel('Head (m)', fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.legend(loc='best', fontsize=9)
ax1.set_xticks(x)
ax1.set_xticklabels(labels, rotation=45, fontsize=8, ha='right')

# Add R² and RMSE
textstr_nat = f'R² = {r2_nat:.3f}\nRMSE = {rmse_nat:.3f} m'
ax1.text(0.02, 0.98, textstr_nat, transform=ax1.transAxes, fontsize=9,
         va='top', ha='left',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.75))

# ---------- Subplot 2: pumping line plot ----------
ax2.plot(x, obs_pmp, marker='o', markersize=5, linewidth=1.5, label='Observed', color='#1f77b4')
ax2.plot(x, sim_pmp, marker='s', markersize=5, linewidth=1.5, label='Simulated', color='#ff7f0e')
ax2.set_title('Heads Comparison - Pumping', fontsize=11, fontweight='bold')
ax2.set_xlabel('Observation Points', fontsize=10)
ax2.set_ylabel('Head (m)', fontsize=10)
ax2.grid(True, alpha=0.3)
ax2.legend(loc='best', fontsize=9)
ax2.set_xticks(x)
ax2.set_xticklabels(labels, rotation=45, fontsize=8, ha='right')

# Add R² and RMSE
textstr_pmp = f'R² = {r2_pmp:.3f}\nRMSE = {rmse_pmp:.3f} m'
ax2.text(0.02, 0.98, textstr_pmp, transform=ax2.transAxes, fontsize=9,
         va='top', ha='left',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.75))

plt.tight_layout()
plt.show()
