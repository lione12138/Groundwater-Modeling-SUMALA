# -*- coding: utf-8 -*-
"""
Script to plot the simulated drawdown of a modflow 6 model (transient)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import flopy
from pathlib import Path

plt.rcParams['figure.dpi'] = 300

# Get the script directory and project root
SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = SCRIPT_DIR.parent.parent  # Model for Exam folder

# ===== Paths and Names =====
directory_S_TS = PROJECT_ROOT / "Output"
directory_O    = PROJECT_ROOT / "Instructions"
model_TS       = 'transient'
csv            = 'calibration_transient'

# Read observed drawdown
obs_ts = pd.read_csv(os.path.join(directory_O, f'{csv}.csv'),
                     index_col='time (s)', delimiter=';')
# Ensure observation time is numeric
obs_ts.index = pd.to_numeric(obs_ts.index, errors='coerce')

# ---------- Read simulation results: use get_data() to avoid datetime overflow ----------
simTS_out = flopy.utils.observationfile.Mf6Obs(
    os.path.join(directory_S_TS, f'{model_TS}.ob_gw_out_head.csv'),
    isBinary=False
)

# get_data returns structured array, convert directly to DataFrame with column names
sim_df = pd.DataFrame(simTS_out.get_data())

# Find time column name (compatible with 'time', 'totim', 'time (s)', etc.)
time_candidates = [c for c in sim_df.columns if c.strip().lower() in ['time', 'totim', 'time (s)']]
if not time_candidates:
    raise RuntimeError(f"Time column not found. Available columns: {list(sim_df.columns)}")
time_col = time_candidates[0]

# Set seconds as time index and sort
sim_df[time_col] = pd.to_numeric(sim_df[time_col], errors='coerce')
sim_df = sim_df.dropna(subset=[time_col]).set_index(time_col).sort_index()

# ---- Smart matching for P1 / Pz12 columns (compatible with possible suffixes) ----
def pick_col(df, target):
    t = target.lower()
    # Exact match
    for c in df.columns:
        if c.strip().lower() == t:
            return c
    # Fallback: prefix or contains match
    for c in df.columns:
        base = c.split('_')[0].strip().lower()
        if base == t or t in c.strip().lower():
            return c
    raise KeyError(f"Cannot find observation column: {target}; Available columns: {list(df.columns)}")

col_p1   = pick_col(sim_df, 'P1')
col_pz12 = pick_col(sim_df, 'Pz12')

# ---- Calculate drawdown (baseline selection: if first time is 0 s, use second row; otherwise use first row) ----
times = sim_df.index.values.astype(float)
baseline_ix = 1 if (len(times) > 1 and abs(times[0]) < 1e-9) else 0
baseline = sim_df.iloc[baseline_ix]

sim_df['P1_DD']   = baseline[col_p1]   - sim_df[col_p1]
sim_df['Pz12_DD'] = baseline[col_pz12] - sim_df[col_pz12]

# ---- Interpolate simulation values to observed time points (consistent with your existing plotting code) ----
obs_ts.index = pd.to_numeric(obs_ts.index, errors='coerce')
obs_times = obs_ts.index.values.astype(float)
sim_times = sim_df.index.values.astype(float)

import numpy as np
p1_sim_interp   = np.interp(obs_times, sim_times, sim_df['P1_DD'].values)
pz12_sim_interp = np.interp(obs_times, sim_times, sim_df['Pz12_DD'].values)


# ---- Plotting (same time axis) ----
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(obs_times, obs_ts['P1'].values,            label='Observed P1')
ax.plot(obs_times, p1_sim_interp,                  label='Simulated P1 (interp)')
ax.plot(obs_times, obs_ts['Pz12'].values,          label='Observed Pz12')
ax.plot(obs_times, pz12_sim_interp,                label='Simulated Pz12 (interp)')
ax.set_title('Observed vs Simulated Drawdown')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Drawdown (m)')
ax.grid(alpha=0.3)
ax.legend()
plt.tight_layout()
plt.show()
