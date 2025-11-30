#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dewatering Results Visualization
Creates figures for:
1. Steady-state head distribution around pit (Section 7.1.4)
2. Transient drawdown evolution at pit corners (Section 7.2.3)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import flopy
from flopy.utils import HeadFile

# ============================================================================
# Configuration
# ============================================================================
# Get the script directory and project root
SCRIPT_DIR = Path(__file__).parent.resolve()
BASE_DIR = SCRIPT_DIR.parent  # Model for Exam folder
REPORT_DIR = BASE_DIR / "Report"

# Correct paths for head files
MINIMUM_RATE_BHD = BASE_DIR / "Minimum Rate" / "Model" / "Steady_cali.bhd"
MINIMUM_TIME_BHD = BASE_DIR / "Minimum Time" / "Minimum_Time.bhd"

# Pit boundaries (model coordinates)
PIT_X_MIN, PIT_X_MAX = 449, 499  # 50m width
PIT_Y_MIN, PIT_Y_MAX = 480, 530  # 50m height
PIT_BOTTOM = 34.0  # m

# Well positions
WELLS = {
    'W1': (461.5, 530), 'W2': (486.5, 530),  # North
    'W3': (461.5, 480), 'W4': (486.5, 480),  # South
    'W5': (449, 492.5), 'W6': (449, 517.5),  # West
    'W7': (499, 492.5), 'W8': (499, 517.5),  # East
}

# Corner positions for monitoring
CORNERS = {
    'SW': (449, 480),
    'SE': (499, 480),
    'NW': (449, 530),
    'NE': (499, 530),
}

# ============================================================================
# Figure 1: Steady-State Head Distribution
# ============================================================================
def create_steady_state_figure():
    """Create steady-state head distribution figure showing dewatering effect"""
    
    # Load head file from Minimum Rate model
    head_file = MINIMUM_RATE_BHD
    
    if head_file.exists():
        try:
            hds = HeadFile(str(head_file))
            head_data = hds.get_data(totim=hds.get_times()[-1])
            # Use Layer 1 data
            head_2d = head_data[0]
            print(f"Successfully loaded head data from: {head_file}")
            print(f"Head data shape: {head_2d.shape}")
            print(f"Head range: {head_2d[head_2d > 0].min():.2f} to {head_2d[head_2d > 0].max():.2f} m")
        except Exception as e:
            print(f"Could not load head file: {e}, using synthetic data")
            head_2d = create_synthetic_steady_heads()
    else:
        print(f"Head file not found: {head_file}, using synthetic data")
        head_2d = create_synthetic_steady_heads()
    
    # Create figure with two subplots: full domain and zoom
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Determine grid coordinates based on model dimensions
    nrow, ncol = head_2d.shape
    # Assuming model extent: X: 0-1000m, Y: 0-1000m
    x = np.linspace(0, 1000, ncol)
    y = np.linspace(1000, 0, nrow)  # Row 0 is at Y=1000 (north)
    X, Y = np.meshgrid(x, y)
    
    # Mask inactive cells (MODFLOW uses 1e30 for inactive, also mask negatives)
    head_masked = np.ma.masked_where((head_2d < 0) | (head_2d > 50), head_2d)
    print(f"Valid head range: {head_masked.min():.2f} to {head_masked.max():.2f} m")
    
    # ===== Left plot: Full domain =====
    ax1 = axes[0]
    levels = np.arange(33, 39, 0.5)
    cs1 = ax1.contourf(X, Y, head_masked, levels=levels, cmap='Blues_r', extend='both')
    plt.colorbar(cs1, ax=ax1, label='Head (m)', shrink=0.8)
    
    # Pit boundary
    pit_rect1 = plt.Rectangle((PIT_X_MIN, PIT_Y_MIN), 50, 50, 
                               fill=False, edgecolor='red', linewidth=2, 
                               linestyle='-', label='Excavation Pit')
    ax1.add_patch(pit_rect1)
    
    # 34m contour
    cs_34 = ax1.contour(X, Y, head_masked, levels=[34.0], colors='red', linewidths=2)
    
    ax1.set_xlabel('X (m)', fontsize=11)
    ax1.set_ylabel('Y (m)', fontsize=11)
    ax1.set_title('Full Domain Head Distribution', fontsize=11)
    ax1.set_xlim(0, 1000)
    ax1.set_ylim(0, 1000)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    
    # ===== Right plot: Zoomed to pit area =====
    ax2 = axes[1]
    levels_zoom = np.arange(33, 36, 0.2)
    cs2 = ax2.contourf(X, Y, head_masked, levels=levels_zoom, cmap='Blues_r', extend='both')
    cbar2 = plt.colorbar(cs2, ax=ax2, label='Hydraulic Head (m)', shrink=0.8)
    
    # Contour lines
    cs_lines = ax2.contour(X, Y, head_masked, levels=np.arange(33, 36, 0.5), 
                           colors='darkblue', linewidths=0.8, alpha=0.7)
    ax2.clabel(cs_lines, inline=True, fontsize=9, fmt='%.1f')
    
    # 34m target contour
    cs_34_zoom = ax2.contour(X, Y, head_masked, levels=[34.0], colors='red', linewidths=2.5)
    
    # Pit boundary
    pit_rect2 = plt.Rectangle((PIT_X_MIN, PIT_Y_MIN), 50, 50, 
                               fill=False, edgecolor='red', linewidth=2.5, 
                               linestyle='--', label='Pit Boundary')
    ax2.add_patch(pit_rect2)
    
    # Plot wells with labels
    for well_name, (wx, wy) in WELLS.items():
        ax2.plot(wx, wy, 'ko', markersize=10, markerfacecolor='yellow', 
                markeredgewidth=2, zorder=10)
        ax2.annotate(well_name, (wx, wy), textcoords="offset points", 
                    xytext=(6, 6), fontsize=8, fontweight='bold')
    
    # Mark pit corners and show head values
    corner_cells = {
        'SW': {'row': 52, 'col': 45, 'pos': (449, 480)},
        'SE': {'row': 52, 'col': 50, 'pos': (499, 480)},
        'NW': {'row': 47, 'col': 45, 'pos': (449, 530)},
        'NE': {'row': 47, 'col': 50, 'pos': (499, 530)},
    }
    for corner, info in corner_cells.items():
        h_val = head_2d[info['row'], info['col']]
        ax2.plot(*info['pos'], 'rs', markersize=8, markerfacecolor='none', 
                markeredgewidth=2, zorder=11)
        ax2.annotate(f"{corner}: {h_val:.2f}m", info['pos'], 
                    textcoords="offset points", xytext=(-25, -15), 
                    fontsize=8, color='darkred', fontweight='bold')
    
    ax2.set_xlabel('X Coordinate (m)', fontsize=11)
    ax2.set_ylabel('Y Coordinate (m)', fontsize=11)
    ax2.set_title('Zoomed View: Pit Area\n(8 wells @ 0.0103 m³/s each)', fontsize=11)
    ax2.set_xlim(400, 550)
    ax2.set_ylim(430, 580)
    ax2.set_aspect('equal')
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(True, alpha=0.3)
    
    # Overall title
    fig.suptitle('Steady-State Dewatering: Minimum Pumping Rate Analysis', 
                 fontsize=13, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    
    output_path = REPORT_DIR / "dewatering_steady_state.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"✓ Steady-state dewatering figure saved to: {output_path}")
    
    plt.show()
    return output_path


def create_synthetic_steady_heads():
    """Create synthetic head data for demonstration"""
    nrow, ncol = 127, 102
    head = np.zeros((nrow, ncol))
    
    # Base gradient from south (high) to north (low) and west to east
    for i in range(nrow):
        for j in range(ncol):
            x = j * 10  # approximate
            y = (nrow - 1 - i) * 10
            # Natural head ~37-38
            head[i, j] = 37.5 - 0.0005 * x + 0.0003 * (y - 500)
    
    # Pit center coordinates: X=474, Y=505 (center of 449-499, 480-530)
    # Convert to row/col: row = (1000 - Y) / 10, col = X / 10
    pit_center_x = 474  # center of pit (449+499)/2
    pit_center_y = 505  # center of pit (480+530)/2
    pit_center_row = int((1000 - pit_center_y) / 10)  # ~49.5 -> 50
    pit_center_col = int(pit_center_x / 10)  # ~47
    
    # Add dewatering drawdown from 8 wells around the pit perimeter
    well_positions = [
        (461.5, 530), (486.5, 530),  # North side
        (461.5, 480), (486.5, 480),  # South side
        (449, 492.5), (449, 517.5),  # West side
        (499, 492.5), (499, 517.5),  # East side
    ]
    
    for i in range(nrow):
        for j in range(ncol):
            cell_x = j * 10
            cell_y = (nrow - 1 - i) * 10
            
            total_drawdown = 0
            for wx, wy in well_positions:
                dist = np.sqrt((cell_x - wx)**2 + (cell_y - wy)**2)
                if dist < 5:
                    dist = 5  # Avoid division issues near well
                # Thiem-like drawdown from each well
                drawdown = 0.6 * np.log(200 / max(dist, 5)) if dist < 200 else 0
                total_drawdown += max(0, drawdown)
            
            head[i, j] -= total_drawdown
    
    return head


# ============================================================================
# Figure 2: Transient Drawdown Evolution
# ============================================================================
def create_transient_figure():
    """Create transient drawdown evolution figure at pit corners using real model data"""
    
    # Try to load actual transient head file
    head_file = MINIMUM_TIME_BHD
    
    # Pit corner cell indices
    corner_cells = {
        'SW Corner': {'row': 52, 'col': 45},
        'SE Corner': {'row': 52, 'col': 50},
        'NW Corner': {'row': 47, 'col': 45},
        'NE Corner': {'row': 47, 'col': 50},
    }
    
    heads = {}
    times_days = None
    crossing_times = {}  # Store when each corner crosses 34m
    
    if head_file.exists():
        try:
            hds = HeadFile(str(head_file))
            times = np.array(hds.get_times())
            times_days = times / 86400.0  # Convert seconds to days
            
            print(f"Successfully loaded transient head data from: {head_file}")
            print(f"Number of time steps: {len(times)}, Total time: {times_days[-1]:.1f} days")
            
            for corner, cell_info in corner_cells.items():
                row, col = cell_info['row'], cell_info['col']
                head_series = []
                for t in times:
                    head_data = hds.get_data(totim=t)
                    head_series.append(head_data[0, row, col])
                heads[corner] = np.array(head_series)
                
                # Find crossing time via linear interpolation
                h_arr = np.array(head_series)
                for i in range(1, len(h_arr)):
                    if h_arr[i-1] >= 34 and h_arr[i] < 34:
                        t_cross = times_days[i-1] + (34 - h_arr[i-1]) / (h_arr[i] - h_arr[i-1]) * (times_days[i] - times_days[i-1])
                        crossing_times[corner] = t_cross
                        print(f"  {corner}: crosses 34m at {t_cross:.1f} days")
                        break
                
        except Exception as e:
            print(f"Could not load transient head file: {e}")
            print("Using synthetic data instead")
            times_days, heads = create_synthetic_transient_data()
            crossing_times = {'SW Corner': 15, 'SE Corner': 11, 'NW Corner': 35, 'NE Corner': 23}
    else:
        print(f"Transient head file not found: {head_file}")
        print("Using synthetic data instead")
        times_days, heads = create_synthetic_transient_data()
        crossing_times = {'SW Corner': 15, 'SE Corner': 11, 'NW Corner': 35, 'NE Corner': 23}
    
    # Get the maximum crossing time (slowest corner)
    max_cross_time = max(crossing_times.values()) if crossing_times else 35
    min_cross_time = min(crossing_times.values()) if crossing_times else 11
    
    # Create figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    
    colors = {'SW Corner': '#2ecc71', 'SE Corner': '#3498db', 
              'NW Corner': '#9b59b6', 'NE Corner': '#e74c3c'}
    markers = {'SW Corner': 'o', 'SE Corner': 's', 'NW Corner': '^', 'NE Corner': 'D'}
    
    # ===== Left plot: Full time range =====
    ax1 = axes[0]
    for corner, head_vals in heads.items():
        ax1.plot(times_days, head_vals, marker=markers[corner], markersize=4,
                color=colors[corner], linewidth=1.5, label=corner, alpha=0.9,
                markevery=max(1, len(times_days)//12))
    
    ax1.axhline(y=34.0, color='red', linestyle='--', linewidth=2, label='Target (34 m)')
    ax1.fill_between([0, max(times_days)], [34, 34], [33, 33], alpha=0.15, color='green')
    
    # Mark crossing times with vertical lines
    for corner, t_cross in crossing_times.items():
        ax1.axvline(x=t_cross, color=colors[corner], linestyle=':', alpha=0.5, linewidth=1)
    
    ax1.set_xlabel('Time (days)', fontsize=11)
    ax1.set_ylabel('Hydraulic Head (m)', fontsize=11)
    ax1.set_title('Full Simulation Period (0-170 days)', fontsize=11)
    ax1.set_xlim(-2, 175)
    ax1.set_ylim(33, 38)
    ax1.legend(loc='upper right', fontsize=8)
    ax1.grid(True, alpha=0.3)
    
    # ===== Right plot: Zoom on critical period =====
    ax2 = axes[1]
    
    # Filter data for zoomed view (0-50 days)
    zoom_mask = times_days <= 50
    t_zoom = times_days[zoom_mask]
    
    for corner, head_vals in heads.items():
        h_zoom = head_vals[zoom_mask]
        ax2.plot(t_zoom, h_zoom, marker=markers[corner], markersize=6,
                color=colors[corner], linewidth=2, label=corner, alpha=0.9,
                markevery=max(1, len(t_zoom)//10))
        
        # Mark the crossing point
        if corner in crossing_times:
            t_c = crossing_times[corner]
            if t_c <= 50:
                ax2.axvline(x=t_c, color=colors[corner], linestyle=':', alpha=0.7, linewidth=1.5)
                ax2.annotate(f'{t_c:.0f}d', (t_c, 34.05), fontsize=8, 
                            color=colors[corner], ha='center', fontweight='bold')
    
    ax2.axhline(y=34.0, color='red', linestyle='--', linewidth=2.5, label='Target (34 m)')
    ax2.fill_between([0, 50], [34, 34], [33.3, 33.3], alpha=0.2, color='green',
                    label='Dewatered Zone')
    
    ax2.set_xlabel('Time (days)', fontsize=11)
    ax2.set_ylabel('Hydraulic Head (m)', fontsize=11)
    ax2.set_title('Zoomed View: Critical Dewatering Period', fontsize=11)
    ax2.set_xlim(-1, 52)
    ax2.set_ylim(33.3, 37.5)
    ax2.legend(loc='upper right', fontsize=8)
    ax2.grid(True, alpha=0.3)
    
    # Add annotation with correct times
    ax2.annotate(f'SE: fastest (~{crossing_times.get("SE Corner", 11):.0f} days)', 
                xy=(12, 33.6), fontsize=9, color='#3498db', style='italic')
    ax2.annotate(f'NW: slowest (~{crossing_times.get("NW Corner", 35):.0f} days)', 
                xy=(35, 34.1), fontsize=9, color='#9b59b6', style='italic')
    
    # Overall title with correct time
    fig.suptitle(f'Transient Dewatering: Minimum Time Analysis\n'
                 f'(All corners < 34 m by day {max_cross_time:.0f}; 8 wells @ 0.0103 m³/s)', 
                 fontsize=12, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    
    output_path = REPORT_DIR / "dewatering_transient.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"✓ Transient dewatering figure saved to: {output_path}")
    print(f"\nSummary: Dewatering complete at all corners by day {max_cross_time:.0f}")
    print(f"  Fastest: SE Corner at {min_cross_time:.0f} days")
    print(f"  Slowest: NW Corner at {max_cross_time:.0f} days")
    
    plt.show()
    return output_path


def create_synthetic_transient_data():
    """Create synthetic transient data for demonstration"""
    times = np.array([0, 0.5, 1, 2, 5, 10, 15, 20, 30, 60, 120, 170])
    
    def drawdown_curve(t, rate=0.25, asymptote=33.5):
        return asymptote + (37.17 - asymptote) * np.exp(-rate * t)
    
    heads = {
        'SW Corner': drawdown_curve(times, rate=0.22, asymptote=33.48),
        'SE Corner': drawdown_curve(times, rate=0.21, asymptote=33.50),
        'NW Corner': drawdown_curve(times, rate=0.22, asymptote=33.49),
        'NE Corner': drawdown_curve(times, rate=0.18, asymptote=33.54),
    }
    
    return times, heads


# ============================================================================
# Main
# ============================================================================
if __name__ == "__main__":
    print("Creating dewatering visualization figures...")
    print("=" * 50)
    
    print("\n1. Creating steady-state head distribution figure...")
    create_steady_state_figure()
    
    print("\n2. Creating transient drawdown evolution figure...")
    create_transient_figure()
    
    print("\n✓ All figures created successfully!")
