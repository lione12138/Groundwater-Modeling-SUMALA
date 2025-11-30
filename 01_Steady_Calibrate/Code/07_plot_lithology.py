#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lithology and Hydraulic Conductivity Distribution Map
Creates a two-panel figure showing KX distribution for Layer 1 and Layer 2
Only shows active model cells based on IDOMAIN

Author: Generated for Groundwater Modeling Report
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
import numpy as np
from pathlib import Path
from flopy.mf6.utils import MfGrdFile

# ============================================================================
# Configuration
# ============================================================================
# Get the script directory and project root
SCRIPT_DIR = Path(__file__).parent.resolve()
BASE_DIR = SCRIPT_DIR.parent  # Steady_Calibrate folder
SHP_PATH = BASE_DIR / "Lithology" / "KX-new.shp"
DIS_GRB = BASE_DIR / "Model" / "Steady_cali.dis.grb"
OUTPUT_DIR = BASE_DIR / "Out"
OUTPUT_DIR.mkdir(exist_ok=True)

# Color scheme matching the reference figure style
LITHOLOGY_COLORS = {
    # Layer 1 (Alluvial deposits)
    'Clayey_Sand': '#808080',      # Gray
    'Loamy_Sand': '#9ACD32',       # Yellow-green
    'Clean_Gravel': '#CD853F',     # Brown/orange
    'Sandy_Gravel': '#ADD8E6',     # Light blue (background)
    # Layer 2 (Bedrock)
    'Sand_Stone': '#9ACD32',       # Yellow-green
    'Lime_Stone': '#4169E1',       # Blue (background)
}

# ============================================================================
# Main Plotting Function
# ============================================================================
def create_lithology_figure():
    """Create two-panel lithology map for Layer 1 and Layer 2 using active cells only"""
    
    # Read shapefile
    gdf = gpd.read_file(str(SHP_PATH))
    print(f"Loaded {len(gdf)} polygons from shapefile")
    print(gdf[['OBJECTNAME', 'KX', 'HIGH_Z', 'LOW_Z']])
    
    # Read grid and IDOMAIN from dis.grb
    grb = MfGrdFile(str(DIS_GRB))
    nlay, nrow, ncol = grb.nlay, grb.nrow, grb.ncol
    idomain = grb.idomain.reshape((nlay, nrow, ncol))
    
    # Get cell centers
    mg = grb.modelgrid
    xcenters = mg.xcellcenters
    ycenters = mg.ycellcenters
    
    print(f"\nGrid: {nlay} layers, {nrow} rows, {ncol} cols")
    print(f"Layer 1 active cells: {np.sum(idomain[0] == 1)}")
    print(f"Layer 2 active cells: {np.sum(idomain[1] == 1)}")
    
    # Separate shapefile by layer
    layer1 = gdf[gdf['HIGH_Z'] == 'Model_Top'].copy()
    layer2 = gdf[gdf['HIGH_Z'] == 'Upper_Aquifer_Bottom'].copy()
    
    print(f"\nLayer 1 lithologies: {layer1['OBJECTNAME'].tolist()}")
    print(f"Layer 2 lithologies: {layer2['OBJECTNAME'].tolist()}")
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))
    
    # ===== Panel a): Layer 1 - Alluvial Deposits =====
    ax1 = axes[0]
    ax1.set_title('a) Layer 1: Alluvial Deposits', fontsize=12, fontweight='bold', loc='left')
    
    # Create active domain mask for Layer 1
    active_l1 = idomain[0] == 1
    
    # Draw background (Sandy Gravel) only for active cells
    # Use pcolormesh with masked array
    background_l1 = np.ones_like(active_l1, dtype=float)
    background_l1[~active_l1] = np.nan
    
    # Get cell edges for pcolormesh
    dx = xcenters[0, 1] - xcenters[0, 0] if ncol > 1 else 10
    dy = ycenters[0, 0] - ycenters[1, 0] if nrow > 1 else 10
    
    x_edges = np.append(xcenters[0, :] - dx/2, xcenters[0, -1] + dx/2)
    y_edges = np.append(ycenters[:, 0] + dy/2, ycenters[-1, 0] - dy/2)
    
    # Plot background for active cells
    ax1.pcolormesh(x_edges, y_edges, background_l1, 
                   cmap=ListedColormap([LITHOLOGY_COLORS['Sandy_Gravel']]),
                   shading='flat', zorder=0)
    
    # Clip lithology polygons to active domain and plot
    for idx, row in layer1.iterrows():
        lith_name = row['OBJECTNAME']
        color = LITHOLOGY_COLORS.get(lith_name, '#CCCCCC')
        layer1[layer1.index == idx].plot(ax=ax1, facecolor=color, edgecolor='black', 
                                          linewidth=0.8, zorder=1)
    
    # Draw active domain boundary
    # Find boundary of active cells
    from scipy import ndimage
    boundary = active_l1.astype(int) - ndimage.binary_erosion(active_l1).astype(int)
    boundary_rows, boundary_cols = np.where(boundary)
    
    # Set axis properties
    ax1.set_xlim(0, 1000)
    ax1.set_ylim(0, 1000)
    ax1.set_aspect('equal')
    ax1.set_xlabel('X (m)', fontsize=11)
    ax1.set_ylabel('Y (m)', fontsize=11)
    
    # North arrow
    ax1.annotate('N', xy=(0.93, 0.93), xycoords='axes fraction', fontsize=14, 
                ha='center', fontweight='bold')
    ax1.annotate('', xy=(0.93, 0.90), xytext=(0.93, 0.82),
                xycoords='axes fraction', 
                arrowprops=dict(arrowstyle='->', lw=2))
    
    # Legend for Layer 1 (only lithology, no hill/river)
    legend_elements_l1 = [
        mpatches.Patch(facecolor=LITHOLOGY_COLORS['Clayey_Sand'], edgecolor='black', 
                      label='Clayey Sand'),
        mpatches.Patch(facecolor=LITHOLOGY_COLORS['Loamy_Sand'], edgecolor='black', 
                      label='Loamy Sand'),
        mpatches.Patch(facecolor=LITHOLOGY_COLORS['Clean_Gravel'], edgecolor='black', 
                      label='Clean Gravel'),
        mpatches.Patch(facecolor=LITHOLOGY_COLORS['Sandy_Gravel'], edgecolor='black', 
                      label='Sandy Gravel'),
    ]
    ax1.legend(handles=legend_elements_l1, loc='lower left', fontsize=9, 
              framealpha=0.95, edgecolor='black')
    
    # Scale bar
    ax1.plot([760, 960], [60, 60], 'k-', linewidth=3)
    ax1.text(860, 25, '200 m', ha='center', fontsize=10, fontweight='bold')
    
    # ===== Panel b): Layer 2 - Bedrock =====
    ax2 = axes[1]
    ax2.set_title('b) Layer 2: Bedrock', fontsize=12, fontweight='bold', loc='left')
    
    # Create active domain mask for Layer 2
    active_l2 = idomain[1] == 1
    
    # Draw background (Limestone) only for active cells
    background_l2 = np.ones_like(active_l2, dtype=float)
    background_l2[~active_l2] = np.nan
    
    ax2.pcolormesh(x_edges, y_edges, background_l2, 
                   cmap=ListedColormap([LITHOLOGY_COLORS['Lime_Stone']]),
                   shading='flat', zorder=0)
    
    # Plot Layer 2 polygons
    for idx, row in layer2.iterrows():
        lith_name = row['OBJECTNAME']
        color = LITHOLOGY_COLORS.get(lith_name, '#CCCCCC')
        layer2[layer2.index == idx].plot(ax=ax2, facecolor=color, edgecolor='black', 
                                          linewidth=0.8, zorder=1)
    
    # Set axis properties
    ax2.set_xlim(0, 1000)
    ax2.set_ylim(0, 1000)
    ax2.set_aspect('equal')
    ax2.set_xlabel('X (m)', fontsize=11)
    ax2.set_ylabel('Y (m)', fontsize=11)
    
    # North arrow
    ax2.annotate('N', xy=(0.93, 0.93), xycoords='axes fraction', fontsize=14, 
                ha='center', fontweight='bold')
    ax2.annotate('', xy=(0.93, 0.90), xytext=(0.93, 0.82),
                xycoords='axes fraction', 
                arrowprops=dict(arrowstyle='->', lw=2))
    
    # Legend for Layer 2 (only lithology)
    legend_elements_l2 = [
        mpatches.Patch(facecolor=LITHOLOGY_COLORS['Sand_Stone'], edgecolor='black', 
                      label='Sandstone'),
        mpatches.Patch(facecolor=LITHOLOGY_COLORS['Lime_Stone'], edgecolor='black', 
                      label='Limestone'),
    ]
    ax2.legend(handles=legend_elements_l2, loc='lower left', fontsize=9, 
              framealpha=0.95, edgecolor='black')
    
    # Scale bar
    ax2.plot([760, 960], [60, 60], 'k-', linewidth=3)
    ax2.text(860, 25, '200 m', ha='center', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    
    # Save figure
    output_path = OUTPUT_DIR / "lithology_kx_distribution.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\n✓ Figure saved to: {output_path}")
    
    plt.show()
    return output_path


# ============================================================================
# Main
# ============================================================================
if __name__ == "__main__":
    print("Creating lithology and hydraulic conductivity distribution map...")
    print("=" * 60)
    create_lithology_figure()
    print("\n✓ Done!")
