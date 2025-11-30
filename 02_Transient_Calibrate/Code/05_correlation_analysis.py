#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Correlation Analysis for Transient Calibration Results
======================================================
Analyze the correlation between storage parameters (Sy and Ss) and model performance metrics (RMSE, R²).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Configuration
ROOT = Path(__file__).resolve().parent
RESULTS_CSV = ROOT.parent / "Output" / "lhs_runs" / "calibration_results.csv"
OUTPUT_DIR = ROOT.parent / "Correlation"

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def load_calibration_results(csv_path):
    """Load and prepare calibration results."""
    df = pd.read_csv(csv_path)
    
    # Filter only successful runs
    df = df[df['success'] == True].copy()
    
    # Drop non-numeric columns
    df = df.drop(columns=['iteration', 'success'], errors='ignore')
    
    print(f"Loaded {len(df)} successful calibration runs")
    print(f"Parameters: {len([c for c in df.columns if 'sy_' in c or 'ss_' in c])}")
    print(f"Metrics: {len([c for c in df.columns if 'RMSE' in c or 'R2' in c])}")
    
    return df


def calculate_correlation_matrix(df):
    """Calculate correlation matrix for all parameters and metrics."""
    corr_matrix = df.corr()
    return corr_matrix


def extract_parameter_metric_correlations(corr_matrix):
    """Extract correlations between parameters and performance metrics."""
    # Define parameter columns (Sy and Ss)
    param_cols = [c for c in corr_matrix.columns if 'sy_' in c or 'ss_' in c]
    
    # Define metric columns
    metric_cols = [c for c in corr_matrix.columns if 'RMSE' in c or 'R2' in c]
    
    # Extract sub-matrix: parameters vs metrics
    param_metric_corr = corr_matrix.loc[param_cols, metric_cols]
    
    return param_metric_corr


def plot_correlation_heatmap(corr_matrix, output_path, title="Correlation Matrix"):
    """Plot correlation heatmap."""
    plt.figure(figsize=(14, 12))
    
    # Create heatmap
    sns.heatmap(corr_matrix, 
                annot=True, 
                fmt='.3f', 
                cmap='coolwarm', 
                center=0,
                vmin=-1, 
                vmax=1,
                square=True,
                linewidths=0.5,
                cbar_kws={'label': 'Correlation Coefficient'})
    
    plt.title(title, fontsize=14, fontweight='bold')
    plt.xlabel('Metrics', fontsize=12)
    plt.ylabel('Parameters', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Saved: {output_path.name}")


def plot_top_correlations(param_metric_corr, output_path, n_top=15):
    """Plot top N parameter-metric correlations."""
    # Flatten correlation matrix and get top correlations
    corr_flat = param_metric_corr.abs().stack().sort_values(ascending=False)
    top_corr = corr_flat.head(n_top)
    
    # Get signed correlations for the top pairs
    top_pairs = [(idx[0], idx[1], param_metric_corr.loc[idx[0], idx[1]]) 
                 for idx in top_corr.index]
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    y_pos = np.arange(len(top_pairs))
    correlations = [pair[2] for pair in top_pairs]
    labels = [f"{pair[0]} vs {pair[1]}" for pair in top_pairs]
    
    colors = ['#d73027' if c < 0 else '#1a9850' for c in correlations]
    
    ax.barh(y_pos, correlations, color=colors, alpha=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel('Correlation Coefficient', fontsize=12)
    ax.set_title(f'Top {n_top} Parameter-Metric Correlations', fontsize=14, fontweight='bold')
    ax.axvline(x=0, color='black', linewidth=0.8, linestyle='--')
    ax.grid(axis='x', alpha=0.3)
    
    # Add value labels
    for i, v in enumerate(correlations):
        ax.text(v + 0.02 if v > 0 else v - 0.02, i, f'{v:.3f}', 
                va='center', ha='left' if v > 0 else 'right', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Saved: {output_path.name}")


def plot_scatter_matrix(df, output_path):
    """Plot scatter matrix for key parameters vs Overall_R2."""
    # Select key parameters with highest correlation to Overall_R2
    if 'Overall_R2' not in df.columns:
        print("Warning: Overall_R2 not found in data")
        return
    
    param_cols = [c for c in df.columns if 'sy_' in c or 'ss_' in c]
    correlations = df[param_cols].corrwith(df['Overall_R2']).abs().sort_values(ascending=False)
    top_params = correlations.head(6).index.tolist()
    
    # Create scatter matrix
    plot_data = df[top_params + ['Overall_R2']].copy()
    plot_data.columns = [c.replace('sy_', 'Sy_').replace('ss_', 'Ss_') for c in plot_data.columns]
    
    fig = plt.figure(figsize=(14, 12))
    axes = pd.plotting.scatter_matrix(plot_data, 
                                      alpha=0.6, 
                                      figsize=(14, 12), 
                                      diagonal='kde',
                                      c=plot_data['Overall_R2'],
                                      cmap='viridis',
                                      s=30)
    
    # Adjust labels
    for ax in axes.flatten():
        ax.set_xlabel(ax.get_xlabel(), fontsize=9, rotation=45, ha='right')
        ax.set_ylabel(ax.get_ylabel(), fontsize=9, rotation=45, ha='right')
    
    plt.suptitle('Parameter Scatter Matrix (Colored by Overall R²)', 
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Saved: {output_path.name}")


def save_correlation_tables(param_metric_corr, output_dir):
    """Save correlation tables to CSV files."""
    # Sort by absolute correlation with Overall_R2
    if 'Overall_R2' in param_metric_corr.columns:
        sorted_by_overall = param_metric_corr.reindex(
            param_metric_corr['Overall_R2'].abs().sort_values(ascending=False).index
        )
        sorted_by_overall.to_csv(output_dir / "correlations_sorted_by_overall_R2.csv")
        print(f"✓ Saved: correlations_sorted_by_overall_R2.csv")
    
    # Sort by absolute correlation with P1_R2
    if 'P1_R2' in param_metric_corr.columns:
        sorted_by_p1 = param_metric_corr.reindex(
            param_metric_corr['P1_R2'].abs().sort_values(ascending=False).index
        )
        sorted_by_p1.to_csv(output_dir / "correlations_sorted_by_P1_R2.csv")
        print(f"✓ Saved: correlations_sorted_by_P1_R2.csv")
    
    # Sort by absolute correlation with Pz12_R2
    if 'Pz12_R2' in param_metric_corr.columns:
        sorted_by_pz12 = param_metric_corr.reindex(
            param_metric_corr['Pz12_R2'].abs().sort_values(ascending=False).index
        )
        sorted_by_pz12.to_csv(output_dir / "correlations_sorted_by_Pz12_R2.csv")
        print(f"✓ Saved: correlations_sorted_by_Pz12_R2.csv")
    
    # Save full correlation matrix
    param_metric_corr.to_csv(output_dir / "parameter_metric_correlations.csv")
    print(f"✓ Saved: parameter_metric_correlations.csv")


def print_correlation_summary(param_metric_corr):
    """Print summary of key correlations."""
    print("\n" + "="*80)
    print("CORRELATION ANALYSIS SUMMARY")
    print("="*80)
    
    for metric in ['Overall_R2', 'P1_R2', 'Pz12_R2']:
        if metric not in param_metric_corr.columns:
            continue
        
        print(f"\n--- Correlations with {metric} ---")
        sorted_corr = param_metric_corr[metric].abs().sort_values(ascending=False)
        
        print("\nTop 5 Positive Correlations:")
        pos_corr = param_metric_corr[metric].sort_values(ascending=False).head(5)
        for param, corr in pos_corr.items():
            print(f"  {param:25s}: {corr:+.4f}")
        
        print("\nTop 5 Negative Correlations:")
        neg_corr = param_metric_corr[metric].sort_values(ascending=True).head(5)
        for param, corr in neg_corr.items():
            print(f"  {param:25s}: {corr:+.4f}")


def plot_parameter_importance(param_metric_corr, output_path):
    """Plot parameter importance based on average absolute correlation."""
    # Calculate average absolute correlation for each parameter
    param_importance = param_metric_corr.abs().mean(axis=1).sort_values(ascending=False)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    y_pos = np.arange(len(param_importance))
    colors = ['#1f77b4' if 'sy_' in p else '#ff7f0e' for p in param_importance.index]
    
    ax.barh(y_pos, param_importance.values, color=colors, alpha=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([p.replace('sy_', 'Sy ').replace('ss_', 'Ss ') 
                        for p in param_importance.index], fontsize=9)
    ax.set_xlabel('Average Absolute Correlation', fontsize=12)
    ax.set_title('Parameter Importance (Average |Correlation| with All Metrics)', 
                 fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#1f77b4', alpha=0.7, label='Specific Yield (Sy)'),
        Patch(facecolor='#ff7f0e', alpha=0.7, label='Specific Storage (Ss)')
    ]
    ax.legend(handles=legend_elements, loc='lower right')
    
    # Add value labels
    for i, v in enumerate(param_importance.values):
        ax.text(v + 0.01, i, f'{v:.3f}', va='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Saved: {output_path.name}")


def main():
    """Main execution."""
    print("\n" + "="*80)
    print("TRANSIENT CALIBRATION - CORRELATION ANALYSIS")
    print("="*80 + "\n")
    
    # Load data
    print("Loading calibration results...")
    df = load_calibration_results(RESULTS_CSV)
    print()
    
    # Calculate correlations
    print("Calculating correlation matrix...")
    corr_matrix = calculate_correlation_matrix(df)
    param_metric_corr = extract_parameter_metric_correlations(corr_matrix)
    print(f"✓ Correlation matrix: {param_metric_corr.shape}")
    print()
    
    # Save correlation tables
    print("Saving correlation tables...")
    save_correlation_tables(param_metric_corr, OUTPUT_DIR)
    print()
    
    # Generate plots
    print("Generating visualization plots...")
    
    # 1. Full correlation heatmap
    plot_correlation_heatmap(
        param_metric_corr,
        OUTPUT_DIR / "correlation_heatmap.png",
        title="Parameter-Metric Correlation Matrix"
    )
    
    # 2. Top correlations bar chart
    plot_top_correlations(
        param_metric_corr,
        OUTPUT_DIR / "top_correlations.png",
        n_top=15
    )
    
    # 3. Parameter importance
    plot_parameter_importance(
        param_metric_corr,
        OUTPUT_DIR / "parameter_importance.png"
    )
    
    # 4. Scatter matrix
    plot_scatter_matrix(
        df,
        OUTPUT_DIR / "scatter_matrix.png"
    )
    
    print()
    
    # Print summary
    print_correlation_summary(param_metric_corr)
    
    print("\n" + "="*80)
    print(f"✓ All results saved to: {OUTPUT_DIR}")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()
