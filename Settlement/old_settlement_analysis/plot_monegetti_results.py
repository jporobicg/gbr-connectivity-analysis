#!/usr/bin/env python3
"""
Create plots of Monegetti piecewise model fits for each family.
Shows data points and fitted model curves.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
from monegetti_piecewise_model import monegetti_competency
from improved_competency_model import (
    load_all_data,
    logistic_competency,
    gompertz_competency,
    weibull_competency
)

sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 10)

def load_latest_estimations(estimations_dir='estimations'):
    """Load the latest estimation file for each family."""
    estimation_files = glob.glob(os.path.join(estimations_dir, 'estimations_*.csv'))
    
    families_data = {}
    for file in estimation_files:
        df = pd.read_csv(file)
        family = df['Family'].iloc[0]
        # Keep only the latest file for each family
        if family not in families_data:
            families_data[family] = df
        else:
            # Compare timestamps
            old_ts = families_data[family]['Timestamp'].iloc[0]
            new_ts = df['Timestamp'].iloc[0]
            if new_ts > old_ts:
                families_data[family] = df
    
    return families_data

def plot_family_monegetti(family, df_family, estimation_data, output_dir='figures'):
    """Create plot for a single family showing data and fitted models."""
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Get Monegetti parameters
    est = estimation_data.iloc[0]
    a, b1, v1, b2, tc, Tcp = est['a'], est['b1'], est['v1'], est['b2'], est['tc'], est['Tcp']
    
    # Aggregate data by age
    age_summary = df_family.groupby('LarvalAge').agg({
        'NoSet': 'sum',
        'NoAlive': 'sum'
    }).reset_index()
    age_summary['PropSettled'] = age_summary['NoSet'] / age_summary['NoAlive']
    age_summary = age_summary.sort_values('LarvalAge')
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Plot data points with error bars (binomial confidence)
    ages = age_summary['LarvalAge'].values
    props = age_summary['PropSettled'].values
    n_total = age_summary['NoAlive'].values
    n_settled = age_summary['NoSet'].values
    
    # Calculate 95% CI using normal approximation
    se = np.sqrt(props * (1 - props) / n_total)
    ci_lower = np.maximum(0, props - 1.96 * se)
    ci_upper = np.minimum(1, props + 1.96 * se)
    
    # Plot data with error bars
    ax.errorbar(ages, props, yerr=[props - ci_lower, ci_upper - props],
                fmt='o', color='#2c3e50', markersize=8, capsize=5,
                capthick=2, label='Observed data', zorder=3, alpha=0.8)
    
    # Plot fitted Monegetti model
    age_range = np.linspace(0, ages.max() + 5, 200)
    monegetti_pred = np.array([monegetti_competency(age, a, b1, v1, b2, tc, Tcp) 
                               for age in age_range])
    
    ax.plot(age_range, monegetti_pred, 'r-', linewidth=3, 
            label='Monegetti piecewise model', zorder=2)
    
    # Mark key points
    ax.axvline(tc, color='orange', linestyle='--', linewidth=2, 
               label=f'tc = {tc:.2f} days', alpha=0.7)
    ax.axvline(Tcp, color='purple', linestyle='--', linewidth=2, 
               label=f'Tcp = {Tcp:.2f} days', alpha=0.7)
    
    # Add TC50 if available
    if not np.isnan(est['TC50']):
        tc50 = est['TC50']
        comp_tc50 = monegetti_competency(tc50, a, b1, v1, b2, tc, Tcp)
        ax.plot(tc50, comp_tc50, 'go', markersize=12, 
                markeredgecolor='black', markeredgewidth=2,
                label=f'TC50 = {tc50:.2f} days', zorder=4)
        ax.axvline(tc50, color='green', linestyle=':', linewidth=1.5, alpha=0.5)
    
    # Formatting
    ax.set_xlabel('Larval Age (days)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Competency (proportion settled)', fontsize=14, fontweight='bold')
    ax.set_title(f'{family} - Monegetti Piecewise Model Fit\n'
                 f'N = {int(est["N_replicates"])} replicates, '
                 f'{int(est["N_larvae"])} larvae, '
                 f'AICc = {est["AICc"]:.1f}',
                 fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, ages.max() + 5)
    ax.set_ylim(-0.05, min(1.1, props.max() * 1.2))
    
    # Add parameter text box
    param_text = (f'Parameters:\n'
                  f'a = {a:.4f}\n'
                  f'b1 = {b1:.4f}\n'
                  f'v1 = {v1:.4f}\n'
                  f'b2 = {b2:.4f}\n'
                  f'tc = {tc:.2f} days\n'
                  f'Tcp = {Tcp:.2f} days')
    
    ax.text(0.98, 0.02, param_text, transform=ax.transAxes,
            fontsize=10, verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    
    # Save figure
    filename = os.path.join(output_dir, f'monegetti_{family}.png')
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {filename}")
    return filename

def create_all_plots():
    """Create plots for all families."""
    print("="*70)
    print("CREATING MONEGETTI MODEL PLOTS")
    print("="*70)
    
    # Load data
    data_path = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity_KERNELS/R_codes/Randal/data/Settlement.csv'
    print("\nLoading settlement data...")
    df = load_all_data(data_path)
    
    # Load estimations
    print("\nLoading estimation results...")
    estimations = load_latest_estimations('estimations')
    print(f"Found results for {len(estimations)} families")
    
    # Create plots
    print("\nCreating plots...")
    for family, est_df in estimations.items():
        print(f"\nProcessing {family}...")
        df_family = df[df['Family'] == family].copy()
        
        if len(df_family) > 0:
            plot_family_monegetti(family, df_family, est_df, 'figures')
        else:
            print(f"  WARNING: No data found for {family}")
    
    print("\n" + "="*70)
    print("All plots created successfully!")
    print("="*70)

if __name__ == "__main__":
    create_all_plots()

