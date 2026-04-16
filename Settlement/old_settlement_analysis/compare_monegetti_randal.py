#!/usr/bin/env python3
"""
Compare Monegetti's A. tenuis data with Randal's Acroporidae data.
Plot both datasets and their respective fitted Monegetti model curves.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from monegetti_piecewise_model import monegetti_competency

sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 8)

# ============================================================================
# LOAD AND PROCESS DATA
# ============================================================================

def load_monegetti_data(filepath):
    """Load and process Monegetti's A. tenuis data."""
    df = pd.read_csv(filepath)
    
    # Rename columns for consistency
    df = df.rename(columns={
        'age (d)': 'LarvalAge',
        'meta': 'NoSet',
        'larvae': 'NoAlive'
    })
    
    # Calculate proportion metamorphosed (settled)
    df['PropSettled'] = df['NoSet'] / df['NoAlive']
    
    # Aggregate by age (sum across replicates)
    age_summary = df.groupby('LarvalAge').agg({
        'NoSet': 'sum',
        'NoAlive': 'sum'
    }).reset_index()
    age_summary['PropSettled'] = age_summary['NoSet'] / age_summary['NoAlive']
    age_summary = age_summary.sort_values('LarvalAge')
    
    # Add dataset label
    age_summary['Dataset'] = "Monegetti (A. tenuis)"
    
    return age_summary

def load_randal_acroporidae(data_path):
    """Load and process Randal's Acroporidae data."""
    from improved_competency_model import load_all_data
    
    df = load_all_data(data_path)
    
    # Filter to Acroporidae
    df_acro = df[df['Family'] == 'Acroporidae'].copy()
    
    # Aggregate by age
    age_summary = df_acro.groupby('LarvalAge').agg({
        'NoSet': 'sum',
        'NoAlive': 'sum'
    }).reset_index()
    age_summary['PropSettled'] = age_summary['NoSet'] / age_summary['NoAlive']
    age_summary = age_summary.sort_values('LarvalAge')
    
    # Add dataset label
    age_summary['Dataset'] = "Randal (Acroporidae)"
    
    return age_summary

# ============================================================================
# MONEGETTI ORIGINAL PARAMETERS
# ============================================================================

# Monegetti's original parameters for A. tenuis (from literature)
MONEGETTI_ORIGINAL_PARAMS = {
    'a': 1.292,
    'b1': 0.001878,
    'v1': 0.3645,
    'b2': 0.3969,
    'tc': 3.332,
    'Tcp': 69.91
}

# Our Monegetti fit for Randal's Acroporidae data
RANDAL_ACRO_MONEGETTI_PARAMS = {
    'a': 0.3316,
    'b1': 0.0017,
    'v1': 0.1000,
    'b2': 0.0010,
    'tc': 3.9038,
    'Tcp': 71.4755
}

# ============================================================================
# PLOTTING
# ============================================================================

def create_comparison_plot(df_monegetti, df_randal, output_file='figures/monegetti_randal_comparison.png'):
    """Create comparison plot with both datasets and fitted curves."""
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Plot Monegetti's data
    ages_mon = df_monegetti['LarvalAge'].values
    props_mon = df_monegetti['PropSettled'].values
    n_total_mon = df_monegetti['NoAlive'].values
    n_settled_mon = df_monegetti['NoSet'].values
    
    # Calculate 95% CI for Monegetti data
    se_mon = np.sqrt(props_mon * (1 - props_mon) / n_total_mon)
    ci_lower_mon = np.maximum(0, props_mon - 1.96 * se_mon)
    ci_upper_mon = np.minimum(1, props_mon + 1.96 * se_mon)
    
    ax.errorbar(ages_mon, props_mon, 
                yerr=[props_mon - ci_lower_mon, ci_upper_mon - props_mon],
                fmt='o', color='#e74c3c', markersize=10, capsize=5,
                capthick=2, label='Monegetti data (A. tenuis)', 
                zorder=3, alpha=0.8, markeredgecolor='darkred', markeredgewidth=1.5)
    
    # Plot Randal's data
    ages_rand = df_randal['LarvalAge'].values
    props_rand = df_randal['PropSettled'].values
    n_total_rand = df_randal['NoAlive'].values
    n_settled_rand = df_randal['NoSet'].values
    
    # Calculate 95% CI for Randal data
    se_rand = np.sqrt(props_rand * (1 - props_rand) / n_total_rand)
    ci_lower_rand = np.maximum(0, props_rand - 1.96 * se_rand)
    ci_upper_rand = np.minimum(1, props_rand + 1.96 * se_rand)
    
    ax.errorbar(ages_rand, props_rand,
                yerr=[props_rand - ci_lower_rand, ci_upper_rand - props_rand],
                fmt='s', color='#3498db', markersize=8, capsize=5,
                capthick=2, label='Randal data (Acroporidae)', 
                zorder=3, alpha=0.8, markeredgecolor='darkblue', markeredgewidth=1.5)
    
    # Plot Monegetti's original fitted curve
    age_range = np.linspace(0, max(ages_mon.max(), ages_rand.max()) + 5, 300)
    mon_original_pred = np.array([monegetti_competency(
        age, 
        MONEGETTI_ORIGINAL_PARAMS['a'],
        MONEGETTI_ORIGINAL_PARAMS['b1'],
        MONEGETTI_ORIGINAL_PARAMS['v1'],
        MONEGETTI_ORIGINAL_PARAMS['b2'],
        MONEGETTI_ORIGINAL_PARAMS['tc'],
        MONEGETTI_ORIGINAL_PARAMS['Tcp']
    ) for age in age_range])
    
    ax.plot(age_range, mon_original_pred, 'r-', linewidth=3, 
            label='Monegetti original fit (A. tenuis)', zorder=2, alpha=0.8)
    
    # Plot our Monegetti fit for Randal's data
    randal_fit_pred = np.array([monegetti_competency(
        age,
        RANDAL_ACRO_MONEGETTI_PARAMS['a'],
        RANDAL_ACRO_MONEGETTI_PARAMS['b1'],
        RANDAL_ACRO_MONEGETTI_PARAMS['v1'],
        RANDAL_ACRO_MONEGETTI_PARAMS['b2'],
        RANDAL_ACRO_MONEGETTI_PARAMS['tc'],
        RANDAL_ACRO_MONEGETTI_PARAMS['Tcp']
    ) for age in age_range])
    
    ax.plot(age_range, randal_fit_pred, 'b--', linewidth=3,
            label='Our Monegetti fit (Randal Acroporidae)', zorder=2, alpha=0.8)
    
    # Mark key parameters
    # Monegetti original
    ax.axvline(MONEGETTI_ORIGINAL_PARAMS['tc'], color='red', linestyle=':', 
               linewidth=2, alpha=0.5, label=f"Monegetti tc = {MONEGETTI_ORIGINAL_PARAMS['tc']:.2f} days")
    ax.axvline(MONEGETTI_ORIGINAL_PARAMS['Tcp'], color='red', linestyle=':', 
               linewidth=2, alpha=0.5, label=f"Monegetti Tcp = {MONEGETTI_ORIGINAL_PARAMS['Tcp']:.2f} days")
    
    # Randal fit
    ax.axvline(RANDAL_ACRO_MONEGETTI_PARAMS['tc'], color='blue', linestyle=':', 
               linewidth=2, alpha=0.5, label=f"Randal tc = {RANDAL_ACRO_MONEGETTI_PARAMS['tc']:.2f} days")
    ax.axvline(RANDAL_ACRO_MONEGETTI_PARAMS['Tcp'], color='blue', linestyle=':', 
               linewidth=2, alpha=0.5, label=f"Randal Tcp = {RANDAL_ACRO_MONEGETTI_PARAMS['Tcp']:.2f} days")
    
    # Formatting
    ax.set_xlabel('Larval Age (days)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Competency (proportion settled/metamorphosed)', fontsize=14, fontweight='bold')
    ax.set_title('Comparison: Monegetti (A. tenuis) vs Randal (Acroporidae)\n'
                 'Monegetti Piecewise Model Fits',
                 fontsize=16, fontweight='bold')
    ax.legend(loc='best', fontsize=10, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, max(ages_mon.max(), ages_rand.max()) + 5)
    ax.set_ylim(-0.05, 1.05)
    
    # Add parameter comparison text box
    param_text = (
        'Monegetti Original (A. tenuis):\n'
        f'a={MONEGETTI_ORIGINAL_PARAMS["a"]:.3f}, '
        f'b1={MONEGETTI_ORIGINAL_PARAMS["b1"]:.4f}, '
        f'v1={MONEGETTI_ORIGINAL_PARAMS["v1"]:.3f}\n'
        f'b2={MONEGETTI_ORIGINAL_PARAMS["b2"]:.3f}, '
        f'tc={MONEGETTI_ORIGINAL_PARAMS["tc"]:.2f}, '
        f'Tcp={MONEGETTI_ORIGINAL_PARAMS["Tcp"]:.2f}\n\n'
        'Our Fit (Randal Acroporidae):\n'
        f'a={RANDAL_ACRO_MONEGETTI_PARAMS["a"]:.3f}, '
        f'b1={RANDAL_ACRO_MONEGETTI_PARAMS["b1"]:.4f}, '
        f'v1={RANDAL_ACRO_MONEGETTI_PARAMS["v1"]:.3f}\n'
        f'b2={RANDAL_ACRO_MONEGETTI_PARAMS["b2"]:.3f}, '
        f'tc={RANDAL_ACRO_MONEGETTI_PARAMS["tc"]:.2f}, '
        f'Tcp={RANDAL_ACRO_MONEGETTI_PARAMS["Tcp"]:.2f}'
    )
    
    ax.text(0.98, 0.02, param_text, transform=ax.transAxes,
            fontsize=10, verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9),
            family='monospace')
    
    # Add data summary
    summary_text = (
        f'Monegetti: N={len(df_monegetti)} ages, '
        f'{df_monegetti["NoAlive"].sum():.0f} larvae\n'
        f'Randal: N={len(df_randal)} ages, '
        f'{df_randal["NoAlive"].sum():.0f} larvae'
    )
    
    ax.text(0.02, 0.98, summary_text, transform=ax.transAxes,
            fontsize=11, verticalalignment='top', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    
    # Save figure
    import os
    os.makedirs('figures', exist_ok=True)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Saved: {output_file}")
    
    return output_file

# ============================================================================
# MAIN
# ============================================================================

def main():
    """Main function to create comparison."""
    print("="*70)
    print("COMPARING MONEGETTI (A. tenuis) vs RANDAL (Acroporidae)")
    print("="*70)
    
    # Load Monegetti data
    print("\nLoading Monegetti's A. tenuis data...")
    monegetti_file = 'monegetti/A.tenuisGBR2012metamorphosis.csv'
    df_monegetti = load_monegetti_data(monegetti_file)
    print(f"  Ages: {df_monegetti['LarvalAge'].min():.1f} - {df_monegetti['LarvalAge'].max():.1f} days")
    print(f"  Total larvae: {df_monegetti['NoAlive'].sum():.0f}")
    print(f"  Total metamorphosed: {df_monegetti['NoSet'].sum():.0f}")
    
    # Load Randal data
    print("\nLoading Randal's Acroporidae data...")
    randal_file = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity_KERNELS/R_codes/Randal/data/Settlement.csv'
    df_randal = load_randal_acroporidae(randal_file)
    print(f"  Ages: {df_randal['LarvalAge'].min():.1f} - {df_randal['LarvalAge'].max():.1f} days")
    print(f"  Total larvae: {df_randal['NoAlive'].sum():.0f}")
    print(f"  Total settled: {df_randal['NoSet'].sum():.0f}")
    
    # Create combined dataset for export
    print("\nCreating combined dataset...")
    df_combined = pd.concat([
        df_monegetti[['LarvalAge', 'NoSet', 'NoAlive', 'PropSettled', 'Dataset']],
        df_randal[['LarvalAge', 'NoSet', 'NoAlive', 'PropSettled', 'Dataset']]
    ], ignore_index=True)
    df_combined = df_combined.sort_values(['Dataset', 'LarvalAge'])
    
    # Save combined dataset
    output_csv = 'monegetti_randal_combined_data.csv'
    df_combined.to_csv(output_csv, index=False)
    print(f"✓ Saved combined data: {output_csv}")
    
    # Create comparison plot
    print("\nCreating comparison plot...")
    create_comparison_plot(df_monegetti, df_randal)
    
    # Print parameter comparison
    print("\n" + "="*70)
    print("PARAMETER COMPARISON")
    print("="*70)
    print("\nMonegetti Original (A. tenuis):")
    for key, val in MONEGETTI_ORIGINAL_PARAMS.items():
        print(f"  {key:3s} = {val:10.4f}")
    
    print("\nOur Monegetti Fit (Randal Acroporidae):")
    for key, val in RANDAL_ACRO_MONEGETTI_PARAMS.items():
        print(f"  {key:3s} = {val:10.4f}")
    
    print("\n" + "="*70)
    print("Comparison complete!")
    print("="*70)

if __name__ == "__main__":
    main()

