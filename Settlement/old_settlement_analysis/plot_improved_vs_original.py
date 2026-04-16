#!/usr/bin/env python3
"""
Compare improved vs original competency models.
Create diagnostic plots and visualizations.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle

sns.set_style("whitegrid")

# Load improved results
with open('improved_competency_results.pkl', 'rb') as f:
    improved_results = pickle.load(f)

# Load original summary
original_df = pd.read_csv('parameter_uncertainty_summary.csv')

# Load data
from improved_competency_model import (
    load_all_data, logistic_competency, gompertz_competency, 
    weibull_competency, estimate_overdispersion
)

data_path = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity_KERNELS/R_codes/Randal/data/Settlement.csv'
df_all = load_all_data(data_path)

# Colors
family_colors = {
    'Acroporidae': '#e74c3c',
    'Poritidae': '#3498db',
    'Merulinidae': '#2ecc71',
    'Diploastreidae': '#f39c12',
    'Agariciidae': '#9b59b6',
    'Lobophylliidae': '#1abc9c',
    'Euphylliidae': '#e91e63'
}

# ============================================================================
# 1. COMPARISON PLOT: Improved vs Original
# ============================================================================

fig, axes = plt.subplots(3, 3, figsize=(18, 15))
axes = axes.flatten()

for idx, (family, result) in enumerate(improved_results.items()):
    ax = axes[idx]
    
    # Get data for this family
    df_fam = df_all[df_all['Family'] == family].copy()
    
    # Aggregate for plotting
    df_agg = df_fam.groupby('LarvalAge').agg({
        'NoSet': 'sum',
        'NoAlive': 'sum'
    }).reset_index()
    df_agg['PropSettled'] = df_agg['NoSet'] / df_agg['NoAlive']
    df_agg['SE'] = np.sqrt(df_agg['PropSettled'] * (1 - df_agg['PropSettled']) / df_agg['NoAlive'])
    
    # Plot observed data with error bars
    ax.errorbar(df_agg['LarvalAge'], df_agg['PropSettled'], 
               yerr=1.96 * df_agg['SE'],
               fmt='o', color=family_colors[family], alpha=0.6,
               capsize=3, markersize=6, label='Observed ± 95% CI',
               zorder=5)
    
    # Plot improved model
    age_range = np.linspace(0, min(60, df_agg['LarvalAge'].max() + 10), 200)
    model_type = result['model']
    params = result['params']
    
    if model_type == 'Logistic':
        pred = logistic_competency(age_range, *params)
    elif model_type == 'Gompertz':
        pred = gompertz_competency(age_range, *params)
    else:  # Weibull
        pred = weibull_competency(age_range, *params)
    
    ax.plot(age_range, pred, linewidth=3, color=family_colors[family],
           label=f'Improved ({model_type})', zorder=4)
    
    # Plot bootstrap uncertainty
    bootstrap_params = result['bootstrap']['bootstrap_params']
    # Sample 50 curves
    for i in np.random.choice(len(bootstrap_params), min(50, len(bootstrap_params)), replace=False):
        p = bootstrap_params[i]
        if model_type == 'Logistic':
            pred_boot = logistic_competency(age_range, *p)
        elif model_type == 'Gompertz':
            pred_boot = gompertz_competency(age_range, *p)
        else:
            pred_boot = weibull_competency(age_range, *p)
        ax.plot(age_range, pred_boot, linewidth=0.5, color=family_colors[family],
               alpha=0.1, zorder=2)
    
    # Formatting
    tc50 = result['tc50']
    phi = result['phi']
    n_rep = result['n_replicates']
    
    ax.set_xlabel('Larval Age (days)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Proportion Competent', fontsize=11, fontweight='bold')
    ax.set_title(f'{family}\n{model_type}, n={n_rep}, TC50={tc50:.1f}d, φ={phi:.2f}',
                fontsize=12, fontweight='bold')
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, age_range[-1])
    ax.set_ylim(0, min(1.0, df_agg['PropSettled'].max() * 1.2))

# Remove empty subplots
for idx in range(len(improved_results), len(axes)):
    fig.delaxes(axes[idx])

plt.tight_layout()
plt.savefig('improved_competency_fits.png', dpi=300, bbox_inches='tight')
print("✓ Saved: improved_competency_fits.png")
plt.close()

# ============================================================================
# 2. DIAGNOSTIC PLOTS: Residuals
# ============================================================================

fig, axes = plt.subplots(3, 3, figsize=(18, 15))
axes = axes.flatten()

for idx, (family, result) in enumerate(improved_results.items()):
    ax = axes[idx]
    
    # Get data
    df_fam = df_all[df_all['Family'] == family].copy()
    
    # Calculate residuals
    model_type = result['model']
    params = result['params']
    
    if model_type == 'Logistic':
        model_func = logistic_competency
    elif model_type == 'Gompertz':
        model_func = gompertz_competency
    else:
        model_func = weibull_competency
    
    phi, pearson_resid = estimate_overdispersion(df_fam, params, model_func)
    
    # Get predicted values for x-axis
    pred_prop = model_func(df_fam['LarvalAge'].values, *params)
    
    # Plot residuals vs fitted
    ax.scatter(pred_prop, pearson_resid, alpha=0.3, s=20,
              color=family_colors[family])
    ax.axhline(0, color='black', linestyle='--', linewidth=1)
    ax.axhline(2, color='red', linestyle=':', linewidth=1, alpha=0.5)
    ax.axhline(-2, color='red', linestyle=':', linewidth=1, alpha=0.5)
    
    ax.set_xlabel('Fitted Proportion', fontsize=10)
    ax.set_ylabel('Pearson Residual', fontsize=10)
    ax.set_title(f'{family} - Residual Plot\n(φ = {phi:.2f})', 
                fontsize=11, fontweight='bold')
    ax.grid(True, alpha=0.3)

# Remove empty subplots
for idx in range(len(improved_results), len(axes)):
    fig.delaxes(axes[idx])

plt.tight_layout()
plt.savefig('improved_competency_diagnostics.png', dpi=300, bbox_inches='tight')
print("✓ Saved: improved_competency_diagnostics.png")
plt.close()

# ============================================================================
# 3. OVERLAY COMPARISON
# ============================================================================

fig, ax = plt.subplots(figsize=(12, 8))

age_range = np.linspace(0, 40, 200)

for family, result in improved_results.items():
    model_type = result['model']
    params = result['params']
    tc50 = result['tc50']
    
    if model_type == 'Logistic':
        pred = logistic_competency(age_range, *params)
    elif model_type == 'Gompertz':
        pred = gompertz_competency(age_range, *params)
    else:
        pred = weibull_competency(age_range, *params)
    
    # Only plot reasonable curves (TC50 < 30 days)
    if tc50 < 30:
        ax.plot(age_range, pred, linewidth=3, color=family_colors[family],
               label=f'{family} (TC50={tc50:.1f}d)', zorder=4)

ax.set_xlabel('Larval Age (days)', fontsize=14, fontweight='bold')
ax.set_ylabel('Proportion Competent', fontsize=14, fontweight='bold')
ax.set_title('Improved Competency Functions (All Treatments Pooled)',
            fontsize=16, fontweight='bold')
ax.legend(loc='lower right', fontsize=11, framealpha=0.9)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 40)
ax.set_ylim(0, 1.0)

plt.tight_layout()
plt.savefig('improved_competency_overlay.png', dpi=300, bbox_inches='tight')
print("✓ Saved: improved_competency_overlay.png")
plt.close()

# ============================================================================
# 4. COMPARISON TABLE
# ============================================================================

print("\n" + "="*80)
print("COMPARISON: Improved vs Original")
print("="*80)

comparison_data = []
for family in improved_results.keys():
    # Improved results
    imp = improved_results[family]
    
    # Original results
    orig = original_df[original_df['Family'] == family]
    if len(orig) > 0:
        orig = orig.iloc[0]
        orig_n = orig['N_samples']
        orig_tc50 = orig['TC50']
    else:
        orig_n = 0
        orig_tc50 = np.nan
    
    comparison_data.append({
        'Family': family,
        'Original_N': orig_n,
        'Improved_N': imp['n_replicates'],
        'Data_Increase': f"{imp['n_replicates']/orig_n:.1f}x" if orig_n > 0 else "N/A",
        'Original_TC50': orig_tc50,
        'Improved_TC50': f"{imp['tc50']:.2f}",
        'Phi': f"{imp['phi']:.2f}",
        'Model': imp['model']
    })

comp_df = pd.DataFrame(comparison_data)
print(comp_df.to_string(index=False))

comp_df.to_csv('comparison_original_vs_improved.csv', index=False)
print("\n✓ Saved: comparison_original_vs_improved.csv")

print("\n" + "="*80)
print("RECOMMENDATIONS")
print("="*80)
print("""
Based on the improved analysis:

1. FAMILIES WITH RELIABLE ESTIMATES (use these):
   - Acroporidae: TC50 = 5.7 days, φ = 5.5, N = 4730 ✓
   - Poritidae: TC50 = 4.3 days, φ = 6.2, N = 1187 ✓
   - Merulinidae: TC50 = 2.0 days, φ = 7.0, N = 3880 ✓
   - Diploastreidae: TC50 = 5.0 days, φ = 7.0, N = 539 ✓
   - Euphylliidae: TC50 = 5.0 days, φ = 5.4, N = 496 ✓

2. FAMILIES WITH ISSUES (use with caution):
   - Agariciidae: Low sample (N=167), very short age range (3-7 days)
   - Lobophylliidae: Unrealistic TC50 (325 days), model may not fit well

3. KEY FINDINGS:
   - All families show overdispersion (φ = 4-7)
   - Most competent by 5-6 days
   - Merulinidae fastest (TC50 = 2 days)
   
4. FOR CONNECTIVITY MODELING:
   - Use the improved equations (binomial likelihood, all data)
   - Consider uncertainty in propagation (bootstrap samples available)
   - May want to use conservative estimates for poorly-sampled families
""")

print("\n✓ All plots and comparisons complete!")

