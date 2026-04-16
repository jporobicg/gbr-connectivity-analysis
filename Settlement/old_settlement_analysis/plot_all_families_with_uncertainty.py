#!/usr/bin/env python3
"""
Create combined visualization of all coral families with uncertainty.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')

sns.set_style("whitegrid")

# Import functions from hybrid model
from hybrid_competency_model import (
    prepare_data_continuous, 
    logistic_competency, 
    gompertz_competency, 
    weibull_competency,
    fit_logistic_continuous,
    fit_gompertz_continuous, 
    fit_weibull_continuous
)

# Load data
data_path = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity_KERNELS/R_codes/Randal/data/Settlement.csv'
family_data = prepare_data_continuous(data_path, treatment='rubble')

# Load summary with best models
summary_df = pd.read_csv('parameter_uncertainty_summary.csv')

# Colors for each family
family_colors = {
    'Acroporidae': '#e74c3c',
    'Poritidae': '#3498db',
    'Merulinidae': '#2ecc71',
    'Diploastreidae': '#f39c12',
    'Agariciidae': '#9b59b6',
    'Lobophylliidae': '#1abc9c',
    'Euphylliidae': '#e91e63'
}

# Create figure
fig, axes = plt.subplots(3, 3, figsize=(18, 15))
axes = axes.flatten()

for idx, (family, df_agg) in enumerate(family_data.items()):
    ax = axes[idx]
    
    # Get best model info
    family_info = summary_df[summary_df['Family'] == family].iloc[0]
    model_type = family_info['Best_Model']
    n_samples = family_info['N_samples']
    tc50 = family_info['TC50']
    tc50_ci = family_info['TC50_CI']
    
    # Fit the best model
    if model_type == 'Logistic':
        params, ll = fit_logistic_continuous(df_agg)
    elif model_type == 'Gompertz':
        params, ll = fit_gompertz_continuous(df_agg)
    elif model_type == 'Weibull':
        params, ll = fit_weibull_continuous(df_agg)
    
    if params is None:
        continue
    
    # Plot observed data
    sizes = (df_agg['NoSet'] + df_agg['NoNotSet']) * 3
    ax.scatter(df_agg['LarvalAge'], df_agg['PropCompetent'],
              s=sizes, alpha=0.6, color=family_colors[family], 
              edgecolors='white', linewidth=1.5, zorder=5, label='Observed')
    
    # Plot model prediction
    age_range = np.linspace(0, min(80, df_agg['LarvalAge'].max() + 10), 200)
    if model_type == 'Logistic':
        pred = logistic_competency(age_range, *params)
    elif model_type == 'Gompertz':
        pred = gompertz_competency(age_range, *params)
    elif model_type == 'Weibull':
        pred = weibull_competency(age_range, *params)
    
    ax.plot(age_range, pred, linewidth=3, color=family_colors[family], 
           label=f'{model_type}', zorder=4)
    
    # Add confidence bands using bootstrap
    # Quick bootstrap for visualization (50 iterations)
    bootstrap_curves = []
    n_obs = len(df_agg)
    
    for _ in range(50):
        try:
            indices = np.random.choice(n_obs, size=n_obs, replace=True)
            df_boot = df_agg.iloc[indices].reset_index(drop=True)
            
            if model_type == 'Logistic':
                params_boot, _ = fit_logistic_continuous(df_boot)
                if params_boot is not None:
                    pred_boot = logistic_competency(age_range, *params_boot)
                    bootstrap_curves.append(pred_boot)
            elif model_type == 'Gompertz':
                params_boot, _ = fit_gompertz_continuous(df_boot)
                if params_boot is not None:
                    pred_boot = gompertz_competency(age_range, *params_boot)
                    bootstrap_curves.append(pred_boot)
            elif model_type == 'Weibull':
                params_boot, _ = fit_weibull_continuous(df_boot)
                if params_boot is not None:
                    pred_boot = weibull_competency(age_range, *params_boot)
                    bootstrap_curves.append(pred_boot)
        except:
            continue
    
    # Plot confidence bands
    if len(bootstrap_curves) > 10:
        bootstrap_curves = np.array(bootstrap_curves)
        lower = np.percentile(bootstrap_curves, 2.5, axis=0)
        upper = np.percentile(bootstrap_curves, 97.5, axis=0)
        ax.fill_between(age_range, lower, upper, alpha=0.3, 
                        color=family_colors[family], zorder=3,
                        label='95% CI')
    
    # Formatting
    ax.set_xlabel('Larval Age (days)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Proportion Competent', fontsize=11, fontweight='bold')
    ax.set_title(f'{family}\n{model_type} (n={n_samples}, TC50={tc50})', 
                fontsize=12, fontweight='bold')
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, age_range[-1])
    ax.set_ylim(0, 1.0)

# Remove empty subplots
for idx in range(len(family_data), len(axes)):
    fig.delaxes(axes[idx])

plt.tight_layout()
plt.savefig('all_families_competency_with_uncertainty.png', dpi=300, bbox_inches='tight')
print("✓ Combined plot saved: all_families_competency_with_uncertainty.png")
plt.close()

# Create a second plot: All families on one panel
fig, ax = plt.subplots(figsize=(12, 8))

for family, df_agg in family_data.items():
    # Get best model info
    family_info = summary_df[summary_df['Family'] == family].iloc[0]
    model_type = family_info['Best_Model']
    n_samples = family_info['N_samples']
    tc50 = family_info['TC50']
    
    # Fit the best model
    if model_type == 'Logistic':
        params, ll = fit_logistic_continuous(df_agg)
    elif model_type == 'Gompertz':
        params, ll = fit_gompertz_continuous(df_agg)
    elif model_type == 'Weibull':
        params, ll = fit_weibull_continuous(df_agg)
    
    if params is None:
        continue
    
    # Plot model prediction
    age_range = np.linspace(0, 80, 200)
    if model_type == 'Logistic':
        pred = logistic_competency(age_range, *params)
    elif model_type == 'Gompertz':
        pred = gompertz_competency(age_range, *params)
    elif model_type == 'Weibull':
        pred = weibull_competency(age_range, *params)
    
    ax.plot(age_range, pred, linewidth=3, color=family_colors[family], 
           label=f'{family} (TC50={tc50}, n={n_samples})', zorder=4)

ax.set_xlabel('Larval Age (days)', fontsize=14, fontweight='bold')
ax.set_ylabel('Proportion Competent', fontsize=14, fontweight='bold')
ax.set_title('Coral Family Competency Functions with Uncertainty', 
            fontsize=16, fontweight='bold')
ax.legend(loc='lower right', fontsize=11, framealpha=0.9)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 80)
ax.set_ylim(0, 1.0)

plt.tight_layout()
plt.savefig('all_families_overlay_competency.png', dpi=300, bbox_inches='tight')
print("✓ Overlay plot saved: all_families_overlay_competency.png")
plt.close()

print("\n✓ All plots created successfully!")
print(f"  - Individual family plots: family_results/<Family>_uncertainty.png")
print(f"  - Combined grid plot: all_families_competency_with_uncertainty.png")
print(f"  - Overlay comparison: all_families_overlay_competency.png")

