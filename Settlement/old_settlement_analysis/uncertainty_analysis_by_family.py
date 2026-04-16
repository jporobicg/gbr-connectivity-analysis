#!/usr/bin/env python3
"""
Family-by-Family Uncertainty Analysis
=====================================

This script performs uncertainty analysis for each family separately,
saving individual results and handling errors gracefully.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize
import warnings
import os
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

# ============================================================================
# BOOTSTRAP UNCERTAINTY
# ============================================================================

def bootstrap_family(df_agg, model_type, n_bootstrap=500):
    """Bootstrap uncertainty for a single family."""
    print(f"    Bootstrapping {model_type}...")
    
    # Original fit
    if model_type == 'Logistic':
        params_orig, ll_orig = fit_logistic_continuous(df_agg)
    elif model_type == 'Gompertz':
        params_orig, ll_orig = fit_gompertz_continuous(df_agg)
    elif model_type == 'Weibull':
        params_orig, ll_orig = fit_weibull_continuous(df_agg)
    
    if params_orig is None:
        return None
    
    # Bootstrap
    bootstrap_params = []
    n_obs = len(df_agg)
    
    for i in range(n_bootstrap):
        try:
            # Resample
            indices = np.random.choice(n_obs, size=n_obs, replace=True)
            df_boot = df_agg.iloc[indices].reset_index(drop=True)
            
            # Fit
            if model_type == 'Logistic':
                params_boot, _ = fit_logistic_continuous(df_boot)
            elif model_type == 'Gompertz':
                params_boot, _ = fit_gompertz_continuous(df_boot)
            elif model_type == 'Weibull':
                params_boot, _ = fit_weibull_continuous(df_boot)
            
            # Validate parameters
            if params_boot is not None:
                valid = True
                if model_type == 'Logistic':
                    L, k, x0 = params_boot
                    valid = 0 < L <= 1 and k > 0 and -100 < x0 < 100
                elif model_type == 'Gompertz':
                    a, b, c = params_boot
                    valid = 0 < a <= 1 and 0 < b < 1e8 and 0 < c < 20
                elif model_type == 'Weibull':
                    tc, a, b, v = params_boot
                    valid = 0 < tc < 30 and 0 < a < 20 and 0 <= b < 1e8 and 0 < v < 20
                
                if valid:
                    bootstrap_params.append(params_boot)
        except:
            continue
    
    print(f"      Successful fits: {len(bootstrap_params)}/{n_bootstrap}")
    
    if len(bootstrap_params) < 50:
        return None
    
    bootstrap_params = np.array(bootstrap_params)
    
    # Calculate statistics with finite value filtering
    n_params = bootstrap_params.shape[1]
    param_means = np.zeros(n_params)
    param_stds = np.zeros(n_params)
    param_cis_lower = np.zeros(n_params)
    param_cis_upper = np.zeros(n_params)
    
    for i in range(n_params):
        values = bootstrap_params[:, i]
        finite_mask = np.isfinite(values)
        
        if np.sum(finite_mask) > 10:
            finite_vals = values[finite_mask]
            param_means[i] = np.mean(finite_vals)
            param_stds[i] = np.std(finite_vals)
            param_cis_lower[i] = np.percentile(finite_vals, 2.5)
            param_cis_upper[i] = np.percentile(finite_vals, 97.5)
        else:
            param_means[i] = np.nan
            param_stds[i] = np.nan
            param_cis_lower[i] = np.nan
            param_cis_upper[i] = np.nan
    
    return {
        'original_params': params_orig,
        'bootstrap_means': param_means,
        'bootstrap_stds': param_stds,
        'bootstrap_cis_lower': param_cis_lower,
        'bootstrap_cis_upper': param_cis_upper,
        'n_successful': len(bootstrap_params),
        'bootstrap_params': bootstrap_params
    }

# ============================================================================
# ANALYSIS FUNCTIONS
# ============================================================================

def analyze_single_family(family, df_agg, output_dir):
    """Analyze a single family and save results."""
    print(f"\n{'='*70}")
    print(f"ANALYZING: {family}")
    print(f"{'='*70}")
    print(f"  Observations: {len(df_agg)}")
    print(f"  Age range: {df_agg['LarvalAge'].min():.0f}-{df_agg['LarvalAge'].max():.0f} days")
    print(f"  Max proportion: {df_agg['PropCompetent'].max():.3f}")
    
    # Determine best model
    params_log, ll_log = fit_logistic_continuous(df_agg)
    params_gomp, ll_gomp = fit_gompertz_continuous(df_agg)
    params_weib, ll_weib = fit_weibull_continuous(df_agg)
    
    aic_log = 2 * 3 + 2 * ll_log if params_log is not None else np.inf
    aic_gomp = 2 * 3 + 2 * ll_gomp if params_gomp is not None else np.inf
    aic_weib = 2 * 4 + 2 * ll_weib if params_weib is not None else np.inf
    
    aics = {'Logistic': aic_log, 'Gompertz': aic_gomp, 'Weibull': aic_weib}
    best_model = min(aics, key=aics.get)
    
    print(f"\n  Model AICs:")
    print(f"    Logistic: {aic_log:.1f}")
    print(f"    Gompertz: {aic_gomp:.1f}")
    print(f"    Weibull: {aic_weib:.1f}")
    print(f"  → Best model: {best_model}")
    
    # Bootstrap uncertainty
    bootstrap_result = bootstrap_family(df_agg, best_model, n_bootstrap=500)
    
    if bootstrap_result is None:
        print(f"  WARNING: Bootstrap failed for {family}")
        return None
    
    # Create summary for this family
    params = bootstrap_result['original_params']
    stds = bootstrap_result['bootstrap_stds']
    cis_lower = bootstrap_result['bootstrap_cis_lower']
    cis_upper = bootstrap_result['bootstrap_cis_upper']
    
    # Format parameters
    param_data = []
    if best_model == 'Logistic':
        param_names = ['L (max competency)', 'k (rate)', 'x0 (midpoint age)']
        tc50 = params[2]  # x0 is the TC50
        tc50_lower = cis_lower[2]
        tc50_upper = cis_upper[2]
    elif best_model == 'Gompertz':
        param_names = ['a (asymptote)', 'b (displacement)', 'c (rate)']
        tc50 = np.log(params[1]) / params[2] if params[2] > 0 and params[1] > 0 else np.nan
        tc50_lower = np.nan
        tc50_upper = np.nan
    elif best_model == 'Weibull':
        param_names = ['tc (onset)', 'a (rate)', 'b (scale)', 'v (shape)']
        tc50 = params[0] + 2  # Approximate
        tc50_lower = cis_lower[0]
        tc50_upper = cis_upper[0] + 2
    
    for i, param_name in enumerate(param_names):
        param_data.append({
            'Parameter': param_name,
            'Estimate': params[i],
            'Std_Error': stds[i],
            'CI_Lower': cis_lower[i],
            'CI_Upper': cis_upper[i]
        })
    
    param_df = pd.DataFrame(param_data)
    
    # Save parameter table
    param_file = os.path.join(output_dir, f"{family}_parameters.csv")
    param_df.to_csv(param_file, index=False)
    print(f"\n  Parameter table saved: {param_file}")
    
    # Print summary
    print(f"\n  Parameter Estimates (±SE, 95% CI):")
    for _, row in param_df.iterrows():
        if not np.isnan(row['Std_Error']):
            print(f"    {row['Parameter']:25s}: {row['Estimate']:.4f} ± {row['Std_Error']:.4f}  [{row['CI_Lower']:.4f}, {row['CI_Upper']:.4f}]")
        else:
            print(f"    {row['Parameter']:25s}: {row['Estimate']:.4f} (CI not available)")
    
    # TC50 summary
    if not np.isnan(tc50):
        if not np.isnan(tc50_lower):
            print(f"\n  TC50 (age at 50% competency): {tc50:.2f} days [{tc50_lower:.2f}, {tc50_upper:.2f}]")
        else:
            print(f"\n  TC50 (age at 50% competency): {tc50:.2f} days")
    
    # Create plot
    create_family_plot(family, df_agg, best_model, params, bootstrap_result, output_dir)
    
    # Save summary dict
    summary = {
        'family': family,
        'best_model': best_model,
        'n_obs': len(df_agg),
        'age_range': f"{df_agg['LarvalAge'].min():.0f}-{df_agg['LarvalAge'].max():.0f}",
        'max_prop': df_agg['PropCompetent'].max(),
        'params': params,
        'param_stds': stds,
        'param_cis_lower': cis_lower,
        'param_cis_upper': cis_upper,
        'tc50': tc50,
        'tc50_ci': f"[{tc50_lower:.2f}, {tc50_upper:.2f}]" if not np.isnan(tc50_lower) else "N/A",
        'n_bootstrap': bootstrap_result['n_successful'],
        'aic': aics[best_model]
    }
    
    return summary

def create_family_plot(family, df_agg, model_type, params, bootstrap_result, output_dir):
    """Create plot for a single family."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Panel 1: Data and model fit with bootstrap CIs
    ax = axes[0]
    
    # Observed data
    sizes = (df_agg['NoSet'] + df_agg['NoNotSet']) * 2
    ax.scatter(df_agg['LarvalAge'], df_agg['PropCompetent'],
              s=sizes, alpha=0.6, color='black', edgecolors='white',
              linewidth=1.5, zorder=5, label='Observed')
    
    # Model prediction
    age_range = np.linspace(0, min(80, df_agg['LarvalAge'].max() + 10), 200)
    if model_type == 'Logistic':
        pred = logistic_competency(age_range, *params)
    elif model_type == 'Gompertz':
        pred = gompertz_competency(age_range, *params)
    elif model_type == 'Weibull':
        pred = weibull_competency(age_range, *params)
    
    ax.plot(age_range, pred, linewidth=3, color='#e74c3c', 
           label=f'{model_type} fit', zorder=4)
    
    # Bootstrap confidence bands (if available)
    bootstrap_params = bootstrap_result['bootstrap_params']
    if len(bootstrap_params) > 100:
        # Sample 100 bootstrap curves
        sample_indices = np.random.choice(len(bootstrap_params), 100, replace=False)
        for idx in sample_indices:
            params_boot = bootstrap_params[idx]
            if model_type == 'Logistic':
                pred_boot = logistic_competency(age_range, *params_boot)
            elif model_type == 'Gompertz':
                pred_boot = gompertz_competency(age_range, *params_boot)
            elif model_type == 'Weibull':
                pred_boot = weibull_competency(age_range, *params_boot)
            ax.plot(age_range, pred_boot, linewidth=0.5, color='#e74c3c', 
                   alpha=0.05, zorder=2)
    
    ax.set_xlabel('Larval Age (days)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Proportion Competent', fontsize=12, fontweight='bold')
    ax.set_title(f'{family} - Model Fit with Bootstrap Uncertainty', 
                fontsize=13, fontweight='bold')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, age_range[-1])
    ax.set_ylim(0, 1.0)
    
    # Panel 2: Parameter distributions
    ax = axes[1]
    
    # Filter to finite, plottable parameters
    if model_type == 'Logistic':
        param_names = ['L', 'k', 'x0']
        colors = ['#e74c3c', '#3498db', '#2ecc71']
    elif model_type == 'Gompertz':
        param_names = ['a', 'c']  # Skip b (too large)
        colors = ['#e74c3c', '#2ecc71']
        bootstrap_params = bootstrap_params[:, [0, 2]]
    elif model_type == 'Weibull':
        param_names = ['tc', 'a', 'v']  # Skip b (may be too large)
        colors = ['#e74c3c', '#3498db', '#f39c12']
        bootstrap_params = bootstrap_params[:, [0, 1, 3]]
    
    positions = np.arange(len(param_names))
    
    for i, (param_name, color) in enumerate(zip(param_names, colors)):
        values = bootstrap_params[:, i]
        finite_vals = values[np.isfinite(values)]
        
        if len(finite_vals) > 10:
            # Violin plot
            parts = ax.violinplot([finite_vals], positions=[i], 
                                 widths=0.7, showmeans=True, showextrema=True)
            for pc in parts['bodies']:
                pc.set_facecolor(color)
                pc.set_alpha(0.6)
    
    ax.set_xticks(positions)
    ax.set_xticklabels(param_names, fontsize=11, fontweight='bold')
    ax.set_ylabel('Parameter Value', fontsize=12, fontweight='bold')
    ax.set_title(f'Bootstrap Parameter Distributions\n(n={bootstrap_result["n_successful"]} samples)', 
                fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plot_file = os.path.join(output_dir, f"{family}_uncertainty.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Plot saved: {plot_file}")

# ============================================================================
# MAIN
# ============================================================================

def main():
    """Main family-by-family analysis."""
    print("=" * 70)
    print("FAMILY-BY-FAMILY UNCERTAINTY ANALYSIS")
    print("=" * 70)
    
    # Setup
    data_path = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity/R_codes/Randal/data/Settlement.csv'
    output_dir = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity/R_codes/family_results/'
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    family_data = prepare_data_continuous(data_path, treatment='rubble')
    
    # Analyze each family
    all_summaries = []
    
    for family, df_agg in family_data.items():
        try:
            summary = analyze_single_family(family, df_agg, output_dir)
            if summary is not None:
                all_summaries.append(summary)
        except Exception as e:
            print(f"\n  ERROR analyzing {family}: {str(e)}")
            import traceback
            traceback.print_exc()
            continue
    
    # Create combined summary table
    print(f"\n{'='*70}")
    print("COMBINED SUMMARY")
    print(f"{'='*70}\n")
    
    summary_data = []
    for s in all_summaries:
        summary_data.append({
            'Family': s['family'],
            'Best_Model': s['best_model'],
            'N_obs': s['n_obs'],
            'Age_Range': s['age_range'],
            'TC50': f"{s['tc50']:.2f}" if not np.isnan(s['tc50']) else "N/A",
            'TC50_CI': s['tc50_ci'],
            'Bootstrap_Samples': s['n_bootstrap'],
            'AIC': f"{s['aic']:.1f}"
        })
    
    summary_df = pd.DataFrame(summary_data)
    print(summary_df.to_string(index=False))
    
    # Save combined summary
    summary_file = os.path.join(output_dir, "all_families_summary.csv")
    summary_df.to_csv(summary_file, index=False)
    print(f"\nCombined summary saved: {summary_file}")
    
    print(f"\n{'='*70}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*70}")
    print(f"\nResults saved in: {output_dir}")
    print(f"  - Individual parameter tables: <Family>_parameters.csv")
    print(f"  - Individual plots: <Family>_uncertainty.png")
    print(f"  - Combined summary: all_families_summary.csv")

if __name__ == "__main__":
    main()

