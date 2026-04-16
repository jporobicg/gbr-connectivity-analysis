#!/usr/bin/env python3
"""
Figure Generation for Larval Competency Modeling
================================================

This module contains all plotting code for generating publication-ready figures.
All figures are saved to the figures_analysis/ directory.

Figures generated:
1. Theoretical model comparison (all models, no data)
2. Family-specific fits (data + model + uncertainty)
3. Model diagnostics
4. All families comparison
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, Any, Optional
import os

# Import model functions
from utility_tools import (
    logistic_competency,
    gompertz_competency,
    weibull_competency,
    monegetti_competency
)

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10


# ============================================================================
# THEORETICAL MODEL COMPARISON
# ============================================================================

def create_theoretical_comparison_figure(output_dir: str) -> str:
    """
    Create figure showing theoretical shapes of all models (no fitted data).
    
    Parameters
    ----------
    output_dir : str
        Directory to save figure
        
    Returns
    -------
    str
        Path to saved figure
    """
    # Create age range
    ages = np.linspace(0, 30, 300)
    
    # Example parameters for illustration
    params = {
        'Logistic': [0.8, 0.5, 5.0],  # L, k, x0
        'Gompertz': [0.8, 1.0, 0.3],  # a, b, c
        'Weibull': [3.0, 0.8, 0.1, 1.0],  # tc, a, b, v
        'Monegetti': [0.3, 0.01, 0.1, 0.001, 3.0, 20.0]  # a, b1, v1, b2, tc, Tcp
    }
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Plot each model
    colors = {
        'Logistic': '#e74c3c',
        'Gompertz': '#3498db',
        'Weibull': '#2ecc71',
        'Monegetti': '#9b59b6'
    }
    
    for model_name, color in colors.items():
        if model_name == 'Logistic':
            pred = logistic_competency(ages, *params[model_name])
        elif model_name == 'Gompertz':
            pred = gompertz_competency(ages, *params[model_name])
        elif model_name == 'Weibull':
            pred = weibull_competency(ages, *params[model_name])
        else:  # Monegetti
            pred = np.array([
                monegetti_competency(age, *params[model_name])
                for age in ages
            ])
        
        ax.plot(ages, pred, linewidth=3, label=model_name, color=color, alpha=0.8)
    
    ax.set_xlabel('Larval Age (days)', fontweight='bold')
    ax.set_ylabel('Competency', fontweight='bold')
    ax.set_title('Theoretical Model Shapes\n(For Illustration Only)', 
                 fontweight='bold', fontsize=14)
    ax.legend(loc='best', frameon=True, fancybox=True, shadow=True)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 30)
    ax.set_ylim(0, 1)
    
    plt.tight_layout()
    
    # Save
    os.makedirs(output_dir, exist_ok=True)
    filename = os.path.join(output_dir, 'theoretical_model_comparison.png')
    plt.savefig(filename, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"  ✓ Saved: {filename}")
    return filename


# ============================================================================
# FAMILY-SPECIFIC FIGURES
# ============================================================================

def create_family_figure(
    family: str,
    df_family: pd.DataFrame,
    result: Dict[str, Any],
    output_dir: str
) -> str:
    """
    Create figure for a single family showing data, fitted model, and diagnostics.
    
    Parameters
    ----------
    family : str
        Family name
    df_family : pd.DataFrame
        Family data
    result : Dict[str, Any]
        Analysis result
    output_dir : str
        Output directory
        
    Returns
    -------
    str
        Path to saved figure
    """
    # Aggregate data by age
    age_summary = df_family.groupby('LarvalAge').agg({
        'NoSet': 'sum',
        'NoAlive': 'sum'
    }).reset_index()
    age_summary['PropSettled'] = age_summary['NoSet'] / age_summary['NoAlive']
    age_summary = age_summary.sort_values('LarvalAge')
    
    # Create figure with subplots
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
    
    # Main plot: Data and fitted model
    ax_main = fig.add_subplot(gs[0, :])
    
    # Plot data with error bars
    ages = age_summary['LarvalAge'].values
    props = age_summary['PropSettled'].values
    n_total = age_summary['NoAlive'].values
    n_settled = age_summary['NoSet'].values
    
    # Calculate 95% CI
    se = np.sqrt(props * (1 - props) / n_total)
    ci_lower = np.maximum(0, props - 1.96 * se)
    ci_upper = np.minimum(1, props + 1.96 * se)
    
    ax_main.errorbar(
        ages, props,
        yerr=[props - ci_lower, ci_upper - props],
        fmt='o', color='#2c3e50', markersize=8, capsize=5,
        capthick=2, label='Observed data', zorder=3, alpha=0.8
    )
    
    # Plot fitted model
    best_model = result['best_model']
    best_result = result['best_model_result']
    params = best_result['params']
    model_func = {
        'Logistic': logistic_competency,
        'Gompertz': gompertz_competency,
        'Weibull': weibull_competency,
        'Monegetti': monegetti_competency
    }[best_model]
    
    age_range = np.linspace(0, ages.max() + 5, 300)
    if best_model == 'Monegetti':
        pred = np.array([model_func(age, *params) for age in age_range])
    else:
        pred = model_func(age_range, *params)
    
    ax_main.plot(
        age_range, pred, 'r-', linewidth=3,
        label=f'{best_model} model (AICc = {best_result["aicc"]:.1f})',
        zorder=2, alpha=0.8
    )
    
    # Mark TC50 if available
    if result['tc50'] is not None:
        tc50 = result['tc50']
        if best_model == 'Monegetti':
            comp_tc50 = monegetti_competency(tc50, *params)
        else:
            comp_tc50 = model_func(np.array([tc50]), *params)[0]
        ax_main.plot(
            tc50, comp_tc50, 'go', markersize=12,
            markeredgecolor='black', markeredgewidth=2,
            label=f'TC50 = {tc50:.2f} days', zorder=4
        )
        ax_main.axvline(tc50, color='green', linestyle=':', linewidth=1.5, alpha=0.5)
    
    ax_main.set_xlabel('Larval Age (days)', fontweight='bold')
    ax_main.set_ylabel('Competency (proportion settled)', fontweight='bold')
    ax_main.set_title(
        f'{family} - {best_model} Model Fit\n'
        f'N = {result["n_replicates"]} replicates, '
        f'{result["n_larvae"]:.0f} larvae',
        fontweight='bold', fontsize=14
    )
    ax_main.legend(loc='best', fontsize=10)
    ax_main.grid(True, alpha=0.3)
    ax_main.set_xlim(0, ages.max() + 5)
    ax_main.set_ylim(-0.05, min(1.1, props.max() * 1.2))
    
    # Residuals plot
    ax_resid = fig.add_subplot(gs[1, 0])
    
    # Recalculate residuals on aggregated data to match plot
    from utility_tools import calculate_residuals
    age_agg = age_summary['LarvalAge'].values
    n_settled_agg = age_summary['NoSet'].values
    n_total_agg = age_summary['NoAlive'].values
    
    # Create temporary dataframe for residual calculation
    df_temp = pd.DataFrame({
        'LarvalAge': age_agg,
        'NoSet': n_settled_agg,
        'NoAlive': n_total_agg
    })
    residuals_agg = calculate_residuals(df_temp, params, model_func)
    
    # Ensure arrays are same size
    if len(age_agg) == len(residuals_agg['pearson']):
        ax_resid.scatter(age_agg, residuals_agg['pearson'], alpha=0.6, color='#3498db')
    else:
        # Fallback: use original residuals if size mismatch
        residuals_orig = best_result['residuals']
        df_orig_ages = df_family['LarvalAge'].values
        if len(df_orig_ages) == len(residuals_orig['pearson']):
            ax_resid.scatter(df_orig_ages, residuals_orig['pearson'], 
                           alpha=0.6, color='#3498db', s=10)
    
    ax_resid.axhline(0, color='red', linestyle='--', linewidth=2)
    ax_resid.axhline(2, color='gray', linestyle=':', linewidth=1)
    ax_resid.axhline(-2, color='gray', linestyle=':', linewidth=1)
    ax_resid.set_xlabel('Larval Age (days)', fontweight='bold')
    ax_resid.set_ylabel('Pearson Residuals', fontweight='bold')
    ax_resid.set_title('Residuals Plot', fontweight='bold')
    ax_resid.grid(True, alpha=0.3)
    
    # Q-Q plot
    ax_qq = fig.add_subplot(gs[1, 1])
    from scipy import stats
    if len(age_agg) == len(residuals_agg['pearson']):
        pearson_resid = residuals_agg['pearson']
    else:
        pearson_resid = best_result['residuals']['pearson']
    
    if len(pearson_resid) > 0:
        stats.probplot(pearson_resid, dist="norm", plot=ax_qq)
    ax_qq.set_title('Q-Q Plot (Normality Check)', fontweight='bold')
    ax_qq.grid(True, alpha=0.3)
    
    # Add parameter text
    param_text = f'{best_model} Parameters:\n'
    param_names = best_result['param_names']
    for name, value in zip(param_names, params):
        param_text += f'{name} = {value:.4f}\n'
    param_text += f'\nAICc = {best_result["aicc"]:.1f}\n'
    param_text += f'φ (overdispersion) = {best_result["phi"]:.2f}'
    
    fig.text(0.98, 0.02, param_text, transform=fig.transFigure,
             fontsize=9, verticalalignment='bottom', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
             family='monospace')
    
    plt.tight_layout()
    
    # Save
    os.makedirs(output_dir, exist_ok=True)
    sanitized = family.replace(' ', '_').replace('/', '_')
    filename = os.path.join(output_dir, f'family_{sanitized}_fit.png')
    plt.savefig(filename, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"  ✓ Saved: {filename}")
    return filename


# ============================================================================
# ALL FAMILIES COMPARISON
# ============================================================================

def create_all_families_figure(
    all_results: Dict[str, Dict[str, Any]],
    all_data: Optional[Dict[str, pd.DataFrame]] = None,
    output_dir: str = 'figures_analysis'
) -> str:
    """
    Create figure comparing all families.
    
    Parameters
    ----------
    all_results : Dict[str, Dict[str, Any]]
        Results for all families
    all_data : Optional[Dict[str, pd.DataFrame]], optional
        Data for each family, by default None
    output_dir : str
        Output directory
        
    Returns
    -------
    str
        Path to saved figure
    """
    n_families = len(all_results)
    n_cols = 3
    n_rows = (n_families + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
    if n_families == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    for idx, (family, result) in enumerate(all_results.items()):
        ax = axes[idx]
        
        best_model = result['best_model']
        best_result = result['best_model_result']
        params = best_result['params']
        
        model_func = {
            'Logistic': logistic_competency,
            'Gompertz': gompertz_competency,
            'Weibull': weibull_competency,
            'Monegetti': monegetti_competency
        }[best_model]
        
        # Plot data if available
        if all_data is not None and family in all_data:
            df_fam = all_data[family]
            age_summary = df_fam.groupby('LarvalAge').agg({
                'NoSet': 'sum',
                'NoAlive': 'sum'
            }).reset_index()
            age_summary['PropSettled'] = age_summary['NoSet'] / age_summary['NoAlive']
            age_summary = age_summary.sort_values('LarvalAge')
            
            ages = age_summary['LarvalAge'].values
            props = age_summary['PropSettled'].values
            n_total = age_summary['NoAlive'].values
            
            # Calculate 95% CI
            se = np.sqrt(props * (1 - props) / n_total)
            ci_lower = np.maximum(0, props - 1.96 * se)
            ci_upper = np.minimum(1, props + 1.96 * se)
            
            ax.errorbar(
                ages, props,
                yerr=[props - ci_lower, ci_upper - props],
                fmt='o', color='#2c3e50', markersize=4, capsize=3,
                alpha=0.6, zorder=2, label='Data'
            )
            max_age = ages.max()
        else:
            max_age = 30
        
        # Plot fitted model
        age_range = np.linspace(0, max_age, 300)
        if best_model == 'Monegetti':
            pred = np.array([model_func(age, *params) for age in age_range])
        else:
            pred = model_func(age_range, *params)
        
        ax.plot(age_range, pred, linewidth=2.5, label=best_model, 
               color='#e74c3c', zorder=3)
        ax.set_title(f'{family}\n{best_model} (AICc={best_result["aicc"]:.1f})', 
                    fontweight='bold', fontsize=11)
        ax.set_xlabel('Age (days)', fontsize=10)
        ax.set_ylabel('Competency', fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, max_age)
        ax.set_ylim(0, 1)
        ax.legend(loc='best', fontsize=8)
    
    # Hide unused subplots
    for idx in range(n_families, len(axes)):
        axes[idx].set_visible(False)
    
    plt.suptitle('All Families - Best Model Fits', fontweight='bold', fontsize=16)
    plt.tight_layout()
    
    # Save
    os.makedirs(output_dir, exist_ok=True)
    filename = os.path.join(output_dir, 'all_families_comparison.png')
    plt.savefig(filename, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"  ✓ Saved: {filename}")
    return filename

