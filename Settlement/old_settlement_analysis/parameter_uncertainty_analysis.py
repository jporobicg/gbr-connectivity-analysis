#!/usr/bin/env python3
"""
Parameter Uncertainty Analysis for Hybrid Competency Model
=========================================================

This script adds uncertainty quantification to our hybrid approach using:
1. Bootstrap resampling
2. Profile likelihood confidence intervals
3. Hessian-based standard errors

This addresses the limitation that our original approach didn't provide
uncertainty estimates for parameters.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize
from scipy.stats import bootstrap
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 10)

# Import our existing functions
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
# UNCERTAINTY QUANTIFICATION METHODS
# ============================================================================

def bootstrap_uncertainty(df_agg, model_type, n_bootstrap=1000, family_name=""):
    """
    Bootstrap uncertainty estimation for model parameters.
    
    Parameters:
    -----------
    df_agg : DataFrame
        Aggregated data for a family
    model_type : str
        'Logistic', 'Gompertz', or 'Weibull'
    n_bootstrap : int
        Number of bootstrap samples
    family_name : str
        Family name for progress reporting
    
    Returns:
    --------
    dict with parameter estimates, standard errors, and confidence intervals
    """
    print(f"  Bootstrapping {model_type} for {family_name}...")
    
    # Original fit
    if model_type == 'Logistic':
        params_orig, ll_orig = fit_logistic_continuous(df_agg)
    elif model_type == 'Gompertz':
        params_orig, ll_orig = fit_gompertz_continuous(df_agg)
    elif model_type == 'Weibull':
        params_orig, ll_orig = fit_weibull_continuous(df_agg)
    
    if params_orig is None:
        return None
    
    # Bootstrap samples
    bootstrap_params = []
    n_obs = len(df_agg)
    
    for i in range(n_bootstrap):
        # Resample with replacement
        bootstrap_indices = np.random.choice(n_obs, size=n_obs, replace=True)
        df_bootstrap = df_agg.iloc[bootstrap_indices].reset_index(drop=True)
        
        # Fit model to bootstrap sample
        try:
            if model_type == 'Logistic':
                params_boot, _ = fit_logistic_continuous(df_bootstrap)
            elif model_type == 'Gompertz':
                params_boot, _ = fit_gompertz_continuous(df_bootstrap)
            elif model_type == 'Weibull':
                params_boot, _ = fit_weibull_continuous(df_bootstrap)
            
            if params_boot is not None:
                # Check for reasonable parameter values
                if model_type == 'Logistic':
                    L, k, x0 = params_boot
                    if 0 < L <= 1 and k > 0 and -50 < x0 < 50:
                        bootstrap_params.append(params_boot)
                elif model_type == 'Gompertz':
                    a, b, c = params_boot
                    if 0 < a <= 1 and 0 < b < 1e6 and 0 < c < 10:
                        bootstrap_params.append(params_boot)
                elif model_type == 'Weibull':
                    tc, a, b, v = params_boot
                    if 0 < tc < 20 and 0 < a < 10 and 0 < b < 1e6 and 0 < v < 10:
                        bootstrap_params.append(params_boot)
        except:
            continue
    
    if len(bootstrap_params) < 50:  # Need sufficient successful fits
        print(f"    Warning: Only {len(bootstrap_params)} successful bootstrap fits")
        return None
    
    bootstrap_params = np.array(bootstrap_params)
    
    # Calculate statistics (only for finite values)
    param_means = []
    param_stds = []
    param_cis = []
    
    for i in range(bootstrap_params.shape[1]):
        param_values = bootstrap_params[:, i]
        finite_mask = np.isfinite(param_values)
        
        if np.sum(finite_mask) > 10:
            finite_values = param_values[finite_mask]
            param_means.append(np.mean(finite_values))
            param_stds.append(np.std(finite_values))
            param_cis.append(np.percentile(finite_values, [2.5, 97.5]))
        else:
            param_means.append(np.nan)
            param_stds.append(np.nan)
            param_cis.append([np.nan, np.nan])
    
    param_means = np.array(param_means)
    param_stds = np.array(param_stds)
    param_cis = np.array(param_cis)
    
    return {
        'original_params': params_orig,
        'bootstrap_means': param_means,
        'bootstrap_stds': param_stds,
        'bootstrap_cis': param_cis,
        'n_successful': len(bootstrap_params),
        'bootstrap_params': bootstrap_params
    }

def profile_likelihood_ci(df_agg, model_type, param_idx, param_range, n_points=20):
    """
    Calculate profile likelihood confidence interval for a single parameter.
    
    Parameters:
    -----------
    df_agg : DataFrame
        Aggregated data
    model_type : str
        Model type
    param_idx : int
        Index of parameter to profile
    param_range : tuple
        (min, max) range for parameter
    n_points : int
        Number of points to evaluate
    
    Returns:
    --------
    dict with parameter values and likelihood profile
    """
    # Get original fit
    if model_type == 'Logistic':
        params_orig, ll_orig = fit_logistic_continuous(df_agg)
    elif model_type == 'Gompertz':
        params_orig, ll_gomp = fit_gompertz_continuous(df_agg)
    elif model_type == 'Weibull':
        params_orig, ll_orig = fit_weibull_continuous(df_agg)
    
    if params_orig is None:
        return None
    
    # Profile likelihood
    param_values = np.linspace(param_range[0], param_range[1], n_points)
    profile_ll = []
    
    for param_val in param_values:
        # Fix parameter and optimize others
        def neg_ll_fixed(param_val_fixed):
            params_test = params_orig.copy()
            params_test[param_idx] = param_val_fixed
            
            if model_type == 'Logistic':
                L, k, x0 = params_test
                if L <= 0 or L > 1 or k <= 0:
                    return 1e10
                pred = logistic_competency(df_agg['LarvalAge'].values, L, k, x0)
            elif model_type == 'Gompertz':
                a, b, c = params_test
                if a <= 0 or a > 1 or b <= 0 or c <= 0:
                    return 1e10
                pred = gompertz_competency(df_agg['LarvalAge'].values, a, b, c)
            elif model_type == 'Weibull':
                tc, a, b, v = params_test
                if tc < 0 or a <= 0 or b <= 0 or v <= 0:
                    return 1e10
                pred = weibull_competency(df_agg['LarvalAge'].values, tc, a, b, v)
            
            pred = np.clip(pred, 1e-10, 1 - 1e-10)
            weights = df_agg['NoSet'] + df_agg['NoNotSet']
            ll = np.sum(weights * (
                df_agg['NoSet'].values * np.log(pred) + 
                df_agg['NoNotSet'].values * np.log(1 - pred)
            ))
            return -ll
        
        try:
            result = minimize(neg_ll_fixed, params_orig[param_idx], method='Nelder-Mead')
            profile_ll.append(result.fun)
        except:
            profile_ll.append(1e10)
    
    return {
        'param_values': param_values,
        'profile_ll': profile_ll,
        'original_ll': ll_orig,
        'original_param': params_orig[param_idx]
    }

def hessian_uncertainty(df_agg, model_type):
    """
    Calculate parameter uncertainties using Hessian matrix.
    
    Parameters:
    -----------
    df_agg : DataFrame
        Aggregated data
    model_type : str
        Model type
    
    Returns:
    --------
    dict with parameter estimates and standard errors
    """
    # Get original fit
    if model_type == 'Logistic':
        params_orig, ll_orig = fit_logistic_continuous(df_agg)
    elif model_type == 'Gompertz':
        params_orig, ll_orig = fit_gompertz_continuous(df_agg)
    elif model_type == 'Weibull':
        params_orig, ll_orig = fit_weibull_continuous(df_agg)
    
    if params_orig is None:
        return None
    
    # Define negative log-likelihood function
    def neg_ll(params):
        if model_type == 'Logistic':
            L, k, x0 = params
            if L <= 0 or L > 1 or k <= 0:
                return 1e10
            pred = logistic_competency(df_agg['LarvalAge'].values, L, k, x0)
        elif model_type == 'Gompertz':
            a, b, c = params
            if a <= 0 or a > 1 or b <= 0 or c <= 0:
                return 1e10
            pred = gompertz_competency(df_agg['LarvalAge'].values, a, b, c)
        elif model_type == 'Weibull':
            tc, a, b, v = params
            if tc < 0 or a <= 0 or b <= 0 or v <= 0:
                return 1e10
            pred = weibull_competency(df_agg['LarvalAge'].values, tc, a, b, v)
        
        pred = np.clip(pred, 1e-10, 1 - 1e-10)
        weights = df_agg['NoSet'] + df_agg['NoNotSet']
        ll = np.sum(weights * (
            df_agg['NoSet'].values * np.log(pred) + 
            df_agg['NoNotSet'].values * np.log(1 - pred)
        ))
        return -ll
    
    # Calculate Hessian using numerical differentiation
    from scipy.optimize import approx_fprime
    from scipy.optimize import approx_fprime
    
    def hessian_approx(f, x, eps=1e-5):
        n = len(x)
        hess = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i == j:
                    # Diagonal elements
                    x_plus = x.copy()
                    x_minus = x.copy()
                    x_plus[i] += eps
                    x_minus[i] -= eps
                    hess[i, j] = (f(x_plus) - 2*f(x) + f(x_minus)) / (eps**2)
                else:
                    # Off-diagonal elements
                    x_pp = x.copy()
                    x_pm = x.copy()
                    x_mp = x.copy()
                    x_mm = x.copy()
                    x_pp[i] += eps; x_pp[j] += eps
                    x_pm[i] += eps; x_pm[j] -= eps
                    x_mp[i] -= eps; x_mp[j] += eps
                    x_mm[i] -= eps; x_mm[j] -= eps
                    hess[i, j] = (f(x_pp) - f(x_pm) - f(x_mp) + f(x_mm)) / (4*eps**2)
        return hess
    
    try:
        hess = hessian_approx(neg_ll, params_orig)
        cov_matrix = np.linalg.inv(hess)
        std_errors = np.sqrt(np.diag(cov_matrix))
        
        return {
            'params': params_orig,
            'std_errors': std_errors,
            'cov_matrix': cov_matrix,
            'hessian': hess
        }
    except:
        return None

# ============================================================================
# ANALYSIS
# ============================================================================

def analyze_uncertainty():
    """Main uncertainty analysis."""
    print("=" * 70)
    print("PARAMETER UNCERTAINTY ANALYSIS")
    print("=" * 70)
    
    # Load data
    data_path = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity/R_codes/Randal/data/Settlement.csv'
    family_data = prepare_data_continuous(data_path, treatment='rubble')
    
    # Results storage
    uncertainty_results = {}
    
    for family, df_agg in family_data.items():
        print(f"\nAnalyzing {family}...")
        print(f"  Data: {len(df_agg)} observations, age range {df_agg['LarvalAge'].min():.0f}-{df_agg['LarvalAge'].max():.0f}")
        
        # Determine best model from original analysis
        params_log, ll_log = fit_logistic_continuous(df_agg)
        params_gomp, ll_gomp = fit_gompertz_continuous(df_agg)
        params_weib, ll_weib = fit_weibull_continuous(df_agg)
        
        aic_log = 2 * 3 + 2 * ll_log if params_log is not None else np.inf
        aic_gomp = 2 * 3 + 2 * ll_gomp if params_gomp is not None else np.inf
        aic_weib = 2 * 4 + 2 * ll_weib if params_weib is not None else np.inf
        
        best_model = min(['Logistic', 'Gompertz', 'Weibull'], 
                        key=lambda x: {'Logistic': aic_log, 'Gompertz': aic_gomp, 'Weibull': aic_weib}[x])
        
        print(f"  Best model: {best_model}")
        
        # Bootstrap uncertainty
        bootstrap_result = bootstrap_uncertainty(df_agg, best_model, n_bootstrap=500, family_name=family)
        
        # Hessian uncertainty
        hessian_result = hessian_uncertainty(df_agg, best_model)
        
        uncertainty_results[family] = {
            'best_model': best_model,
            'bootstrap': bootstrap_result,
            'hessian': hessian_result,
            'data': df_agg
        }
    
    return uncertainty_results

def create_uncertainty_table(uncertainty_results):
    """Create summary table with parameter uncertainties."""
    summary_data = []
    
    for family, results in uncertainty_results.items():
        best_model = results['best_model']
        bootstrap = results['bootstrap']
        hessian = results['hessian']
        
        if bootstrap is not None and bootstrap['n_successful'] >= 50:
            params = bootstrap['original_params']
            stds = bootstrap['bootstrap_stds']
            cis = bootstrap['bootstrap_cis']
            
            # Create parameter string with proper formatting
            param_parts = []
            if best_model == 'Logistic':
                L, k, x0 = params
                if not np.isnan(stds[0]): param_parts.append(f"L={L:.3f}±{stds[0]:.3f}")
                if not np.isnan(stds[1]): param_parts.append(f"k={k:.3f}±{stds[1]:.3f}")
                if not np.isnan(stds[2]): param_parts.append(f"x0={x0:.2f}±{stds[2]:.2f}")
                param_str = ", ".join(param_parts)
                tc50 = x0
                if not np.isnan(cis[0,2]) and not np.isnan(cis[1,2]):
                    tc50_ci = f"{cis[0,2]:.1f}-{cis[1,2]:.1f}"
                else:
                    tc50_ci = "N/A"
            elif best_model == 'Gompertz':
                a, b, c = params
                if not np.isnan(stds[0]): param_parts.append(f"a={a:.3f}±{stds[0]:.3f}")
                if not np.isnan(stds[1]): param_parts.append(f"b={b:.1e}±{stds[1]:.1e}")
                if not np.isnan(stds[2]): param_parts.append(f"c={c:.3f}±{stds[2]:.3f}")
                param_str = ", ".join(param_parts)
                tc50 = np.log(b) / c if c > 0 and b > 0 else np.nan
                tc50_ci = "N/A"
            elif best_model == 'Weibull':
                tc, a, b, v = params
                if not np.isnan(stds[0]): param_parts.append(f"tc={tc:.2f}±{stds[0]:.2f}")
                if not np.isnan(stds[1]): param_parts.append(f"a={a:.3f}±{stds[1]:.3f}")
                if not np.isnan(stds[2]): param_parts.append(f"b={b:.3f}±{stds[2]:.3f}")
                if not np.isnan(stds[3]): param_parts.append(f"v={v:.3f}±{stds[3]:.3f}")
                param_str = ", ".join(param_parts)
                tc50 = tc + 2  # Rough approximation
                if not np.isnan(cis[0,0]) and not np.isnan(cis[1,0]):
                    tc50_ci = f"{cis[0,0]:.1f}-{cis[1,0]:.1f}"
                else:
                    tc50_ci = "N/A"
            
            summary_data.append({
                'Family': family,
                'Best_Model': best_model,
                'Parameters': param_str,
                'TC50': f"{tc50:.1f}" if not np.isnan(tc50) else "N/A",
                'TC50_CI': tc50_ci,
                'Bootstrap_N': bootstrap['n_successful'],
                'Hessian_Available': 'Yes' if hessian is not None else 'No'
            })
        else:
            summary_data.append({
                'Family': family,
                'Best_Model': best_model,
                'Parameters': 'Bootstrap failed',
                'TC50': 'N/A',
                'TC50_CI': 'N/A',
                'Bootstrap_N': bootstrap['n_successful'] if bootstrap is not None else 0,
                'Hessian_Available': 'No'
            })
    
    return pd.DataFrame(summary_data)

def plot_uncertainty_results(uncertainty_results, save_path):
    """Create visualization of uncertainty results."""
    n_families = len(uncertainty_results)
    ncols = min(3, n_families)
    nrows = int(np.ceil(n_families / ncols))
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows))
    if n_families == 1:
        axes = np.array([axes])
    axes = axes.flatten() if n_families > 1 else axes
    
    fig.suptitle('Parameter Uncertainty Analysis\n(Bootstrap Confidence Intervals)', 
                 fontsize=16, fontweight='bold', y=0.995)
    
    for idx, (family, results) in enumerate(uncertainty_results.items()):
        if idx >= len(axes):
            break
        
        ax = axes[idx]
        bootstrap = results['bootstrap']
        
        if bootstrap is not None:
            # Plot bootstrap parameter distributions with finite value filtering
            bootstrap_params = bootstrap['bootstrap_params']
            
            # Filter out infinite and NaN values
            finite_params = []
            param_names = []
            
            if results['best_model'] == 'Logistic':
                param_names = ['L (max)', 'k (rate)', 'x0 (midpoint)']
            elif results['best_model'] == 'Gompertz':
                param_names = ['a (max)', 'b (rate)', 'c (shape)']
            elif results['best_model'] == 'Weibull':
                param_names = ['tc (onset)', 'a (rate)', 'b (scale)', 'v (shape)']
            
            for i in range(bootstrap_params.shape[1]):
                param_values = bootstrap_params[:, i]
                # Filter finite values
                finite_mask = np.isfinite(param_values)
                if np.sum(finite_mask) > 10:  # Need at least 10 finite values
                    finite_params.append(param_values[finite_mask])
                else:
                    finite_params.append([])
            
            # Plot only parameters with sufficient finite values
            colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12']
            for i, (param_values, param_name) in enumerate(zip(finite_params, param_names)):
                if len(param_values) > 10:
                    ax.hist(param_values, bins=min(20, len(param_values)//5), 
                           alpha=0.6, label=param_name, density=True, color=colors[i % len(colors)])
            
            ax.set_xlabel('Parameter Value', fontsize=11, fontweight='bold')
            ax.set_ylabel('Density', fontsize=11, fontweight='bold')
            ax.set_title(f'{family}\n{results["best_model"]}', fontsize=12, fontweight='bold')
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, 'Bootstrap\nFailed', ha='center', va='center', 
                   transform=ax.transAxes, fontsize=12)
            ax.set_title(f'{family}', fontsize=12, fontweight='bold')
    
    # Hide extra subplots
    for idx in range(n_families, len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"\nUncertainty plot saved: {save_path}")
    plt.show()

# ============================================================================
# MAIN
# ============================================================================

def main():
    """Main uncertainty analysis."""
    # Run analysis
    uncertainty_results = analyze_uncertainty()
    
    # Create summary table
    summary_df = create_uncertainty_table(uncertainty_results)
    print(f"\n{'=' * 70}")
    print("UNCERTAINTY SUMMARY TABLE")
    print(f"{'=' * 70}\n")
    print(summary_df.to_string(index=False))
    
    # Save results
    output_dir = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity/R_codes/'
    summary_df.to_csv(output_dir + 'parameter_uncertainty_summary.csv', index=False)
    
    # Create visualization
    plot_uncertainty_results(uncertainty_results, 
                           save_path=output_dir + 'parameter_uncertainty_plots.png')
    
    print(f"\n{'=' * 70}")
    print("UNCERTAINTY ANALYSIS COMPLETE")
    print("=" * 70)
    print("\nOutput files:")
    print("  - parameter_uncertainty_summary.csv")
    print("  - parameter_uncertainty_plots.png")
    print("\nNote: Bootstrap method provides robust uncertainty estimates")
    print("      for our hybrid approach parameters.")

if __name__ == "__main__":
    main()
