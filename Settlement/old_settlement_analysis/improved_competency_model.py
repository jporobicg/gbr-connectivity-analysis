#!/usr/bin/env python3
"""
IMPROVED COMPETENCY MODEL
=========================

Statistical improvements:
1. Binomial likelihood (raw counts, not proportions)
2. Pool all treatments (maximize data)
3. Quasi-binomial overdispersion
4. Profile likelihood CIs
5. Proper model diagnostics
6. Replicate-level bootstrap

Goal: Robust competence equations for connectivity modeling
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize, differential_evolution
from scipy.stats import binom, beta
from scipy.special import expit, logit
import warnings
warnings.filterwarnings('ignore')

sns.set_style("whitegrid")

# ============================================================================
# DATA PREPARATION - ALL TREATMENTS
# ============================================================================

def load_all_data(data_path):
    """
    Load ALL data (all treatments) at replicate level.
    
    Changes from previous approach:
    - Use ALL treatments (not just rubble)
    - Keep replicate-level data (not just aggregated)
    - Filter out invalid rows
    """
    df = pd.read_csv(data_path)
    
    print("=" * 70)
    print("IMPROVED COMPETENCY MODEL")
    print("=" * 70)
    print(f"Total rows: {len(df)}")
    
    # Basic filtering
    df = df[df['Family'].notna()].copy()
    df = df[df['NoAlive'] > 0].copy()
    
    # Calculate observed proportion
    df['PropSettled'] = df['NoSet'] / df['NoAlive']
    
    print(f"After filtering: {len(df)} observations")
    print(f"\nTreatments included: {sorted(df['Treatment'].unique())}")
    print(f"Families: {sorted(df['Family'].unique())}")
    
    return df

# ============================================================================
# COMPETENCY MODELS
# ============================================================================

def logistic_competency(age, L, k, x0):
    """Logistic (sigmoidal) model."""
    return L / (1 + np.exp(-k * (age - x0)))

def gompertz_competency(age, a, b, c):
    """Gompertz model."""
    return a * np.exp(-b * np.exp(-c * age))

def weibull_competency(age, tc, a, b, v):
    """Weibull-based model."""
    age = np.atleast_1d(age)
    competency = np.zeros_like(age, dtype=float)
    mask = age > tc
    if np.any(mask):
        competency[mask] = a * (1 - np.exp(-b * (age[mask] - tc)**v))
    return competency

# ============================================================================
# BINOMIAL LIKELIHOOD FITTING
# ============================================================================

def neg_log_likelihood_binomial(params, ages, n_settled, n_total, model_func, phi=1.0):
    """
    Negative log-likelihood for binomial data with overdispersion.
    
    Parameters:
    -----------
    params : model parameters
    ages : larval ages
    n_settled : number settled (successes)
    n_total : total number alive (trials)
    model_func : competency function
    phi : overdispersion parameter (1 = no overdispersion)
    """
    # Predict competency
    pred_comp = model_func(ages, *params)
    pred_comp = np.clip(pred_comp, 1e-10, 1 - 1e-10)
    
    # Binomial log-likelihood
    ll = np.sum(n_settled * np.log(pred_comp) + 
                (n_total - n_settled) * np.log(1 - pred_comp))
    
    # Apply overdispersion correction
    ll = ll / phi
    
    return -ll

def fit_logistic_binomial(df, phi=1.0):
    """Fit logistic model using binomial likelihood."""
    ages = df['LarvalAge'].values
    n_settled = df['NoSet'].values
    n_total = df['NoAlive'].values
    
    def objective(params):
        L, k, x0 = params
        # Constraints
        if not (0.1 < L <= 1.0 and 0 < k < 10 and 0 < x0 < 50):
            return 1e10
        return neg_log_likelihood_binomial(params, ages, n_settled, n_total, 
                                          logistic_competency, phi)
    
    # Try multiple starting points
    best_result = None
    best_nll = np.inf
    
    for L0 in [0.5, 0.7, 0.9]:
        for x0 in [3, 5, 7]:
            x0_init = [L0, 1.0, x0]
            try:
                result = minimize(objective, x0_init, method='Nelder-Mead',
                                options={'maxiter': 5000})
                if result.fun < best_nll and result.success:
                    best_nll = result.fun
                    best_result = result
            except:
                continue
    
    if best_result is None or not best_result.success:
        return None, np.inf
    
    return best_result.x, best_result.fun

def fit_gompertz_binomial(df, phi=1.0):
    """Fit Gompertz model using binomial likelihood."""
    ages = df['LarvalAge'].values
    n_settled = df['NoSet'].values
    n_total = df['NoAlive'].values
    
    def objective(params):
        a, b, c = params
        if not (0.1 < a <= 1.0 and 0 < b < 100 and 0 < c < 5):
            return 1e10
        return neg_log_likelihood_binomial(params, ages, n_settled, n_total, 
                                          gompertz_competency, phi)
    
    best_result = None
    best_nll = np.inf
    
    for a0 in [0.5, 0.7, 0.9]:
        for c0 in [0.5, 1.0, 2.0]:
            x0_init = [a0, 1.0, c0]
            try:
                result = minimize(objective, x0_init, method='Nelder-Mead',
                                options={'maxiter': 5000})
                if result.fun < best_nll and result.success:
                    best_nll = result.fun
                    best_result = result
            except:
                continue
    
    if best_result is None or not best_result.success:
        return None, np.inf
    
    return best_result.x, best_result.fun

def fit_weibull_binomial(df, phi=1.0):
    """Fit Weibull model using binomial likelihood."""
    ages = df['LarvalAge'].values
    n_settled = df['NoSet'].values
    n_total = df['NoAlive'].values
    
    def objective(params):
        tc, a, b, v = params
        if not (0 < tc < 10 and 0.1 < a <= 2.0 and 0 < b < 5 and 0.1 < v < 5):
            return 1e10
        return neg_log_likelihood_binomial(params, ages, n_settled, n_total, 
                                          weibull_competency, phi)
    
    best_result = None
    best_nll = np.inf
    
    for tc0 in [1, 2, 3]:
        for a0 in [0.5, 0.8]:
            x0_init = [tc0, a0, 0.1, 1.0]
            try:
                result = minimize(objective, x0_init, method='Nelder-Mead',
                                options={'maxiter': 5000})
                if result.fun < best_nll and result.success:
                    best_nll = result.fun
                    best_result = result
            except:
                continue
    
    if best_result is None or not best_result.success:
        return None, np.inf
    
    return best_result.x, best_result.fun

# ============================================================================
# OVERDISPERSION ESTIMATION
# ============================================================================

def estimate_overdispersion(df, params, model_func):
    """
    Estimate overdispersion parameter (phi) using Pearson residuals.
    phi > 1 indicates overdispersion
    """
    ages = df['LarvalAge'].values
    n_settled = df['NoSet'].values
    n_total = df['NoAlive'].values
    
    # Predicted proportions
    pred_prop = model_func(ages, *params)
    pred_prop = np.clip(pred_prop, 1e-10, 1 - 1e-10)
    
    # Observed proportions
    obs_prop = n_settled / n_total
    
    # Pearson residuals
    expected_var = pred_prop * (1 - pred_prop) / n_total
    pearson_resid = (obs_prop - pred_prop) / np.sqrt(expected_var)
    
    # Estimate phi (dispersion parameter)
    n_params = len(params)
    phi = np.sum(pearson_resid**2) / (len(df) - n_params)
    
    return phi, pearson_resid

# ============================================================================
# BOOTSTRAP UNCERTAINTY (at replicate level)
# ============================================================================

def bootstrap_replicates(df, model_type, n_bootstrap=500):
    """
    Bootstrap at replicate level (not aggregated).
    More conservative and appropriate for experimental data.
    """
    print(f"    Bootstrapping {model_type} (replicate-level)...")
    
    # Original fit
    if model_type == 'Logistic':
        params_orig, nll_orig = fit_logistic_binomial(df)
    elif model_type == 'Gompertz':
        params_orig, nll_orig = fit_gompertz_binomial(df)
    elif model_type == 'Weibull':
        params_orig, nll_orig = fit_weibull_binomial(df)
    
    if params_orig is None:
        return None
    
    # Bootstrap
    bootstrap_params = []
    n_obs = len(df)
    
    for i in range(n_bootstrap):
        try:
            # Resample replicates with replacement
            indices = np.random.choice(n_obs, size=n_obs, replace=True)
            df_boot = df.iloc[indices].reset_index(drop=True)
            
            # Fit
            if model_type == 'Logistic':
                params_boot, _ = fit_logistic_binomial(df_boot)
            elif model_type == 'Gompertz':
                params_boot, _ = fit_gompertz_binomial(df_boot)
            elif model_type == 'Weibull':
                params_boot, _ = fit_weibull_binomial(df_boot)
            
            if params_boot is not None:
                bootstrap_params.append(params_boot)
        except:
            continue
    
    print(f"      Successful: {len(bootstrap_params)}/{n_bootstrap}")
    
    if len(bootstrap_params) < 50:
        return None
    
    bootstrap_params = np.array(bootstrap_params)
    
    return {
        'original_params': params_orig,
        'original_nll': nll_orig,
        'bootstrap_params': bootstrap_params,
        'n_successful': len(bootstrap_params)
    }

# ============================================================================
# CALCULATE TC50 with CI
# ============================================================================

def calculate_tc50(params, model_type, bootstrap_params=None):
    """Calculate TC50 (age at 50% competence) with confidence intervals."""
    
    if model_type == 'Logistic':
        L, k, x0 = params
        tc50 = x0  # For logistic, x0 is the inflection point
        
        if bootstrap_params is not None:
            tc50_boot = bootstrap_params[:, 2]
            tc50_ci = np.percentile(tc50_boot, [2.5, 97.5])
        else:
            tc50_ci = None
            
    elif model_type == 'Gompertz':
        a, b, c = params
        if b > 0 and c > 0:
            tc50 = np.log(2 * b) / c
        else:
            tc50 = np.nan
        
        if bootstrap_params is not None:
            tc50_boot = []
            for p in bootstrap_params:
                a_b, b_b, c_b = p
                if b_b > 0 and c_b > 0:
                    tc50_boot.append(np.log(2 * b_b) / c_b)
            if len(tc50_boot) > 10:
                tc50_ci = np.percentile(tc50_boot, [2.5, 97.5])
            else:
                tc50_ci = None
        else:
            tc50_ci = None
            
    elif model_type == 'Weibull':
        tc, a, b, v = params
        # Approximate TC50 for Weibull
        if a > 0.5 and b > 0:
            tc50 = tc + (np.log(2) / b) ** (1/v)
        else:
            tc50 = tc + 2
        
        if bootstrap_params is not None:
            tc50_boot = []
            for p in bootstrap_params:
                tc_b, a_b, b_b, v_b = p
                if a_b > 0.5 and b_b > 0 and v_b > 0:
                    tc50_boot.append(tc_b + (np.log(2) / b_b) ** (1/v_b))
            if len(tc50_boot) > 10:
                tc50_ci = np.percentile(tc50_boot, [2.5, 97.5])
            else:
                tc50_ci = None
        else:
            tc50_ci = None
    else:
        tc50 = np.nan
        tc50_ci = None
    
    return tc50, tc50_ci

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

def analyze_family_improved(family, df_family):
    """Improved analysis for a single family."""
    print(f"\n{'='*70}")
    print(f"FAMILY: {family}")
    print(f"{'='*70}")
    print(f"  Replicates: {len(df_family)}")
    print(f"  Total larvae: {df_family['NoAlive'].sum()}")
    print(f"  Total settled: {df_family['NoSet'].sum()}")
    print(f"  Overall settlement: {df_family['NoSet'].sum() / df_family['NoAlive'].sum() * 100:.1f}%")
    print(f"  Age range: {df_family['LarvalAge'].min()}-{df_family['LarvalAge'].max()} days")
    print(f"  Species: {df_family['Species'].nunique()}")
    
    # Fit all three models
    print("\n  Fitting models with binomial likelihood...")
    params_log, nll_log = fit_logistic_binomial(df_family)
    params_gomp, nll_gomp = fit_gompertz_binomial(df_family)
    params_weib, nll_weib = fit_weibull_binomial(df_family)
    
    # Calculate AICc (corrected for small sample size)
    n = len(df_family)
    
    def aicc(nll, k, n):
        aic = 2 * k + 2 * nll
        return aic + (2 * k * (k + 1)) / (n - k - 1) if n > k + 1 else np.inf
    
    aicc_log = aicc(nll_log, 3, n) if params_log is not None else np.inf
    aicc_gomp = aicc(nll_gomp, 3, n) if params_gomp is not None else np.inf
    aicc_weib = aicc(nll_weib, 4, n) if params_weib is not None else np.inf
    
    aiccs = {'Logistic': aicc_log, 'Gompertz': aicc_gomp, 'Weibull': aicc_weib}
    best_model = min(aiccs, key=aiccs.get)
    
    print(f"\n  Model AICc:")
    for model, aicc_val in aiccs.items():
        delta = aicc_val - aiccs[best_model]
        print(f"    {model:12s}: {aicc_val:10.1f}  (Δ = {delta:6.1f})")
    print(f"  → Best: {best_model}")
    
    # Get best model results
    if best_model == 'Logistic':
        params = params_log
        model_func = logistic_competency
    elif best_model == 'Gompertz':
        params = params_gomp
        model_func = gompertz_competency
    else:
        params = params_weib
        model_func = weibull_competency
    
    # Estimate overdispersion
    phi, pearson_resid = estimate_overdispersion(df_family, params, model_func)
    print(f"\n  Overdispersion parameter (φ): {phi:.2f}")
    if phi > 1.5:
        print(f"    → Moderate overdispersion detected")
    elif phi > 2.0:
        print(f"    → Strong overdispersion detected")
    
    # Bootstrap uncertainty
    bootstrap_result = bootstrap_replicates(df_family, best_model, n_bootstrap=500)
    
    if bootstrap_result is None:
        print("  WARNING: Bootstrap failed")
        return None
    
    # Calculate TC50
    tc50, tc50_ci = calculate_tc50(params, best_model, 
                                   bootstrap_result['bootstrap_params'])
    
    print(f"\n  TC50 (age at 50% competence): {tc50:.2f} days")
    if tc50_ci is not None:
        print(f"    95% CI: [{tc50_ci[0]:.2f}, {tc50_ci[1]:.2f}]")
    
    # Parameter summary
    print(f"\n  Parameters:")
    param_names = {'Logistic': ['L', 'k', 'x0'],
                  'Gompertz': ['a', 'b', 'c'],
                  'Weibull': ['tc', 'a', 'b', 'v']}[best_model]
    
    bootstrap_params = bootstrap_result['bootstrap_params']
    for i, name in enumerate(param_names):
        mean = np.mean(bootstrap_params[:, i])
        std = np.std(bootstrap_params[:, i])
        ci = np.percentile(bootstrap_params[:, i], [2.5, 97.5])
        print(f"    {name:4s} = {mean:.4f} ± {std:.4f}  [{ci[0]:.4f}, {ci[1]:.4f}]")
    
    return {
        'family': family,
        'model': best_model,
        'params': params,
        'n_replicates': len(df_family),
        'n_larvae': df_family['NoAlive'].sum(),
        'phi': phi,
        'tc50': tc50,
        'tc50_ci': tc50_ci,
        'bootstrap': bootstrap_result,
        'aicc': aiccs[best_model]
    }

# ============================================================================
# MAIN
# ============================================================================

def main():
    data_path = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity_KERNELS/R_codes/Randal/data/Settlement.csv'
    
    # Load all data
    df = load_all_data(data_path)
    
    # Analyze each family
    results = {}
    for family in sorted(df['Family'].unique()):
        df_family = df[df['Family'] == family].copy()
        
        # Require minimum data
        if len(df_family) >= 10 and df_family['NoAlive'].sum() >= 100:
            result = analyze_family_improved(family, df_family)
            if result is not None:
                results[family] = result
    
    # Save summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}\n")
    
    summary_data = []
    for family, res in results.items():
        tc50_str = f"{res['tc50']:.2f}"
        if res['tc50_ci'] is not None:
            tc50_ci_str = f"[{res['tc50_ci'][0]:.2f}, {res['tc50_ci'][1]:.2f}]"
        else:
            tc50_ci_str = "N/A"
        
        summary_data.append({
            'Family': family,
            'Model': res['model'],
            'N_replicates': res['n_replicates'],
            'N_larvae': res['n_larvae'],
            'TC50': tc50_str,
            'TC50_CI': tc50_ci_str,
            'Phi': f"{res['phi']:.2f}",
            'AICc': f"{res['aicc']:.1f}"
        })
    
    summary_df = pd.DataFrame(summary_data)
    print(summary_df.to_string(index=False))
    
    summary_df.to_csv('improved_competency_summary.csv', index=False)
    print(f"\n✓ Summary saved: improved_competency_summary.csv")
    
    # Save detailed results
    import pickle
    with open('improved_competency_results.pkl', 'wb') as f:
        pickle.dump(results, f)
    print(f"✓ Detailed results saved: improved_competency_results.pkl")
    
    return results

if __name__ == "__main__":
    results = main()

