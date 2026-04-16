#!/usr/bin/env python3
"""
Monegetti et al. Piecewise Weibull-Exponential Competency Model
================================================================

Based on:
Moneghetti et al. (2019) "High-frequency sampling and piecewise models 
reshape dispersal kernels of a common reef coral"

Model Description:
- Precompetency period (t < tc): No settlement possible
- Early competency (tc < t < Tcp): Weibull loss of competency
- Late competency (t > Tcp): Exponential loss of competency

Competency is calculated through integration:
- Probability of becoming competent at time tau
- Probability of maintaining competency from tau to t
- Integrated over all possible tau values

Parameters:
- a: Rate of acquisition of competency
- b1, v1: Weibull loss parameters (early period)
- b2: Exponential loss parameter (late period)
- tc: Precompetency period (onset)
- Tcp: Change point between early and late phases
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize, differential_evolution
from scipy.integrate import quad
import warnings
import os
from datetime import datetime
import re
warnings.filterwarnings('ignore')

sns.set_style("whitegrid")

# ============================================================================
# MONEGETTI PIECEWISE MODEL
# ============================================================================

def monegetti_competency(age, a, b1, v1, b2, tc, Tcp):
    """
    Monegetti et al. piecewise Weibull-exponential competency model.
    
    Parameters:
    -----------
    age : float or array
        Larval age in days
    a : float
        Rate of acquisition of competency (when t > tc)
    b1 : float
        Weibull loss parameter (early period, scale)
    v1 : float
        Weibull loss parameter (early period, shape)
    b2 : float
        Exponential loss parameter (late period)
    tc : float
        Precompetency period (onset age)
    Tcp : float
        Change point between early and late phases
        
    Returns:
    --------
    competency : float or array
        Probability of being competent at given age
    """
    age = np.atleast_1d(age)
    competency = np.zeros_like(age, dtype=float)
    
    for i, t in enumerate(age):
        if t < tc:
            # Precompetent
            competency[i] = 0.0
            
        elif t <= Tcp:
            # Early phase: integrate Weibull loss
            def integrand_early(tau):
                return a * np.exp(-a * (tau - tc)) * np.exp(-((b1 * (t - tau))**v1))
            
            try:
                result, error = quad(integrand_early, tc, t, limit=50, epsabs=1e-4, epsrel=1e-4)
                competency[i] = max(0, min(1, result))
            except:
                competency[i] = 0.0
                
        else:  # t > Tcp
            # Late phase: integrate both phases
            # Part 1: Early phase up to Tcp
            def integrand_latter_1(tau):
                return (a * np.exp(-a * (tau - tc)) * 
                       np.exp(-((b1 * (Tcp - tau))**v1)) * 
                       np.exp(-b2 * (t - Tcp)))
            
            # Part 2: Late phase from Tcp to t
            def integrand_latter_2(tau):
                return a * np.exp(-a * (tau - tc)) * np.exp(-b2 * (t - tau))
            
            try:
                result1, _ = quad(integrand_latter_1, tc, Tcp, limit=50, epsabs=1e-4, epsrel=1e-4)
                result2, _ = quad(integrand_latter_2, Tcp, t, limit=50, epsabs=1e-4, epsrel=1e-4)
                competency[i] = max(0, min(1, result1 + result2))
            except:
                competency[i] = 0.0
    
    return competency if len(competency) > 1 else competency[0]

# ============================================================================
# BINOMIAL LIKELIHOOD FITTING
# ============================================================================

def neg_log_likelihood_monegetti(params, ages, n_settled, n_total):
    """
    Negative log-likelihood for Monegetti model with binomial data.
    """
    a, b1, v1, b2, tc, Tcp = params
    
    # Predict competency for each observation
    pred_comp = np.zeros_like(ages, dtype=float)
    for i, age in enumerate(ages):
        pred_comp[i] = monegetti_competency(age, a, b1, v1, b2, tc, Tcp)
    
    # Clip to avoid log(0)
    pred_comp = np.clip(pred_comp, 1e-10, 1 - 1e-10)
    
    # Binomial log-likelihood
    ll = np.sum(n_settled * np.log(pred_comp) + 
                (n_total - n_settled) * np.log(1 - pred_comp))
    
    return -ll

def fit_monegetti_model(df):
    """
    Fit Monegetti piecewise model to settlement data.
    
    Parameters:
    -----------
    df : DataFrame
        Must have columns: LarvalAge, NoSet, NoAlive
        
    Returns:
    --------
    params : array or None
        [a, b1, v1, b2, tc, Tcp] if successful, None otherwise
    nll : float
        Negative log-likelihood
    """
    ages = df['LarvalAge'].values
    n_settled = df['NoSet'].values
    n_total = df['NoAlive'].values
    
    # Determine constraints from data
    # tc: must be before first settlement
    first_settlement = df[df['NoSet'] > 0]['LarvalAge'].min()
    tc_min = max(0, first_settlement - 2)
    tc_max = first_settlement
    
    # Tcp: must be after tc and before max age
    max_age = ages.max()
    Tcp_min = first_settlement
    Tcp_max = max_age
    
    # Parameter bounds
    # a: rate of acquisition (0.1 to 5)
    # b1: Weibull scale (0.001 to 1)
    # v1: Weibull shape (0.1 to 5)
    # b2: Exponential rate (0.001 to 1)
    # tc: precompetency (tc_min to tc_max)
    # Tcp: change point (Tcp_min to Tcp_max)
    bounds = [
        (0.1, 5.0),      # a
        (0.001, 1.0),    # b1
        (0.1, 5.0),      # v1
        (0.001, 1.0),    # b2
        (tc_min, tc_max), # tc
        (Tcp_min, Tcp_max) # Tcp
    ]
    
    print(f"    Constraints: tc ∈ [{tc_min:.1f}, {tc_max:.1f}], Tcp ∈ [{Tcp_min:.1f}, {Tcp_max:.1f}]")
    
    def objective(params):
        # Additional constraint: tc < Tcp
        if params[4] >= params[5]:
            return 1e10
        return neg_log_likelihood_monegetti(params, ages, n_settled, n_total)
    
    # Try local optimization with multiple starting points (faster than global)
    print(f"    Running optimization with multiple starting points...")
    
    # Generate starting points
    n_starts = 5
    best_result = None
    best_nll = np.inf
    
    # Create starting points
    starts = []
    for _ in range(n_starts):
        start = [
            np.random.uniform(0.5, 2.0),      # a
            np.random.uniform(0.01, 0.5),    # b1
            np.random.uniform(0.5, 2.0),     # v1
            np.random.uniform(0.01, 0.5),    # b2
            np.random.uniform(tc_min, tc_max),  # tc
            np.random.uniform(Tcp_min, Tcp_max)  # Tcp
        ]
        # Ensure tc < Tcp
        if start[4] >= start[5]:
            start[5] = start[4] + 1
        starts.append(start)
    
    for i, start in enumerate(starts):
        try:
            result = minimize(
                objective,
                start,
                method='L-BFGS-B',
                bounds=bounds,
                options={'maxiter': 50, 'ftol': 1e-3}
            )
            
            if result.fun < best_nll and result.success:
                best_nll = result.fun
                best_result = result
                print(f"    Start {i+1}/{n_starts}: NLL = {result.fun:.2f} ✓")
            else:
                print(f"    Start {i+1}/{n_starts}: Failed")
        except Exception as e:
            print(f"    Start {i+1}/{n_starts}: Error - {str(e)}")
            continue
    
    if best_result is not None and best_nll < 1e9:
        print(f"    ✓ Best result: NLL = {best_nll:.2f}")
        return best_result.x, best_nll
    else:
        print(f"    ✗ All optimizations failed")
        return None, np.inf

# ============================================================================
# ANALYSIS FUNCTIONS
# ============================================================================

def analyze_family_monegetti(family, df_family):
    """
    Fit Monegetti piecewise model to a single family.
    """
    print(f"\n{'='*70}")
    print(f"FAMILY: {family}")
    print(f"{'='*70}")
    print(f"  Replicates: {len(df_family)}")
    print(f"  Age range: {df_family['LarvalAge'].min()}-{df_family['LarvalAge'].max()} days")
    print(f"  Total larvae: {df_family['NoAlive'].sum()}")
    print(f"  Total settled: {df_family['NoSet'].sum()}")
    
    # Fit Monegetti model
    print(f"\n  Fitting Monegetti piecewise model (6 parameters)...")
    params, nll = fit_monegetti_model(df_family)
    
    if params is None:
        print(f"  WARNING: Model fitting failed")
        return None
    
    a, b1, v1, b2, tc, Tcp = params
    
    print(f"\n  Parameters:")
    print(f"    a   (acquisition rate)   = {a:.4f}")
    print(f"    b1  (Weibull scale)      = {b1:.4f}")
    print(f"    v1  (Weibull shape)      = {v1:.4f}")
    print(f"    b2  (exponential rate)   = {b2:.4f}")
    print(f"    tc  (precompetency)      = {tc:.4f} days")
    print(f"    Tcp (change point)       = {Tcp:.4f} days")
    
    # Calculate AICc
    n = len(df_family)
    k = 6  # number of parameters
    aicc = 2 * k + 2 * nll + (2 * k * (k + 1)) / (n - k - 1) if n > k + 1 else np.inf
    
    print(f"\n  AICc: {aicc:.1f}")
    print(f"  Negative Log-Likelihood: {nll:.1f}")
    
    # Calculate TC50 (approximate)
    # Find age where competency reaches 50% of maximum
    test_ages = np.linspace(tc, df_family['LarvalAge'].max(), 100)
    pred_comp = np.array([monegetti_competency(age, a, b1, v1, b2, tc, Tcp) 
                         for age in test_ages])
    max_comp = pred_comp.max()
    
    if max_comp > 0.1:
        idx_50 = np.argmin(np.abs(pred_comp - 0.5 * max_comp))
        tc50 = test_ages[idx_50]
        print(f"\n  TC50 (age at 50% of max competency): {tc50:.2f} days")
        print(f"  Maximum competency: {max_comp:.3f}")
    else:
        tc50 = np.nan
        print(f"\n  WARNING: Maximum competency too low ({max_comp:.3f})")
    
    return {
        'family': family,
        'model': 'Monegetti',
        'params': params,
        'param_names': ['a', 'b1', 'v1', 'b2', 'tc', 'Tcp'],
        'n_replicates': len(df_family),
        'n_larvae': df_family['NoAlive'].sum(),
        'nll': nll,
        'aicc': aicc,
        'tc50': tc50,
        'max_competency': pred_comp.max() if len(pred_comp) > 0 else np.nan
    }

# ============================================================================
# COMPARISON WITH SIMPLE MODELS
# ============================================================================

def compare_with_simple_models(family, df_family, monegetti_result):
    """
    Compare Monegetti model with simple Weibull/Logistic models.
    """
    from improved_competency_model import (
        fit_logistic_binomial, fit_gompertz_binomial, fit_weibull_binomial
    )
    
    print(f"\n  Comparing with simple models...")
    
    # Fit simple models
    params_log, nll_log = fit_logistic_binomial(df_family)
    params_gomp, nll_gomp = fit_gompertz_binomial(df_family)
    params_weib, nll_weib = fit_weibull_binomial(df_family)
    
    n = len(df_family)
    
    def aicc(nll, k, n):
        aic = 2 * k + 2 * nll
        return aic + (2 * k * (k + 1)) / (n - k - 1) if n > k + 1 else np.inf
    
    aicc_log = aicc(nll_log, 3, n) if params_log is not None else np.inf
    aicc_gomp = aicc(nll_gomp, 3, n) if params_gomp is not None else np.inf
    aicc_weib = aicc(nll_weib, 4, n) if params_weib is not None else np.inf
    aicc_mon = monegetti_result['aicc']
    
    print(f"\n  Model Comparison (AICc):")
    print(f"    Logistic (3 params): {aicc_log:10.1f}")
    print(f"    Gompertz (3 params): {aicc_gomp:10.1f}")
    print(f"    Weibull  (4 params): {aicc_weib:10.1f}")
    print(f"    Monegetti (6 params): {aicc_mon:10.1f}")
    
    # Best model
    models = {'Logistic': aicc_log, 'Gompertz': aicc_gomp, 
              'Weibull': aicc_weib, 'Monegetti': aicc_mon}
    best_model = min(models, key=models.get)
    delta_aicc = aicc_mon - models[best_model]
    
    print(f"\n  → Best model: {best_model}")
    print(f"  → Monegetti Δ AICc: {delta_aicc:+.1f}")
    
    if delta_aicc < -10:
        print(f"  → Monegetti is STRONGLY preferred (Δ AICc < -10)")
    elif delta_aicc < -2:
        print(f"  → Monegetti is preferred (Δ AICc < -2)")
    elif delta_aicc < 2:
        print(f"  → Models are equivalent (|Δ AICc| < 2)")
    else:
        print(f"  → Simpler model preferred (Δ AICc > 2)")
    
    return {
        'aicc_logistic': aicc_log,
        'aicc_gompertz': aicc_gomp,
        'aicc_weibull': aicc_weib,
        'aicc_monegetti': aicc_mon,
        'best_model': best_model,
        'delta_aicc': delta_aicc
    }

# ============================================================================
# FILE UTILITIES
# ============================================================================

def sanitize_filename(name):
    """
    Sanitize family name for use in filename.
    Removes or replaces special characters that are problematic in filenames.
    
    Parameters:
    -----------
    name : str
        Family name to sanitize
        
    Returns:
    --------
    sanitized : str
        Sanitized name safe for use in filenames
    """
    # Replace spaces and special characters with underscores
    sanitized = re.sub(r'[^\w\s-]', '', str(name))
    sanitized = re.sub(r'[-\s]+', '_', sanitized)
    # Remove leading/trailing underscores
    sanitized = sanitized.strip('_')
    return sanitized

def save_family_result(family, result, comparison, output_dir='estimations'):
    """
    Save individual family result to CSV file with timestamp.
    
    Parameters:
    -----------
    family : str
        Family name
    result : dict
        Result dictionary from analyze_family_monegetti
    comparison : dict
        Comparison dictionary from compare_with_simple_models
    output_dir : str
        Directory to save results (default: 'estimations')
        
    Returns:
    --------
    filename : str
        Path to saved file
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Sanitize family name for filename
    sanitized_family = sanitize_filename(family)
    
    # Generate timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    # Create filename
    filename = os.path.join(output_dir, f'estimations_{sanitized_family}_{timestamp}.csv')
    
    # Prepare data for CSV
    if result is not None:
        a, b1, v1, b2, tc, Tcp = result['params']
        
        data = {
            'Family': [family],
            'Timestamp': [timestamp],
            'Model': [result['model']],
            # Parameters
            'a': [a],
            'b1': [b1],
            'v1': [v1],
            'b2': [b2],
            'tc': [tc],
            'Tcp': [Tcp],
            # Metrics
            'N_replicates': [result['n_replicates']],
            'N_larvae': [result['n_larvae']],
            'NLL': [result['nll']],
            'AICc': [result['aicc']],
            'TC50': [result['tc50'] if not np.isnan(result['tc50']) else np.nan],
            'Max_competency': [result['max_competency']],
        }
        
        # Add comparison results if available
        if comparison:
            data['AICc_Logistic'] = [comparison.get('aicc_logistic', np.nan)]
            data['AICc_Gompertz'] = [comparison.get('aicc_gompertz', np.nan)]
            data['AICc_Weibull'] = [comparison.get('aicc_weibull', np.nan)]
            data['AICc_Monegetti'] = [comparison.get('aicc_monegetti', np.nan)]
            data['Best_model'] = [comparison.get('best_model', 'N/A')]
            data['Delta_AICc'] = [comparison.get('delta_aicc', np.nan)]
        else:
            data['AICc_Logistic'] = [np.nan]
            data['AICc_Gompertz'] = [np.nan]
            data['AICc_Weibull'] = [np.nan]
            data['AICc_Monegetti'] = [np.nan]
            data['Best_model'] = ['N/A']
            data['Delta_AICc'] = [np.nan]
    else:
        # If result is None (fitting failed), save what we can
        data = {
            'Family': [family],
            'Timestamp': [timestamp],
            'Model': ['Monegetti'],
            'Status': ['FAILED'],
            'a': [np.nan],
            'b1': [np.nan],
            'v1': [np.nan],
            'b2': [np.nan],
            'tc': [np.nan],
            'Tcp': [np.nan],
            'N_replicates': [np.nan],
            'N_larvae': [np.nan],
            'NLL': [np.nan],
            'AICc': [np.nan],
            'TC50': [np.nan],
            'Max_competency': [np.nan],
            'AICc_Logistic': [np.nan],
            'AICc_Gompertz': [np.nan],
            'AICc_Weibull': [np.nan],
            'AICc_Monegetti': [np.nan],
            'Best_model': ['N/A'],
            'Delta_AICc': [np.nan],
        }
    
    # Create DataFrame and save
    df_result = pd.DataFrame(data)
    df_result.to_csv(filename, index=False)
    
    return filename

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

def main(test_mode=False):
    """
    Fit Monegetti piecewise model to all families.
    
    Parameters:
    -----------
    test_mode : bool
        If True, only analyze first family (for testing)
    """
    from improved_competency_model import load_all_data
    
    print("="*70)
    print("MONEGETTI PIECEWISE WEIBULL-EXPONENTIAL MODEL")
    print("="*70)
    if test_mode:
        print("TEST MODE: Analyzing first family only")
    print("="*70)
    
    # Load data
    data_path = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity_KERNELS/R_codes/Randal/data/Settlement.csv'
    df = load_all_data(data_path)
    
    # Analyze each family
    results = {}
    comparisons = {}
    
    families_to_analyze = sorted(df['Family'].unique())
    if test_mode:
        families_to_analyze = families_to_analyze[:1]
        print(f"\nTEST MODE: Analyzing only {families_to_analyze[0]}")
    
    total_families = len(families_to_analyze)
    print(f"\nTotal families to analyze: {total_families}")
    print(f"Results will be saved to: estimations/estimations_*Family*_*timestamp*.csv\n")
    
    for idx, family in enumerate(families_to_analyze, 1):
        print(f"\n[{idx}/{total_families}] Processing family: {family}")
        
        df_family = df[df['Family'] == family].copy()
        
        # Require minimum data (reduced for testing)
        if len(df_family) >= 30 and df_family['NoAlive'].sum() >= 200:
            result = analyze_family_monegetti(family, df_family)
            
            comparison = None
            if result is not None:
                results[family] = result
                
                # Compare with simple models
                try:
                    comparison = compare_with_simple_models(family, df_family, result)
                    comparisons[family] = comparison
                except Exception as e:
                    print(f"  WARNING: Model comparison failed: {str(e)}")
                    comparison = None
            
            # Save result immediately after processing this family
            try:
                filename = save_family_result(family, result, comparison)
                print(f"\n  ✓ Saved: {filename}")
            except Exception as e:
                print(f"\n  ✗ ERROR saving result: {str(e)}")
        else:
            print(f"  SKIPPED: Insufficient data (replicates: {len(df_family)}, larvae: {df_family['NoAlive'].sum()})")
            # Save a record even for skipped families
            try:
                filename = save_family_result(family, None, None)
                print(f"  ✓ Saved skip record: {filename}")
            except Exception as e:
                print(f"  ✗ ERROR saving skip record: {str(e)}")
    
    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}\n")
    
    if len(results) > 0:
        summary_data = []
        for family, res in results.items():
            comp = comparisons.get(family, {})
            
            summary_data.append({
                'Family': family,
                'N': res['n_replicates'],
                'tc': f"{res['params'][4]:.2f}",
                'Tcp': f"{res['params'][5]:.2f}",
                'TC50': f"{res['tc50']:.2f}" if not np.isnan(res['tc50']) else "N/A",
                'Max_Comp': f"{res['max_competency']:.3f}",
                'AICc': f"{res['aicc']:.1f}",
                'Best_Simple': comp.get('best_model', 'N/A'),
                'Δ_AICc': f"{comp.get('delta_aicc', np.nan):+.1f}"
            })
        
        summary_df = pd.DataFrame(summary_data)
        print(summary_df.to_string(index=False))
        
        # Save results
        summary_df.to_csv('monegetti_model_summary.csv', index=False)
        print(f"\n✓ Summary saved: monegetti_model_summary.csv")
        
        import pickle
        with open('monegetti_model_results.pkl', 'wb') as f:
            pickle.dump({'results': results, 'comparisons': comparisons}, f)
        print(f"✓ Detailed results saved: monegetti_model_results.pkl")
    else:
        print("No families were successfully analyzed.")
    
    print(f"\n✓ Individual family results saved in: estimations/")
    print(f"  Format: estimations_*Family*_*timestamp*.csv")
    
    return results, comparisons

if __name__ == "__main__":
    import sys
    test_mode = '--test' in sys.argv
    results, comparisons = main(test_mode=test_mode)

