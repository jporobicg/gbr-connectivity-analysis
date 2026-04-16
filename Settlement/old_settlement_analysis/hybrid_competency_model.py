#!/usr/bin/env python3
"""
Hybrid Competency Model: Continuous Proportion + Monegetti-Style Fitting
=========================================================================

This approach combines:
1. Randal's data processing (rubble treatment, family-level)
2. Continuous proportion settled (NOT binary cumulative)
3. Monegetti-style model fitting (Weibull, Logistic, Gompertz)

Key advantage: Avoids binary threshold artifacts while allowing natural
fluctuations in settlement proportion over time.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize, differential_evolution
from scipy.integrate import quad
import warnings
warnings.filterwarnings('ignore')

sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 10)

# ============================================================================
# DATA PREPARATION
# ============================================================================

def prepare_data_continuous(data_path, treatment='rubble'):
    """
    Prepare data using CONTINUOUS proportions (not binary cumulative).
    
    Key differences from binary approach:
    - Uses actual proportion settled at each age
    - Allows proportion to fluctuate (increase or decrease)
    - Aggregates to family level but keeps age structure
    """
    df = pd.read_csv(data_path)
    
    print("=" * 70)
    print("HYBRID APPROACH: Continuous Proportion + Monegetti Models")
    print("=" * 70)
    print(f"Total observations: {len(df)}")
    
    # Rename if needed
    if 'Treatment' in df.columns and 'SpecificTreatment' not in df.columns:
        df['SpecificTreatment'] = df['Treatment']
    
    # Filter to rubble treatment
    df_rubble = df[df['SpecificTreatment'] == treatment].copy()
    print(f"After filtering to {treatment}: {len(df_rubble)}")
    
    # Calculate proportion settled (continuous, not binary)
    df_rubble['PropSettled'] = df_rubble['NoSet'] / (df_rubble['NoSet'] + df_rubble['NoNotSet'])
    
    # Aggregate by family and age
    family_data = {}
    for family in df_rubble['Family'].dropna().unique():
        df_fam = df_rubble[df_rubble['Family'] == family].copy()
        
        # Aggregate by age - sum larvae, calculate proportion
        agg = df_fam.groupby('LarvalAge').agg({
            'NoSet': 'sum',
            'NoNotSet': 'sum',
            'NoAlive': 'sum'
        }).reset_index()
        
        # Calculate continuous proportion
        agg['PropCompetent'] = agg['NoSet'] / (agg['NoSet'] + agg['NoNotSet'])
        agg['Family'] = family
        agg['N_Species'] = df_fam['Species'].nunique()
        
        # Only include families with sufficient data
        if len(agg) >= 5 and agg['NoAlive'].sum() > 50:
            family_data[family] = agg
    
    print(f"\nFamilies with sufficient data: {len(family_data)}")
    print(f"Families: {list(family_data.keys())}")
    
    return family_data

# ============================================================================
# MONEGETTI-STYLE MODELS (from first analysis)
# ============================================================================

def logistic_competency(age, L, k, x0):
    """Logistic (sigmoidal) model."""
    return L / (1 + np.exp(-k * (age - x0)))

def gompertz_competency(age, a, b, c):
    """Gompertz model."""
    return a * np.exp(-b * np.exp(-c * age))

def weibull_competency(age, tc, a, b, v):
    """Simplified Weibull model."""
    age = np.atleast_1d(age)
    competency = np.zeros_like(age, dtype=float)
    
    for i, t in enumerate(age):
        if t < tc:
            competency[i] = 0
        else:
            def integrand(tau):
                return a * np.exp(-a * (tau - tc)) * np.exp(-((b * (t - tau)) ** v))
            
            try:
                result, _ = quad(integrand, tc, t)
                competency[i] = np.clip(result, 0, 1)
            except:
                competency[i] = 0
    
    return competency

# ============================================================================
# MODEL FITTING TO CONTINUOUS PROPORTIONS
# ============================================================================

def fit_logistic_continuous(df_agg, family_name=""):
    """
    Fit logistic model to continuous proportion data.
    Uses weighted least squares (weight by sample size).
    """
    def neg_weighted_likelihood(params):
        L, k, x0 = params
        if L <= 0 or L > 1 or k <= 0:
            return 1e10
        
        pred = logistic_competency(df_agg['LarvalAge'].values, L, k, x0)
        pred = np.clip(pred, 1e-10, 1 - 1e-10)
        
        # Weighted binomial likelihood
        weights = df_agg['NoSet'] + df_agg['NoNotSet']  # Sample size as weight
        ll = np.sum(weights * (
            df_agg['NoSet'].values * np.log(pred) + 
            df_agg['NoNotSet'].values * np.log(1 - pred)
        ))
        
        return -ll
    
    # Initial guess
    max_comp = df_agg['PropCompetent'].max()
    x0_init = df_agg.loc[df_agg['PropCompetent'] >= max_comp/2, 'LarvalAge'].min()
    if np.isnan(x0_init):
        x0_init = df_agg['LarvalAge'].median()
    
    initial = [max(max_comp, 0.5), 0.2, x0_init]
    
    try:
        result = minimize(neg_weighted_likelihood, initial, method='Nelder-Mead',
                         options={'maxiter': 5000})
        
        if family_name:
            L_opt, k_opt, x0_opt = result.x
            print(f"  Logistic: L={L_opt:.3f}, k={k_opt:.3f}, x0={x0_opt:.2f}, AIC={2*3 + 2*result.fun:.1f}")
        
        return result.x, result.fun
    except:
        return None, 1e10

def fit_gompertz_continuous(df_agg, family_name=""):
    """Fit Gompertz model to continuous proportion data."""
    def neg_weighted_likelihood(params):
        a, b, c = params
        if a <= 0 or a > 1 or b <= 0 or c <= 0:
            return 1e10
        
        pred = gompertz_competency(df_agg['LarvalAge'].values, a, b, c)
        pred = np.clip(pred, 1e-10, 1 - 1e-10)
        
        weights = df_agg['NoSet'] + df_agg['NoNotSet']
        ll = np.sum(weights * (
            df_agg['NoSet'].values * np.log(pred) + 
            df_agg['NoNotSet'].values * np.log(1 - pred)
        ))
        
        return -ll
    
    max_comp = df_agg['PropCompetent'].max()
    initial = [max(max_comp, 0.5), 2.0, 0.15]
    
    try:
        result = minimize(neg_weighted_likelihood, initial, method='Nelder-Mead',
                         options={'maxiter': 5000})
        
        if family_name:
            a_opt, b_opt, c_opt = result.x
            print(f"  Gompertz: a={a_opt:.3f}, b={b_opt:.3f}, c={c_opt:.3f}, AIC={2*3 + 2*result.fun:.1f}")
        
        return result.x, result.fun
    except:
        return None, 1e10

def fit_weibull_continuous(df_agg, family_name=""):
    """Fit Weibull model to continuous proportion data."""
    def neg_weighted_likelihood(params):
        tc, a, b, v = np.exp(params)
        
        min_age = df_agg['LarvalAge'].min()
        if tc >= min_age:
            return 1e10
        
        pred = weibull_competency(df_agg['LarvalAge'].values, tc, a, b, v)
        pred = np.clip(pred, 1e-10, 1 - 1e-10)
        
        weights = df_agg['NoSet'] + df_agg['NoNotSet']
        ll = np.sum(weights * (
            df_agg['NoSet'].values * np.log(pred) + 
            df_agg['NoNotSet'].values * np.log(1 - pred)
        ))
        
        return -ll
    
    min_age = df_agg['LarvalAge'].min()
    initial = [np.log(max(1.0, min_age - 1)), np.log(1.0), np.log(0.01), np.log(0.5)]
    
    try:
        result = minimize(neg_weighted_likelihood, initial, method='Nelder-Mead',
                         options={'maxiter': 5000})
        params_opt = np.exp(result.x)
        
        if family_name:
            print(f"  Weibull: tc={params_opt[0]:.2f}, a={params_opt[1]:.3f}, b={params_opt[2]:.3f}, v={params_opt[3]:.3f}, AIC={2*4 + 2*result.fun:.1f}")
        
        return params_opt, result.fun
    except:
        return None, 1e10

# ============================================================================
# MONEGETTI MODEL FOR COMPARISON
# ============================================================================

def monegetti_model(age):
    """Monegetti's fitted model."""
    tc = 3.332
    Tcp = 69.91
    a = 1.292
    b1 = 0.001878
    v1 = 0.3645
    b2 = 0.3969
    
    age = np.atleast_1d(age)
    competency = np.zeros_like(age, dtype=float)
    
    for i, t in enumerate(age):
        if t < tc:
            competency[i] = 0
        elif t <= Tcp:
            def int_early(tau):
                return a * np.exp(-a * (tau - tc)) * np.exp(-((b1 * (t - tau)) ** v1))
            try:
                result, _ = quad(int_early, tc, t)
                competency[i] = np.clip(result, 0, 1)
            except:
                competency[i] = 0
        else:
            def int_latter_1(tau):
                return a * np.exp(-a * (tau - tc)) * np.exp(-((b1 * (Tcp - tau)) ** v1)) * np.exp(-b2 * (t - Tcp))
            def int_latter_2(tau):
                return a * np.exp(-a * (tau - tc)) * np.exp(-b2 * (t - tau))
            try:
                result1, _ = quad(int_latter_1, tc, Tcp)
                result2, _ = quad(int_latter_2, Tcp, t)
                competency[i] = np.clip(result1 + result2, 0, 1)
            except:
                competency[i] = 0
    
    return competency

# ============================================================================
# ANALYSIS
# ============================================================================

def analyze_families(family_data):
    """Fit all three models to each family."""
    results = {}
    
    print(f"\n{'=' * 70}")
    print("FITTING MODELS TO CONTINUOUS PROPORTIONS")
    print(f"{'=' * 70}\n")
    
    for family, df_agg in family_data.items():
        print(f"\n{family}")
        print(f"  Observations: {len(df_agg)}")
        print(f"  Age range: {df_agg['LarvalAge'].min():.0f}-{df_agg['LarvalAge'].max():.0f} days")
        print(f"  Max proportion: {df_agg['PropCompetent'].max():.3f}")
        
        # Fit all three models
        params_log, ll_log = fit_logistic_continuous(df_agg, family)
        params_gomp, ll_gomp = fit_gompertz_continuous(df_agg, family)
        params_weib, ll_weib = fit_weibull_continuous(df_agg, family)
        
        # Calculate AICs
        aic_log = 2 * 3 + 2 * ll_log if params_log is not None else np.inf
        aic_gomp = 2 * 3 + 2 * ll_gomp if params_gomp is not None else np.inf
        aic_weib = 2 * 4 + 2 * ll_weib if params_weib is not None else np.inf
        
        # Determine best
        aics = {'Logistic': aic_log, 'Gompertz': aic_gomp, 'Weibull': aic_weib}
        best_model = min(aics, key=aics.get)
        
        print(f"  → Best model: {best_model}")
        
        results[family] = {
            'data': df_agg,
            'models': {
                'Logistic': {'params': params_log, 'll': ll_log, 'aic': aic_log},
                'Gompertz': {'params': params_gomp, 'll': ll_gomp, 'aic': aic_gomp},
                'Weibull': {'params': params_weib, 'll': ll_weib, 'aic': aic_weib}
            },
            'best_model': best_model
        }
    
    return results

# ============================================================================
# VISUALIZATION
# ============================================================================

def plot_family_results(results, save_path):
    """Create comprehensive visualization."""
    n_families = len(results)
    ncols = min(3, n_families)
    nrows = int(np.ceil(n_families / ncols))
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows))
    if n_families == 1:
        axes = np.array([axes])
    axes = axes.flatten() if n_families > 1 else axes
    
    fig.suptitle('Hybrid Approach: Continuous Proportion + Monegetti Models\n(Rubble Treatment)', 
                 fontsize=16, fontweight='bold', y=0.995)
    
    age_range = np.linspace(0, 80, 400)
    colors = {'Logistic': '#e74c3c', 'Gompertz': '#3498db', 'Weibull': '#2ecc71'}
    
    for idx, (family, result) in enumerate(results.items()):
        if idx >= len(axes):
            break
        
        ax = axes[idx]
        df_agg = result['data']
        best_model = result['best_model']
        
        # Plot observed proportions (size by sample size)
        sizes = (df_agg['NoSet'] + df_agg['NoNotSet']) * 2
        ax.scatter(df_agg['LarvalAge'], df_agg['PropCompetent'],
                  s=sizes, alpha=0.6, color='black', edgecolors='white',
                  linewidth=1.5, zorder=5, label='Observed')
        
        # Plot best model (thick line)
        params = result['models'][best_model]['params']
        if params is not None:
            if best_model == 'Logistic':
                pred = logistic_competency(age_range, *params)
            elif best_model == 'Gompertz':
                pred = gompertz_competency(age_range, *params)
            elif best_model == 'Weibull':
                pred = weibull_competency(age_range, *params)
            
            ax.plot(age_range, pred, linewidth=3.5, 
                   label=f'{best_model} (best)', 
                   color=colors[best_model], zorder=4)
        
        # Plot other models (thin lines)
        for model_name, model_data in result['models'].items():
            if model_name != best_model and model_data['params'] is not None:
                params = model_data['params']
                if model_name == 'Logistic':
                    pred = logistic_competency(age_range, *params)
                elif model_name == 'Gompertz':
                    pred = gompertz_competency(age_range, *params)
                elif model_name == 'Weibull':
                    pred = weibull_competency(age_range, *params)
                
                ax.plot(age_range, pred, linewidth=1.5, linestyle='--',
                       label=model_name, color=colors[model_name], 
                       alpha=0.6, zorder=2)
        
        # Add Monegetti
        monegetti_pred = monegetti_model(age_range)
        ax.plot(age_range, monegetti_pred, linewidth=2, linestyle=':',
               label='Monegetti', color='#9b59b6', alpha=0.8, zorder=3)
        
        ax.set_xlabel('Larval Age (days)', fontsize=11, fontweight='bold')
        ax.set_ylabel('Proportion Competent', fontsize=11, fontweight='bold')
        ax.set_title(f'{family}', fontsize=12, fontweight='bold')
        ax.legend(loc='lower right', fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, min(80, df_agg['LarvalAge'].max() + 10))
        ax.set_ylim(0, 1.0)
    
    # Hide extra subplots
    for idx in range(n_families, len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved: {save_path}")
    plt.show()

def create_summary_table(results):
    """Create summary table."""
    summary_data = []
    
    for family, result in results.items():
        best_model = result['best_model']
        params = result['models'][best_model]['params']
        aic = result['models'][best_model]['aic']
        df_agg = result['data']
        
        if params is not None:
            if best_model == 'Logistic':
                L, k, x0 = params
                param_str = f"L={L:.3f}, k={k:.3f}, x0={x0:.1f}"
                tc50 = x0  # Age at 50% competency
            elif best_model == 'Gompertz':
                a, b, c = params
                param_str = f"a={a:.3f}, b={b:.2f}, c={c:.3f}"
                # Approximate TC50 for Gompertz
                tc50 = np.log(b) / c if c > 0 else np.nan
            elif best_model == 'Weibull':
                tc, a, b, v = params
                param_str = f"tc={tc:.1f}, a={a:.3f}, b={b:.3f}, v={v:.3f}"
                tc50 = tc + 2  # Rough approximation
            else:
                param_str = "N/A"
                tc50 = np.nan
        else:
            param_str = "Failed"
            tc50 = np.nan
        
        summary_data.append({
            'Family': family,
            'N_ages': len(df_agg),
            'Age_range': f"{df_agg['LarvalAge'].min():.0f}-{df_agg['LarvalAge'].max():.0f}",
            'Max_Prop': f"{df_agg['PropCompetent'].max():.3f}",
            'Best_Model': best_model,
            'TC50_approx': f"{tc50:.1f}" if not np.isnan(tc50) else "N/A",
            'Parameters': param_str,
            'AIC': f"{aic:.1f}"
        })
    
    return pd.DataFrame(summary_data)

# ============================================================================
# MAIN
# ============================================================================

def main():
    """Main analysis."""
    print("\n" + "=" * 70)
    print("HYBRID COMPETENCY MODEL")
    print("Continuous Proportions + Monegetti-Style Fitting")
    print("=" * 70)
    
    # Paths
    data_path = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity/R_codes/Randal/data/Settlement.csv'
    output_dir = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity/R_codes/'
    
    # Prepare data (continuous proportions)
    family_data = prepare_data_continuous(data_path, treatment='rubble')
    
    if len(family_data) == 0:
        print("\nNo families with sufficient data!")
        return
    
    # Fit models
    results = analyze_families(family_data)
    
    # Summary table
    print(f"\n{'=' * 70}")
    print("SUMMARY TABLE")
    print(f"{'=' * 70}\n")
    summary_df = create_summary_table(results)
    print(summary_df.to_string(index=False))
    summary_df.to_csv(output_dir + 'hybrid_model_summary.csv', index=False)
    
    # Visualization
    plot_family_results(results, 
                       save_path=output_dir + 'hybrid_competency_models.png')
    
    # Generate predictions
    print(f"\n{'=' * 70}")
    print("PREDICTIONS AT KEY AGES")
    print(f"{'=' * 70}\n")
    
    key_ages = np.array([3, 5, 7, 10, 15, 20, 30, 40, 50])
    all_predictions = {'Age (days)': key_ages}
    
    for family, result in results.items():
        best_model = result['best_model']
        params = result['models'][best_model]['params']
        
        if params is not None:
            if best_model == 'Logistic':
                pred = logistic_competency(key_ages, *params)
            elif best_model == 'Gompertz':
                pred = gompertz_competency(key_ages, *params)
            elif best_model == 'Weibull':
                pred = weibull_competency(key_ages, *params)
            
            all_predictions[family] = pred
    
    # Add Monegetti
    all_predictions['Monegetti'] = monegetti_model(key_ages)
    
    pred_df = pd.DataFrame(all_predictions)
    print(pred_df.to_string(index=False, float_format=lambda x: f'{x:.3f}' if x < 10 else f'{x:.0f}'))
    pred_df.to_csv(output_dir + 'hybrid_predictions.csv', index=False)
    
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print("\nKey Advantages of This Approach:")
    print("  1. Uses CONTINUOUS proportions (not binary)")
    print("  2. Allows natural fluctuations in settlement")
    print("  3. Avoids threshold artifacts")
    print("  4. Fits Monegetti-style models")
    print("  5. Weighted by sample size")
    print("\nOutput files:")
    print("  - hybrid_competency_models.png")
    print("  - hybrid_model_summary.csv")
    print("  - hybrid_predictions.csv")

if __name__ == "__main__":
    main()

