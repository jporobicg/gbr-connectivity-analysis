#!/usr/bin/env python3
"""
Monegetti Piecewise Weibull-Exponential Competency Model - Main Orchestrator
============================================================================

This is the main entry point for the larval competency modeling analysis.
It orchestrates the full workflow: data loading → model fitting → diagnostics → outputs.

Based on:
Moneghetti et al. (2019) "High-frequency sampling and piecewise models 
reshape dispersal kernels of a common reef coral"

Model Description:
- Precompetency period (t < tc): No settlement possible
- Early competency (tc < t < Tcp): Weibull loss of competency
- Late competency (t > Tcp): Exponential loss of competency

All helper functions are in utility_tools.py
All plotting functions are in figures_analysis.py
"""

import pandas as pd
import numpy as np
import warnings
import os
import sys
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any

# Import utility functions
from utility_tools import (
    load_settlement_data,
    logistic_competency,
    gompertz_competency,
    weibull_competency,
    monegetti_competency,
    fit_logistic_binomial,
    fit_gompertz_binomial,
    fit_weibull_binomial,
    fit_monegetti_binomial,
    neg_log_likelihood_monegetti,
    calculate_aicc,
    estimate_overdispersion,
    calculate_residuals,
    calculate_tc50_monegetti,
    sanitize_filename,
    ensure_directory,
    get_monegetti_bounds
)

warnings.filterwarnings('ignore')

# Import figures module (will be created)
try:
    import figures_analysis
except ImportError:
    figures_analysis = None
    print("WARNING: figures_analysis module not found. Figures will not be generated.")


# ============================================================================
# MODEL FORMULATIONS (All models defined here)
# ============================================================================

MODEL_DEFINITIONS = {
    'Logistic': {
        'function': logistic_competency,
        'n_params': 3,
        'param_names': ['L', 'k', 'x0'],
        'description': 'Logistic (sigmoidal) model',
        'origin': 'Standard sigmoid function'
    },
    'Gompertz': {
        'function': gompertz_competency,
        'n_params': 3,
        'param_names': ['a', 'b', 'c'],
        'description': 'Gompertz growth model',
        'origin': 'Gompertz (1825)'
    },
    'Weibull': {
        'function': weibull_competency,
        'n_params': 4,
        'param_names': ['tc', 'a', 'b', 'v'],
        'description': 'Weibull-based competency model',
        'origin': 'Weibull distribution adaptation'
    },
    'Monegetti': {
        'function': monegetti_competency,
        'n_params': 6,
        'param_names': ['a', 'b1', 'v1', 'b2', 'tc', 'Tcp'],
        'description': 'Piecewise Weibull-Exponential model',
        'origin': 'Moneghetti et al. (2019)'
    }
}


# ============================================================================
# MODEL FITTING WORKFLOW
# ============================================================================

def fit_all_models(df: pd.DataFrame) -> Dict[str, Dict[str, Any]]:
    """
    Fit all candidate models to a dataset.
    
    Parameters
    ----------
    df : pd.DataFrame
        Data with columns: LarvalAge, NoSet, NoAlive
        
    Returns
    -------
    Dict[str, Dict[str, Any]]
        Dictionary with results for each model
    """
    results = {}
    n_obs = len(df)
    
    # Fit each model
    for model_name, model_def in MODEL_DEFINITIONS.items():
        if model_name == 'Monegetti':
            params, nll = fit_monegetti_binomial(df, model_def['function'])
        else:
            fit_func = {
                'Logistic': fit_logistic_binomial,
                'Gompertz': fit_gompertz_binomial,
                'Weibull': fit_weibull_binomial
            }[model_name]
            params, nll = fit_func(df)
        
        if params is not None:
            aicc = calculate_aicc(nll, model_def['n_params'], n_obs)
            
            # Calculate diagnostics
            phi, pearson_resid = estimate_overdispersion(
                df, params, model_def['function']
            )
            residuals = calculate_residuals(df, params, model_def['function'])
            
            results[model_name] = {
                'params': params,
                'nll': nll,
                'aicc': aicc,
                'phi': phi,
                'pearson_residuals': pearson_resid,
                'residuals': residuals,
                'n_params': model_def['n_params'],
                'param_names': model_def['param_names'],
                'success': True
            }
        else:
            results[model_name] = {
                'params': None,
                'nll': np.inf,
                'aicc': np.inf,
                'success': False
            }
    
    return results


def select_best_model(model_results: Dict[str, Dict[str, Any]]) -> Tuple[str, Dict[str, Any]]:
    """
    Select best model based on AICc.
    
    Parameters
    ----------
    model_results : Dict[str, Dict[str, Any]]
        Results from fit_all_models
        
    Returns
    -------
    Tuple[str, Dict[str, Any]]
        (best_model_name, best_model_results)
    """
    # Filter to successful models
    successful = {
        name: res for name, res in model_results.items()
        if res.get('success', False)
    }
    
    if not successful:
        return None, None
    
    # Select by AICc
    best_model = min(successful, key=lambda x: successful[x]['aicc'])
    return best_model, successful[best_model]


# ============================================================================
# FAMILY ANALYSIS
# ============================================================================

def analyze_family(
    family: str,
    df_family: pd.DataFrame,
    min_replicates: int = 30,
    min_larvae: int = 200
) -> Optional[Dict[str, Any]]:
    """
    Analyze a single family: fit all models, select best, calculate diagnostics.
    
    Parameters
    ----------
    family : str
        Family name
    df_family : pd.DataFrame
        Family data
    min_replicates : int, optional
        Minimum number of replicates required, by default 30
    min_larvae : int, optional
        Minimum number of larvae required, by default 200
        
    Returns
    -------
    Optional[Dict[str, Any]]
        Analysis results or None if insufficient data
    """
    # Check data requirements
    n_replicates = len(df_family)
    n_larvae = df_family['NoAlive'].sum()
    
    if n_replicates < min_replicates or n_larvae < min_larvae:
        return None
    
    print(f"\n{'='*70}")
    print(f"FAMILY: {family}")
    print(f"{'='*70}")
    print(f"  Replicates: {n_replicates}")
    print(f"  Age range: {df_family['LarvalAge'].min():.1f}-{df_family['LarvalAge'].max():.1f} days")
    print(f"  Total larvae: {n_larvae:.0f}")
    print(f"  Total settled: {df_family['NoSet'].sum():.0f}")
    
    # Fit all models
    print(f"\n  Fitting all candidate models...")
    model_results = fit_all_models(df_family)
    
    # Select best model
    best_model_name, best_model_result = select_best_model(model_results)
    
    if best_model_name is None:
        print(f"  WARNING: All model fits failed")
        return None
    
    print(f"\n  Model Comparison (AICc):")
    for model_name in ['Logistic', 'Gompertz', 'Weibull', 'Monegetti']:
        if model_name in model_results and model_results[model_name]['success']:
            aicc = model_results[model_name]['aicc']
            delta = aicc - best_model_result['aicc']
            marker = "✓" if model_name == best_model_name else " "
            print(f"    {marker} {model_name:12s}: {aicc:10.1f}  (Δ = {delta:6.1f})")
    
    print(f"\n  → Best model: {best_model_name}")
    
    # Calculate TC50 for Monegetti model
    tc50 = None
    max_competency = None
    if best_model_name == 'Monegetti':
        params = best_model_result['params']
        max_age = df_family['LarvalAge'].max()
        tc50, max_competency = calculate_tc50_monegetti(
            params, max_age, monegetti_competency
        )
        if tc50 is not None:
            print(f"  TC50: {tc50:.2f} days")
            print(f"  Max competency: {max_competency:.3f}")
    
    # Compile results
    result = {
        'family': family,
        'best_model': best_model_name,
        'best_model_result': best_model_result,
        'all_model_results': model_results,
        'n_replicates': n_replicates,
        'n_larvae': n_larvae,
        'age_range': (df_family['LarvalAge'].min(), df_family['LarvalAge'].max()),
        'tc50': tc50,
        'max_competency': max_competency
    }
    
    return result


# ============================================================================
# MAIN WORKFLOW ORCHESTRATOR
# ============================================================================

def run_full_analysis(
    data_path: str,
    output_dir: str = 'results',
    figures_dir: str = 'figures_analysis',
    test_mode: bool = False,
    treatments: Optional[List[str]] = None
) -> Dict[str, Any]:
    """
    Run the complete analysis pipeline.
    
    Parameters
    ----------
    data_path : str
        Path to settlement data CSV
    output_dir : str, optional
        Directory for results, by default 'results'
    figures_dir : str, optional
        Directory for figures, by default 'figures_analysis'
    test_mode : bool, optional
        If True, only analyze first family, by default False
        
    Returns
    -------
    Dict[str, Any]
        Complete analysis results
    """
    print("="*70)
    print("MONEGETTI PIECEWISE WEIBULL-EXPONENTIAL MODEL ANALYSIS")
    print("="*70)
    if test_mode:
        print("TEST MODE: Analyzing first family only")
    print("="*70)
    
    # Create output directories
    ensure_directory(output_dir)
    ensure_directory(figures_dir)
    
    # Load data
    print("\nLoading settlement data...")
    if treatments:
        print(f"  Filtering to treatments: {treatments}")
    df = load_settlement_data(data_path, treatments=treatments)
    print(f"  Total observations: {len(df)}")
    print(f"  Families: {sorted(df['Family'].unique())}")
    if 'Treatment' in df.columns:
        print(f"  Treatments included: {sorted(df['Treatment'].unique())}")
    
    # Analyze each family
    families_to_analyze = sorted(df['Family'].unique())
    if test_mode:
        families_to_analyze = families_to_analyze[:1]
    
    all_results = {}
    all_data = {}  # Store data for each family for figures
    summary_data = []
    
    print(f"\nAnalyzing {len(families_to_analyze)} families...")
    
    for idx, family in enumerate(families_to_analyze, 1):
        print(f"\n[{idx}/{len(families_to_analyze)}] Processing: {family}")
        
        df_family = df[df['Family'] == family].copy()
        result = analyze_family(family, df_family)
        
        if result is not None:
            all_results[family] = result
            all_data[family] = df_family  # Store for figures
            
            # Add to summary
            best = result['best_model_result']
            summary_data.append({
                'Family': family,
                'Best_Model': result['best_model'],
                'N_Replicates': result['n_replicates'],
                'N_Larvae': result['n_larvae'],
                'AICc': f"{best['aicc']:.1f}",
                'NLL': f"{best['nll']:.1f}",
                'Phi': f"{best['phi']:.2f}",
                'TC50': f"{result['tc50']:.2f}" if result['tc50'] is not None else "N/A"
            })
            
            # Save individual family result
            save_family_result(family, result, output_dir)
            
            # Generate figures if module available
            if figures_analysis is not None:
                try:
                    figures_analysis.create_family_figure(
                        family, df_family, result, figures_dir
                    )
                except Exception as e:
                    print(f"  WARNING: Figure generation failed: {e}")
    
    # Create summary table
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_path = os.path.join(output_dir, 'model_comparison_summary.csv')
        summary_df.to_csv(summary_path, index=False)
        print(f"\n{'='*70}")
        print("SUMMARY")
        print(f"{'='*70}")
        print(summary_df.to_string(index=False))
        print(f"\n✓ Summary saved: {summary_path}")
    
    # Generate comparison figures if module available
    if figures_analysis is not None:
        try:
            figures_analysis.create_theoretical_comparison_figure(figures_dir)
            figures_analysis.create_all_families_figure(
                all_results, all_data, figures_dir
            )
        except Exception as e:
            print(f"WARNING: Comparison figure generation failed: {e}")
    
    # Generate LaTeX tables if module available
    try:
        from generate_latex_tables import (
            generate_model_comparison_table,
            generate_final_selection_table
        )
        
        tables_dir = os.path.join(output_dir, 'tables')
        ensure_directory(tables_dir)
        
        generate_model_comparison_table(
            all_results,
            os.path.join(tables_dir, 'model_comparison_table.tex')
        )
        generate_final_selection_table(
            all_results,
            os.path.join(tables_dir, 'final_model_selection.tex')
        )
        print(f"\n✓ LaTeX tables saved to: {tables_dir}/")
    except ImportError:
        print("  Note: LaTeX table generation skipped (module not available)")
    
    return {
        'results': all_results,
        'summary': summary_df if summary_data else None,
        'n_families': len(all_results)
    }


def save_family_result(
    family: str,
    result: Dict[str, Any],
    output_dir: str
) -> str:
    """
    Save individual family result to CSV.
    
    Parameters
    ----------
    family : str
        Family name
    result : Dict[str, Any]
        Analysis result
    output_dir : str
        Output directory
        
    Returns
    -------
    str
        Path to saved file
    """
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    sanitized = sanitize_filename(family)
    filename = os.path.join(output_dir, f'{sanitized}_{timestamp}.csv')
    
    best = result['best_model_result']
    all_models = result['all_model_results']
    
    # Create comprehensive result row
    row = {
        'Family': family,
        'Timestamp': timestamp,
        'Best_Model': result['best_model'],
        'N_Replicates': result['n_replicates'],
        'N_Larvae': result['n_larvae'],
        'Best_AICc': best['aicc'],
        'Best_NLL': best['nll'],
        'Best_Phi': best['phi'],
        'TC50': result['tc50'] if result['tc50'] is not None else np.nan
    }
    
    # Add AICc for all models
    for model_name in ['Logistic', 'Gompertz', 'Weibull', 'Monegetti']:
        if model_name in all_models and all_models[model_name]['success']:
            row[f'{model_name}_AICc'] = all_models[model_name]['aicc']
            row[f'{model_name}_NLL'] = all_models[model_name]['nll']
        else:
            row[f'{model_name}_AICc'] = np.nan
            row[f'{model_name}_NLL'] = np.nan
    
    # Add parameters for best model
    param_names = best['param_names']
    params = best['params']
    for name, value in zip(param_names, params):
        row[f'Best_{name}'] = value
    
    df_result = pd.DataFrame([row])
    df_result.to_csv(filename, index=False)
    
    return filename


# ============================================================================
# MAIN ENTRY POINT
# ============================================================================

def main():
    """Main entry point for the analysis."""
    # Default data path
    default_data_path = (
        '/home/por07g/Documents/Projects/GBR_modeling/'
        'Connectivity_KERNELS/R_codes/Randal/data/Settlement.csv'
    )
    
    # Check command line arguments
    test_mode = '--test' in sys.argv
    
    # Filter to ecologically meaningful treatments: rubble and CCA
    treatments = ['rubble', 'CCA']
    
    # Run analysis
    results = run_full_analysis(
        data_path=default_data_path,
        output_dir='results',
        figures_dir='figures_analysis',
        test_mode=test_mode,
        treatments=treatments
    )
    
    print(f"\n{'='*70}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*70}")
    print(f"  Families analyzed: {results['n_families']}")
    print(f"  Results saved to: results/")
    print(f"  Figures saved to: figures_analysis/")
    print(f"{'='*70}")
    
    return results


if __name__ == "__main__":
    results = main()

