#!/usr/bin/env python3
"""
Generate LaTeX tables from analysis results.
"""

import pandas as pd
import numpy as np
import os
from typing import Dict, Any


def generate_model_comparison_table(results: Dict[str, Dict[str, Any]], 
                                    output_path: str) -> str:
    """
    Generate LaTeX table with model comparison (all models, all families).
    
    Parameters
    ----------
    results : Dict[str, Dict[str, Any]]
        Analysis results
    output_path : str
        Path to save LaTeX table
        
    Returns
    -------
    str
        LaTeX table code
    """
    rows = []
    
    for family, result in results.items():
        all_models = result['all_model_results']
        
        row = {'Family': family}
        
        # Add AICc and NLL for each model
        for model_name in ['Logistic', 'Gompertz', 'Weibull', 'Monegetti']:
            if model_name in all_models and all_models[model_name]['success']:
                row[f'{model_name}_AICc'] = f"{all_models[model_name]['aicc']:.1f}"
                row[f'{model_name}_NLL'] = f"{all_models[model_name]['nll']:.1f}"
            else:
                row[f'{model_name}_AICc'] = '--'
                row[f'{model_name}_NLL'] = '--'
        
        rows.append(row)
    
    df = pd.DataFrame(rows)
    
    # Generate LaTeX
    latex = "\\begin{tabular}{lccccccccc}\n"
    latex += "\\toprule\n"
    latex += "Family & \\multicolumn{2}{c}{Logistic} & \\multicolumn{2}{c}{Gompertz} & "
    latex += "\\multicolumn{2}{c}{Weibull} & \\multicolumn{2}{c}{Monegetti} \\\\\n"
    latex += "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}\n"
    latex += " & AICc & NLL & AICc & NLL & AICc & NLL & AICc & NLL \\\\\n"
    latex += "\\midrule\n"
    
    for _, row in df.iterrows():
        latex += f"{row['Family']} & "
        for model in ['Logistic', 'Gompertz', 'Weibull', 'Monegetti']:
            aicc = row[f'{model}_AICc']
            nll = row[f'{model}_NLL']
            latex += f"{aicc} & {nll} & "
        latex = latex.rstrip(' & ') + " \\\\\n"
    
    latex += "\\bottomrule\n"
    latex += "\\end{tabular}\n"
    
    # Save
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(latex)
    
    return latex


def generate_final_selection_table(results: Dict[str, Dict[str, Any]],
                                   output_path: str) -> str:
    """
    Generate LaTeX table with final model selection.
    
    Parameters
    ----------
    results : Dict[str, Dict[str, Any]]
        Analysis results
    output_path : str
        Path to save LaTeX table
        
    Returns
    -------
    str
        LaTeX table code
    """
    rows = []
    
    for family, result in results.items():
        best_model = result['best_model']
        best_result = result['best_model_result']
        all_models = result['all_model_results']
        
        # Calculate delta AICc
        best_aicc = best_result['aicc']
        deltas = {}
        for model_name in ['Logistic', 'Gompertz', 'Weibull', 'Monegetti']:
            if model_name in all_models and all_models[model_name]['success']:
                deltas[model_name] = all_models[model_name]['aicc'] - best_aicc
            else:
                deltas[model_name] = np.nan
        
        # Get parameters for best model
        params = best_result['params']
        param_names = best_result['param_names']
        
        row = {
            'Family': family,
            'Best_Model': best_model,
            'N_Replicates': result['n_replicates'],
            'N_Larvae': result['n_larvae'],
            'AICc': f"{best_aicc:.1f}",
            'TC50': f"{result['tc50']:.2f}" if result['tc50'] is not None else "N/A"
        }
        
        # Add parameter values
        for name, value in zip(param_names, params):
            row[f'Param_{name}'] = f"{value:.4f}"
        
        rows.append(row)
    
    df = pd.DataFrame(rows)
    
    # Generate LaTeX
    latex = "\\begin{tabular}{llrrrr}\n"
    latex += "\\toprule\n"
    latex += "Family & Best Model & N Replicates & N Larvae & AICc & TC50 (days) \\\\\n"
    latex += "\\midrule\n"
    
    for _, row in df.iterrows():
        latex += f"{row['Family']} & {row['Best_Model']} & "
        latex += f"{row['N_Replicates']} & {row['N_Larvae']} & "
        latex += f"{row['AICc']} & {row['TC50']} \\\\\n"
    
    latex += "\\bottomrule\n"
    latex += "\\end{tabular}\n"
    
    # Save
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(latex)
    
    return latex

