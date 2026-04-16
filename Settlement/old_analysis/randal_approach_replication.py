#!/usr/bin/env python3
"""
Replication of Randal's Larval Competency Approach
===================================================
This script replicates Randal's methodology from his GitHub repository:
https://github.com/open-AIMS/larval_settlement_competency

Key features of Randal's approach:
1. Add age=0 observations (NoSet=0, NoNotSet=10)
2. Threshold-based binary settlement scoring (0.1 to 0.9)
3. Cumulative settlement tracking (once exceeds threshold, stays 1)
4. Species-level analysis (not family)
5. Logistic regression: Settlement ~ Treatment * LarvalAge
6. Calculate LD50 (age at 50% competency)

Adaptation: Uses frequentist GLM instead of Bayesian brms
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from scipy.optimize import minimize_scalar
import warnings
warnings.filterwarnings('ignore')

sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 10)

# ============================================================================
# DATA PROCESSING (Following Randal's 10_processData.R)
# ============================================================================

def add_age_zero_observations(df):
    """
    Add age=0 observations for each species and treatment combination.
    At birth, larvae cannot settle (biological constraint).
    
    This replicates Randal's approach in 10_processData.R lines 22-44
    """
    # Get unique combinations of species, family, treatment, and spawn date
    unique_combos = df[['Species', 'Family', 'SpecificTreatment', 'SpawnDate']].drop_duplicates()
    
    age_zero_rows = []
    for _, row in unique_combos.iterrows():
        age_zero_rows.append({
            'Species': row['Species'],
            'Family': row['Family'],
            'SpecificTreatment': row['SpecificTreatment'],
            'SpawnDate': row['SpawnDate'],
            'ReadDate': row['SpawnDate'],
            'LarvalAge': 0,
            'NoSet': 0,
            'NoNotSet': 10,
            'NoAlive': 10,
            'Plate': f"{row['SpawnDate']}_0_AGE0"
        })
    
    df_age_zero = pd.DataFrame(age_zero_rows)
    df_combined = pd.concat([df, df_age_zero], ignore_index=True)
    
    print(f"Added {len(age_zero_rows)} age=0 observations")
    return df_combined

def calculate_cumulative_settlement(df, threshold):
    """
    Calculate cumulative binary settlement for a given threshold.
    
    This replicates Randal's approach in 20_fitModels.R lines 32-38:
    1. Score as 1 if proportion settled > threshold, else 0
    2. Sort by LarvalAge
    3. Calculate cumsum - once it exceeds threshold, stays competent (1)
    
    Parameters:
    -----------
    df : DataFrame
        Data for one species and treatment
    threshold : float
        Settlement proportion threshold (e.g., 0.3 means 30%)
    
    Returns:
    --------
    df : DataFrame
        With added 'Settle' and 'Settlement' columns
    """
    df = df.copy()
    
    # Calculate fraction settled
    df['FractionSettled'] = df['NoSet'] / (df['NoSet'] + df['NoNotSet'])
    
    # Binary: exceeds threshold?
    df['Settle'] = (df['FractionSettled'] > threshold).astype(int)
    
    # Sort by age and calculate cumulative
    df = df.sort_values('LarvalAge')
    df['cumsumSettle'] = df['Settle'].cumsum()
    
    # Once cumsum > 0, competent forever
    df['Settlement'] = (df['cumsumSettle'] > 0).astype(int)
    
    return df

def prepare_data_randal_style(data_path, threshold=0.3, treatment='rubble'):
    """
    Prepare data following Randal's methodology.
    
    Steps:
    1. Load data
    2. Filter to exclude age=0 (will add back properly)
    3. Add age=0 observations
    4. Filter to specified treatment (default: rubble)
    5. Calculate cumulative settlement
    """
    # Load data
    df = pd.read_csv(data_path)
    
    print("=" * 70)
    print("DATA PREPARATION (Randal's Method - Modified for Family)")
    print("=" * 70)
    print(f"Total raw observations: {len(df)}")
    
    # Rename column if needed
    if 'Treatment' in df.columns and 'SpecificTreatment' not in df.columns:
        df['SpecificTreatment'] = df['Treatment']
    
    # Remove any existing age=0 (will add back properly)
    df = df[df['LarvalAge'] != 0].copy()
    print(f"After removing age=0: {len(df)}")
    
    # Add age=0 observations
    df = add_age_zero_observations(df)
    print(f"After adding age=0: {len(df)}")
    
    # Filter to specified treatment only
    df_treatment = df[df['SpecificTreatment'] == treatment].copy()
    print(f"After filtering to {treatment}: {len(df_treatment)}")
    
    print(f"\nTreatment: {treatment}")
    print(f"Species: {df_treatment['Species'].nunique()}")
    print(f"Families: {df_treatment['Family'].nunique()}")
    
    # Calculate cumulative settlement for each FAMILY (not species)
    results = []
    for family, group in df_treatment.groupby(['Family']):
        group_processed = calculate_cumulative_settlement(group, threshold)
        results.append(group_processed)
    
    df_final = pd.concat(results, ignore_index=True)
    
    # Remove observations with no larvae
    df_final = df_final[(df_final['NoSet'] + df_final['NoNotSet']) > 0]
    
    print(f"\nFinal dataset size: {len(df_final)}")
    print(f"Threshold used: {threshold}")
    
    return df_final

# ============================================================================
# MODEL FITTING (Simple Logistic Regression)
# ============================================================================

def fit_logistic_model(df_family):
    """
    Fit logistic regression: Settlement ~ LarvalAge
    
    This is a simplified version of Randal's Bayesian model (20_fitModels.R line 96)
    Using sklearn's LogisticRegression instead of brms
    Now fits for a single family (no treatment effects since we're using only rubble)
    
    Returns:
    --------
    model_info : dict
        Dictionary with fitted model information
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    
    if len(df_family) < 5:  # Need minimum data
        return None
    
    X = df_family[['LarvalAge']].values
    y = df_family['Settlement'].values
    
    # Standardize age for better convergence
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Fit logistic regression
    model = LogisticRegression(max_iter=1000, random_state=42)
    
    try:
        model.fit(X_scaled, y)
        model_info = {
            'model': model,
            'scaler': scaler,
            'n_obs': len(df_family),
            'ages': df_family['LarvalAge'].values,
            'settlement': y
        }
        return model_info
    except:
        print(f"  Warning: Could not fit model")
        return None

def calculate_ld50(model_info):
    """
    Calculate LD50: age at which 50% probability of competency.
    
    This is Randal's key metric (see 40_summariseModels.R)
    """
    model = model_info['model']
    scaler = model_info['scaler']
    
    def prob_at_age(age):
        age_scaled = scaler.transform([[age]])
        prob = model.predict_proba(age_scaled)[0, 1]
        return abs(prob - 0.5)  # Distance from 0.5
    
    # Find age where probability = 0.5
    min_age = model_info['ages'].min()
    max_age = model_info['ages'].max()
    
    try:
        result = minimize_scalar(prob_at_age, bounds=(min_age, max_age), method='bounded')
        return result.x
    except:
        return np.nan

def predict_competency_curve(model_info, age_range):
    """
    Generate predicted competency curve for visualization.
    """
    model = model_info['model']
    scaler = model_info['scaler']
    
    ages = np.array(age_range).reshape(-1, 1)
    ages_scaled = scaler.transform(ages)
    probs = model.predict_proba(ages_scaled)[:, 1]
    
    return probs

# ============================================================================
# ANALYSIS FOR ALL SPECIES
# ============================================================================

def analyze_all_families(df, threshold=0.3, treatment='rubble'):
    """
    Analyze all families separately (modified from Randal's approach).
    """
    results = {}
    
    family_list = df['Family'].unique()
    print(f"\n{'=' * 70}")
    print(f"ANALYZING {len(family_list)} FAMILIES (Threshold = {threshold}, Treatment = {treatment})")
    print(f"{'=' * 70}\n")
    
    for family in sorted(family_list):
        if pd.isna(family):  # Skip NaN families
            continue
            
        df_family = df[df['Family'] == family].copy()
        
        print(f"\n{family}")
        print(f"  Observations: {len(df_family)}")
        print(f"  Age range: {df_family['LarvalAge'].min():.0f}-{df_family['LarvalAge'].max():.0f} days")
        print(f"  Species in family: {df_family['Species'].nunique()}")
        
        # Fit model
        model_info = fit_logistic_model(df_family)
        
        if model_info is None:
            print(f"  Warning: Model not fitted (insufficient data)")
            continue
        
        # Calculate LD50
        ld50 = calculate_ld50(model_info)
        if not np.isnan(ld50):
            print(f"  LD50: {ld50:.1f} days")
        else:
            print(f"  LD50: unable to calculate")
        
        results[family] = {
            'data': df_family,
            'model': model_info,
            'ld50': ld50,
            'n_species': df_family['Species'].nunique()
        }
    
    return results

# ============================================================================
# VISUALIZATION
# ============================================================================

def plot_family_competency(results, threshold, treatment, save_dir='./'):
    """
    Create visualization for each family (modified from Randal's approach).
    """
    color = '#00BA38'  # Green for rubble
    
    for family, result in results.items():
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        fig.suptitle(f'{family} - Threshold = {threshold} - Treatment = {treatment}', 
                     fontsize=14, fontweight='bold')
        
        df_family = result['data']
        model_info = result['model']
        
        # Plot 1: Raw settlement proportions
        ax1 = axes[0]
        
        # Plot observed data
        ax1.scatter(df_family['LarvalAge'], df_family['FractionSettled'],
                   alpha=0.5, s=60, label=f'{treatment} (n={len(df_family)})',
                   color=color, edgecolors='white', linewidth=0.5)
        
        # Plot model prediction
        if model_info is not None:
            age_range = np.linspace(df_family['LarvalAge'].min(),
                                   df_family['LarvalAge'].max(), 100)
            pred = predict_competency_curve(model_info, age_range)
            ax1.plot(age_range, pred, linewidth=2.5, color=color, alpha=0.7)
        
        ax1.set_xlabel('Larval Age (days)', fontsize=11, fontweight='bold')
        ax1.set_ylabel('Settlement Proportion', fontsize=11, fontweight='bold')
        ax1.set_title('Raw Settlement Data', fontsize=11, fontweight='bold')
        ax1.legend(loc='lower right')
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(-0.05, 1.05)
        
        # Plot 2: Binary cumulative settlement
        ax2 = axes[1]
        
        # Plot binary settlement
        colors_map = df_family['Settlement'].map({0: 'lightgray', 1: color})
        ax2.scatter(df_family['LarvalAge'], df_family['Settlement'],
                   alpha=0.6, s=60, c=colors_map,
                   edgecolors='white', linewidth=0.5)
        
        # Plot model prediction
        if model_info is not None:
            age_range = np.linspace(df_family['LarvalAge'].min(),
                                   df_family['LarvalAge'].max(), 100)
            pred = predict_competency_curve(model_info, age_range)
            ax2.plot(age_range, pred, linewidth=2.5,
                    color=color, label=treatment, alpha=0.8)
            
            # Add LD50 marker
            if not np.isnan(result['ld50']):
                ld50 = result['ld50']
                ax2.axvline(ld50, color=color,
                           linestyle='--', alpha=0.5, linewidth=1.5)
                ax2.text(ld50, 0.95, f'LD50={ld50:.1f}d',
                        rotation=90, va='top', ha='right',
                        fontsize=9, color=color)
        
        ax2.set_xlabel('Larval Age (days)', fontsize=11, fontweight='bold')
        ax2.set_ylabel('Cumulative Settlement (0/1)', fontsize=11, fontweight='bold')
        ax2.set_title(f'Binary Settlement (threshold={threshold})', fontsize=11, fontweight='bold')
        ax2.legend(loc='lower right')
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim(-0.05, 1.05)
        
        plt.tight_layout()
        
        # Save
        filename = f"{save_dir}randal_family_{family.replace(' ', '_')}_thresh{threshold}_{treatment}.png"
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"  Saved: {filename}")

def create_summary_comparison(results_by_threshold, treatment, save_path):
    """
    Create summary plot comparing LD50 across thresholds for families.
    """
    summary_data = []
    
    for threshold, results in results_by_threshold.items():
        for family, result in results.items():
            ld50 = result['ld50']
            if not np.isnan(ld50):
                summary_data.append({
                    'Family': family,
                    'Threshold': threshold,
                    'LD50': ld50,
                    'N_Species': result['n_species']
                })
    
    df_summary = pd.DataFrame(summary_data)
    
    if len(df_summary) == 0:
        print("No LD50 values to plot")
        return df_summary
    
    # Create visualization
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(f'LD50 Summary by Family: Age at 50% Competency\n(Treatment = {treatment}, Randal\'s Approach)',
                 fontsize=14, fontweight='bold')
    
    color = '#00BA38'  # Green for rubble
    
    # Plot 1: LD50 by family and threshold
    ax1 = axes[0]
    family_order = df_summary.groupby('Family')['LD50'].mean().sort_values().index
    
    thresholds = sorted(df_summary['Threshold'].unique())
    width = 0.25
    
    for i, thresh in enumerate(thresholds):
        df_thresh = df_summary[df_summary['Threshold'] == thresh]
        df_thresh_family = df_thresh.groupby('Family')['LD50'].mean().reindex(family_order)
        
        x_pos = np.arange(len(family_order))
        offset = (i - len(thresholds)/2 + 0.5) * width
        
        ax1.bar(x_pos + offset, df_thresh_family.values, width=width,
               label=f'Threshold={thresh}', alpha=0.7)
    
    ax1.set_xticks(np.arange(len(family_order)))
    ax1.set_xticklabels(family_order, rotation=45, ha='right')
    ax1.set_ylabel('LD50 (days)', fontsize=11, fontweight='bold')
    ax1.set_xlabel('Family', fontsize=11, fontweight='bold')
    ax1.set_title('LD50 by Family and Threshold', fontsize=11, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Plot 2: LD50 comparison (threshold = 0.3)
    ax2 = axes[1]
    df_main = df_summary[df_summary['Threshold'] == 0.3].copy()
    
    if len(df_main) > 0:
        df_main = df_main.sort_values('LD50')
        
        bars = ax2.barh(range(len(df_main)), df_main['LD50'].values,
                       color=color, alpha=0.7, edgecolor='black', linewidth=1)
        
        ax2.set_yticks(range(len(df_main)))
        ax2.set_yticklabels(df_main['Family'].values)
        ax2.set_xlabel('LD50 (days)', fontsize=11, fontweight='bold')
        ax2.set_title('LD50 by Family (Threshold = 0.3)', fontsize=11, fontweight='bold')
        ax2.grid(True, alpha=0.3, axis='x')
        
        # Add labels with number of species
        for i, (idx, row) in enumerate(df_main.iterrows()):
            ax2.text(row['LD50'] + 0.5, i, f"n={row['N_Species']}", 
                    va='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"\nSummary plot saved: {save_path}")
    plt.close()
    
    return df_summary

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

def main():
    """
    Main analysis replicating Randal's approach - FAMILY LEVEL, RUBBLE ONLY.
    """
    print("\n" + "=" * 70)
    print("RANDAL'S APPROACH REPLICATION")
    print("Family-level cumulative binary settlement analysis")
    print("Treatment: RUBBLE ONLY")
    print("=" * 70)
    
    # Paths
    data_path = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity/R_codes/Randal/data/Settlement.csv'
    output_dir = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity/R_codes/'
    
    # Treatment to analyze
    treatment = 'rubble'
    
    # Test multiple thresholds (like Randal does in 20_fitModels.R line 10)
    thresholds = [0.3, 0.5, 0.7]  # Randal tests 0.1 to 0.9, we'll do a subset
    
    # Option to create individual plots
    create_individual_plots = True  # Set to False to skip
    
    results_by_threshold = {}
    
    for threshold in thresholds:
        print(f"\n{'=' * 70}")
        print(f"THRESHOLD = {threshold}")
        print(f"{'=' * 70}")
        
        # Prepare data
        df = prepare_data_randal_style(data_path, threshold=threshold, treatment=treatment)
        
        # Analyze all families
        results = analyze_all_families(df, threshold=threshold, treatment=treatment)
        results_by_threshold[threshold] = results
        
        # Create plots for each family
        if create_individual_plots and len(results) > 0:
            print(f"\nGenerating individual family plots...")
            plot_family_competency(results, threshold, treatment, save_dir=output_dir)
        elif not create_individual_plots:
            print(f"\nSkipping individual family plots")
    
    # Create summary comparison across thresholds
    print(f"\n{'=' * 70}")
    print("CREATING SUMMARY COMPARISON")
    print(f"{'=' * 70}")
    
    df_summary = create_summary_comparison(
        results_by_threshold,
        treatment,
        save_path=output_dir + f'randal_ld50_summary_{treatment}.png'
    )
    
    # Save summary table
    if len(df_summary) > 0:
        summary_file = output_dir + f'randal_ld50_summary_{treatment}.csv'
        df_summary.to_csv(summary_file, index=False)
        print(f"Summary table saved: {summary_file}")
        
        print("\n" + "=" * 70)
        print("LD50 SUMMARY TABLE")
        print("=" * 70)
        
        # Pivot table for easy viewing
        pivot = df_summary.pivot_table(
            values='LD50',
            index='Family',
            columns='Threshold',
            aggfunc='mean'
        )
        print(pivot.to_string())
    
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print("\nModifications from original Randal approach:")
    print("  ✓ Data processing: Same (age=0, cumulative binary)")
    print("  ✓ Threshold testing: Subset (0.3, 0.5, 0.7 instead of 0.1-0.9)")
    print("  ✗ Taxonomic level: FAMILY (not species)")
    print("  ✗ Treatment: RUBBLE ONLY (not CCA/rubble/disc)")
    print("  ✗ Statistical model: Frequentist GLM instead of Bayesian brms")
    print("  ✗ Random effects: Not included (simplified)")
    print("  ✓ LD50 calculation: Same concept")
    print("  ✓ Visualizations: Similar style")
    
    print("\nOutput files:")
    print(f"  - randal_family_<family>_thresh<X>_{treatment}.png (per family)")
    print(f"  - randal_ld50_summary_{treatment}.png (comparison)")
    print(f"  - randal_ld50_summary_{treatment}.csv (LD50 values)")

if __name__ == "__main__":
    main()

