"""
Analyze settlement data to create family-specific competence functions.

This script processes empirical settlement data to derive competence curves
for different coral families.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

# Read the settlement data
print("Loading settlement data...")
df = pd.read_csv('data/Settlement.csv')

print(f"\nDataset shape: {df.shape}")
print(f"\nColumns: {df.columns.tolist()}")
print(f"\nFamilies in dataset: {df['Family'].unique()}")
print(f"\nSpecies in dataset: {df['Species'].unique()}")

# Filter for control treatments only (to get baseline settlement patterns)
df_control = df[df['Treatment'] == 'control'].copy()

print(f"\nControl treatment data: {len(df_control)} rows")

# Calculate settlement proportion
df_control['PropSettled'] = df_control['NoSet'] / df_control['NoAlive']

# Group by Family and LarvalAge to get average settlement
family_age_summary = df_control.groupby(['Family', 'LarvalAge']).agg({
    'NoSet': 'sum',
    'NoAlive': 'sum',
    'PropSettled': 'mean'
}).reset_index()

family_age_summary['CumulativeProp'] = family_age_summary.groupby('Family')['PropSettled'].cumsum()

print("\n" + "="*80)
print("FAMILY-LEVEL SETTLEMENT SUMMARY")
print("="*80)

for family in family_age_summary['Family'].unique():
    family_data = family_age_summary[family_age_summary['Family'] == family]
    
    print(f"\n{family}:")
    print(f"  Age range: {family_data['LarvalAge'].min()} - {family_data['LarvalAge'].max()} days")
    print(f"  Total larvae: {family_data['NoAlive'].sum()}")
    print(f"  Total settled: {family_data['NoSet'].sum()}")
    print(f"  Overall settlement rate: {family_data['NoSet'].sum() / family_data['NoAlive'].sum() * 100:.1f}%")
    
    # Find first settlement
    first_settlement = family_data[family_data['NoSet'] > 0]['LarvalAge'].min()
    print(f"  First settlement: day {first_settlement}")
    
    # Find peak settlement age (highest daily settlement)
    peak_age = family_data.loc[family_data['PropSettled'].idxmax(), 'LarvalAge']
    peak_prop = family_data['PropSettled'].max()
    print(f"  Peak settlement: day {peak_age} ({peak_prop*100:.1f}% daily rate)")

# Create detailed plot
families = family_age_summary['Family'].unique()
n_families = len(families)

fig, axes = plt.subplots(n_families, 2, figsize=(14, 4*n_families))
if n_families == 1:
    axes = axes.reshape(1, -1)

for idx, family in enumerate(families):
    family_data = family_age_summary[family_age_summary['Family'] == family].sort_values('LarvalAge')
    
    ages = family_data['LarvalAge'].values
    daily_prop = family_data['PropSettled'].values
    
    # Calculate cumulative settlement
    cumulative_settled = family_data['NoSet'].cumsum().values
    cumulative_total = family_data['NoAlive'].iloc[0] * len(family_data)  # Approximate
    cumulative_prop = cumulative_settled / family_data['NoAlive'].sum()
    
    # Plot 1: Daily settlement proportion
    ax1 = axes[idx, 0]
    ax1.bar(ages, daily_prop, alpha=0.7, color=f'C{idx}', edgecolor='black', linewidth=1)
    ax1.set_xlabel('Larval Age (days)', fontsize=11)
    ax1.set_ylabel('Proportion Settled (daily)', fontsize=11)
    ax1.set_title(f'{family} - Daily Settlement', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Cumulative settlement
    ax2 = axes[idx, 1]
    ax2.plot(ages, cumulative_prop, 'o-', color=f'C{idx}', linewidth=2, markersize=6)
    ax2.set_xlabel('Larval Age (days)', fontsize=11)
    ax2.set_ylabel('Cumulative Proportion Settled', fontsize=11)
    ax2.set_title(f'{family} - Cumulative Settlement', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, max(cumulative_prop)*1.1)

plt.tight_layout()
plt.savefig('settlement_data_analysis.png', dpi=300, bbox_inches='tight')
print("\n✓ Plot saved as 'settlement_data_analysis.png'")
plt.close()

# Now create a comparison of all families
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

for idx, family in enumerate(families):
    family_data = family_age_summary[family_age_summary['Family'] == family].sort_values('LarvalAge')
    
    ages = family_data['LarvalAge'].values
    daily_prop = family_data['PropSettled'].values
    
    # Calculate cumulative
    cumulative_settled = family_data['NoSet'].cumsum().values
    cumulative_prop = cumulative_settled / family_data['NoAlive'].sum()
    
    # Plot daily settlement
    ax1.plot(ages, daily_prop, 'o-', label=family, linewidth=2, markersize=5)
    
    # Plot cumulative settlement
    ax2.plot(ages, cumulative_prop, 'o-', label=family, linewidth=2, markersize=5)

ax1.set_xlabel('Larval Age (days)', fontsize=12)
ax1.set_ylabel('Daily Settlement Proportion', fontsize=12)
ax1.set_title('Daily Settlement by Family', fontsize=14, fontweight='bold')
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

ax2.set_xlabel('Larval Age (days)', fontsize=12)
ax2.set_ylabel('Cumulative Settlement Proportion', fontsize=12)
ax2.set_title('Cumulative Settlement by Family', fontsize=14, fontweight='bold')
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('family_settlement_comparison.png', dpi=300, bbox_inches='tight')
print("✓ Plot saved as 'family_settlement_comparison.png'")
plt.close()

# Export summary statistics for each family
print("\n" + "="*80)
print("EXPORTING FAMILY PARAMETERS")
print("="*80)

family_params = {}

for family in families:
    family_data = family_age_summary[family_age_summary['Family'] == family].sort_values('LarvalAge')
    
    ages = family_data['LarvalAge'].values
    daily_prop = family_data['PropSettled'].values
    
    # Find key parameters
    first_settlement = family_data[family_data['NoSet'] > 0]['LarvalAge'].min()
    peak_age = family_data.loc[family_data['PropSettled'].idxmax(), 'LarvalAge']
    peak_prop = family_data['PropSettled'].max()
    
    # Find age at 50% cumulative settlement
    cumulative_settled = family_data['NoSet'].cumsum().values
    cumulative_prop = cumulative_settled / family_data['NoAlive'].sum()
    
    if max(cumulative_prop) >= 0.5:
        age_50_idx = np.argmax(cumulative_prop >= 0.5)
        age_50 = ages[age_50_idx]
    else:
        age_50 = ages[-1]  # Use last age if never reaches 50%
    
    family_params[family] = {
        't_precomp': float(first_settlement - 0.5),  # Slightly before first settlement
        't_50': float(age_50),
        't_peak': float(peak_age),
        'peak_prop': float(peak_prop),
        'age_range': (float(ages.min()), float(ages.max()))
    }
    
    print(f"\n{family}:")
    print(f"  t_precomp (precompetency end): {family_params[family]['t_precomp']:.1f} days")
    print(f"  t_50 (50% settled): {family_params[family]['t_50']:.1f} days")
    print(f"  t_peak (peak settlement): {family_params[family]['t_peak']:.1f} days")
    print(f"  Peak daily proportion: {family_params[family]['peak_prop']:.3f}")

# Save parameters to file
import json
with open('family_settlement_parameters.json', 'w') as f:
    json.dump(family_params, f, indent=2)
print("\n✓ Parameters saved to 'family_settlement_parameters.json'")

print("\n" + "="*80)
print("ANALYSIS COMPLETE!")
print("="*80)


