"""
Create detailed plots of settlement patterns by family and age.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read the settlement data
print("Loading settlement data...")
df = pd.read_csv('data/Settlement.csv')

# Filter for control treatments only
df_control = df[df['Treatment'] == 'control'].copy()

# Calculate settlement proportion
df_control['PropSettled'] = df_control['NoSet'] / df_control['NoAlive']

# Get unique families (excluding NaN)
families = df_control['Family'].dropna().unique()
families = sorted(families)

print(f"Families to plot: {families}")

# Create a large figure with subplots for each family
n_families = len(families)
fig, axes = plt.subplots(n_families, 1, figsize=(14, 4*n_families))

# If only one family, make axes iterable
if n_families == 1:
    axes = [axes]

for idx, family in enumerate(families):
    ax = axes[idx]
    
    # Get data for this family
    family_data = df_control[df_control['Family'] == family].copy()
    
    # Group by age and calculate mean settlement proportion and total counts
    age_summary = family_data.groupby('LarvalAge').agg({
        'NoSet': 'sum',
        'NoAlive': 'sum',
        'PropSettled': 'mean'
    }).reset_index()
    
    # Calculate actual proportion from sums
    age_summary['ActualProp'] = age_summary['NoSet'] / age_summary['NoAlive']
    
    # Sort by age
    age_summary = age_summary.sort_values('LarvalAge')
    
    ages = age_summary['LarvalAge'].values
    prop_settled = age_summary['ActualProp'].values
    total_larvae = age_summary['NoAlive'].values
    num_settled = age_summary['NoSet'].values
    
    # Create bar plot with settlement proportion
    bars = ax.bar(ages, prop_settled, width=0.8, alpha=0.7, 
                   color=f'C{idx}', edgecolor='black', linewidth=1.5)
    
    # Add text labels on bars showing number settled / number alive
    for i, (age, prop, n_set, n_alive) in enumerate(zip(ages, prop_settled, num_settled, total_larvae)):
        if prop > 0.001:  # Only label bars with some settlement
            ax.text(age, prop, f'{int(n_set)}/{int(n_alive)}', 
                   ha='center', va='bottom', fontsize=8, fontweight='bold')
    
    # Styling
    ax.set_xlabel('Larval Age (days)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Proportion Settled', fontsize=12, fontweight='bold')
    ax.set_title(f'{family} - Settlement by Larval Age (Control Treatment)', 
                fontsize=14, fontweight='bold', pad=10)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_xlim(ages.min()-1, ages.max()+1)
    
    # Add summary statistics box
    total_settled = age_summary['NoSet'].sum()
    total_alive = age_summary['NoAlive'].sum()
    overall_rate = total_settled / total_alive * 100 if total_alive > 0 else 0
    first_settlement = age_summary[age_summary['NoSet'] > 0]['LarvalAge'].min()
    
    stats_text = (
        f'Total: {int(total_settled)}/{int(total_alive)} settled ({overall_rate:.1f}%)\n'
        f'First settlement: Day {first_settlement:.0f}\n'
        f'Age range: {ages.min():.0f}-{ages.max():.0f} days'
    )
    
    ax.text(0.98, 0.97, stats_text, transform=ax.transAxes, 
           fontsize=10, verticalalignment='top', horizontalalignment='right',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8, edgecolor='black', linewidth=1.5))

plt.tight_layout()
plt.savefig('settlement_by_family_detailed.png', dpi=300, bbox_inches='tight')
print("✓ Plot saved as 'settlement_by_family_detailed.png'")
plt.close()

# Create a second figure with just the main families (those with settlement)
main_families = []
for family in families:
    family_data = df_control[df_control['Family'] == family]
    total_settled = family_data['NoSet'].sum()
    if total_settled > 0:
        main_families.append(family)

print(f"\nFamilies with settlement data: {main_families}")

# Create comparison plot
fig, axes = plt.subplots(len(main_families), 2, figsize=(16, 5*len(main_families)))
if len(main_families) == 1:
    axes = axes.reshape(1, -1)

for idx, family in enumerate(main_families):
    family_data = df_control[df_control['Family'] == family].copy()
    
    # Group by age
    age_summary = family_data.groupby('LarvalAge').agg({
        'NoSet': 'sum',
        'NoAlive': 'sum'
    }).reset_index()
    age_summary['ActualProp'] = age_summary['NoSet'] / age_summary['NoAlive']
    age_summary = age_summary.sort_values('LarvalAge')
    
    ages = age_summary['LarvalAge'].values
    daily_prop = age_summary['ActualProp'].values
    
    # Calculate cumulative settlement
    cumulative_settled = age_summary['NoSet'].cumsum()
    cumulative_total = age_summary['NoAlive'].sum()
    cumulative_prop = cumulative_settled / cumulative_total
    
    # Plot 1: Daily settlement with actual counts
    ax1 = axes[idx, 0]
    bars = ax1.bar(ages, daily_prop, alpha=0.7, color=f'C{idx}', 
                   edgecolor='black', linewidth=1.5, width=0.8)
    
    # Add labels on significant bars
    for age, prop, n_set, n_alive in zip(ages, daily_prop, age_summary['NoSet'], age_summary['NoAlive']):
        if prop > 0.02:  # Label bars above 2%
            ax1.text(age, prop, f'{int(n_set)}/{int(n_alive)}', 
                    ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    ax1.set_xlabel('Larval Age (days)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Daily Settlement Proportion', fontsize=12, fontweight='bold')
    ax1.set_title(f'{family} - Daily Settlement', fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.set_xlim(ages.min()-1, ages.max()+1)
    
    # Plot 2: Cumulative settlement
    ax2 = axes[idx, 1]
    ax2.plot(ages, cumulative_prop, 'o-', color=f'C{idx}', 
            linewidth=3, markersize=8, markeredgecolor='black', markeredgewidth=1.5)
    
    # Add shaded area under curve
    ax2.fill_between(ages, 0, cumulative_prop, alpha=0.3, color=f'C{idx}')
    
    # Mark key points
    first_settle_idx = np.where(cumulative_prop > 0)[0][0]
    ax2.axvline(ages[first_settle_idx], color='red', linestyle='--', 
               linewidth=2, alpha=0.7, label=f'First settlement (day {ages[first_settle_idx]:.0f})')
    
    # Find 50% cumulative if reached
    if cumulative_prop.max() >= 0.5:
        idx_50 = np.where(cumulative_prop >= 0.5)[0][0]
        ax2.axhline(0.5, color='green', linestyle='--', linewidth=2, alpha=0.7)
        ax2.axvline(ages[idx_50], color='green', linestyle='--', 
                   linewidth=2, alpha=0.7, label=f'50% settled (day {ages[idx_50]:.0f})')
    
    ax2.set_xlabel('Larval Age (days)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Cumulative Proportion Settled', fontsize=12, fontweight='bold')
    ax2.set_title(f'{family} - Cumulative Settlement', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10, loc='best')
    ax2.set_xlim(ages.min()-1, ages.max()+1)
    ax2.set_ylim(0, cumulative_prop.max()*1.1)
    
    # Add final cumulative percentage
    final_cum = cumulative_prop.iloc[-1] * 100
    ax2.text(0.98, 0.02, f'Final: {final_cum:.2f}% settled', 
            transform=ax2.transAxes, fontsize=11, fontweight='bold',
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8, 
                     edgecolor='black', linewidth=1.5))

plt.tight_layout()
plt.savefig('settlement_daily_vs_cumulative.png', dpi=300, bbox_inches='tight')
print("✓ Plot saved as 'settlement_daily_vs_cumulative.png'")
plt.close()

print("\n✓ All plots created successfully!")


