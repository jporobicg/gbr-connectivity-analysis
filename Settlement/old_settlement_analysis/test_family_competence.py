"""
Test script to visualize family-based competence functions.

This script generates plots comparing competence curves across different coral families
based on the Randall et al. (2024) paper.
"""

import numpy as np
import matplotlib.pyplot as plt
from ecological_processes import family_based_competence

# Define age range (0-30 days)
ages = np.linspace(0, 30, 300)

# Define families to compare
families = ['acroporidae', 'merulinidae', 'pocilloporidae', 'poritidae', 'fungiidae']

# Create figure
fig, axes = plt.subplots(2, 1, figsize=(12, 10))

# Plot 1: All families together
ax1 = axes[0]
for family in families:
    competence = family_based_competence(ages, family=family)
    ax1.plot(ages, competence, linewidth=2, label=family.capitalize())

ax1.set_xlabel('Larval Age (days)', fontsize=12)
ax1.set_ylabel('Competence (0-1)', fontsize=12)
ax1.set_title('Comparison of Family-Based Competence Functions\n(Based on Randall et al. 2024)', 
              fontsize=14, fontweight='bold')
ax1.legend(loc='upper right', fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 30)
ax1.set_ylim(0, 1.05)

# Plot 2: Individual family subplots
families_subset = ['acroporidae', 'merulinidae', 'pocilloporidae']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

ax2 = axes[1]
for i, family in enumerate(families_subset):
    competence = family_based_competence(ages, family=family)
    ax2.plot(ages, competence, linewidth=2.5, label=family.capitalize(), 
             color=colors[i], alpha=0.8)
    
    # Mark key points
    # Find precompetency and peak
    first_nonzero = np.where(competence > 0.01)[0]
    if len(first_nonzero) > 0:
        precomp_idx = first_nonzero[0]
        ax2.axvline(ages[precomp_idx], color=colors[i], linestyle='--', alpha=0.3)
    
    peak_idx = np.argmax(competence)
    ax2.plot(ages[peak_idx], competence[peak_idx], 'o', color=colors[i], 
             markersize=8, markeredgecolor='black', markeredgewidth=1)

ax2.set_xlabel('Larval Age (days)', fontsize=12)
ax2.set_ylabel('Competence (0-1)', fontsize=12)
ax2.set_title('Key GBR Families: Acroporidae, Merulinidae, and Pocilloporidae', 
              fontsize=14, fontweight='bold')
ax2.legend(loc='upper right', fontsize=10)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, 25)
ax2.set_ylim(0, 1.05)

# Add annotations
ax2.text(0.98, 0.05, 
         'Dashed lines: Precompetency period\nCircles: Peak competence', 
         transform=ax2.transAxes, fontsize=9, verticalalignment='bottom',
         horizontalalignment='right', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('family_competence_comparison.png', dpi=300, bbox_inches='tight')
print("✓ Plot saved as 'family_competence_comparison.png'")
plt.close()

# Print summary statistics
print("\n" + "="*70)
print("FAMILY-BASED COMPETENCE SUMMARY")
print("="*70)
print("\nBased on: Randall et al. (2024) Communications Biology")
print("Paper: Larval precompetency and settlement behaviour in 25 Indo-Pacific coral species")
print("\n" + "-"*70)

for family in families:
    competence = family_based_competence(ages, family=family)
    
    # Find key statistics
    first_competent = np.where(competence > 0.01)[0]
    precomp_age = ages[first_competent[0]] if len(first_competent) > 0 else np.nan
    
    peak_idx = np.argmax(competence)
    peak_age = ages[peak_idx]
    peak_comp = competence[peak_idx]
    
    # Find age at 50% of peak
    half_peak = peak_comp * 0.5
    above_half = np.where(competence >= half_peak)[0]
    if len(above_half) > 0:
        age_50_idx = above_half[0]
        age_50 = ages[age_50_idx]
    else:
        age_50 = np.nan
    
    # Find competence window (competence > 10% of peak)
    window_mask = competence > (peak_comp * 0.1)
    window_ages = ages[window_mask]
    window_duration = window_ages[-1] - window_ages[0] if len(window_ages) > 0 else 0
    
    print(f"\n{family.upper()}")
    print(f"  Precompetency period:  {precomp_age:.1f} days")
    print(f"  Age at 50% competence: {age_50:.1f} days")
    print(f"  Peak competence age:   {peak_age:.1f} days (value: {peak_comp:.3f})")
    print(f"  Competence window:     {window_duration:.1f} days (>10% of peak)")

print("\n" + "="*70)
print("\nECOLOGICAL INTERPRETATION:")
print("-"*70)
print("• Acroporidae:    Early settlers, short dispersal (fast-growing, competitive)")
print("• Merulinidae:    Late settlers, longer dispersal (stress-tolerant, hardy)")
print("• Pocilloporidae: Very early settlers (often brooders, short dispersal)")
print("• Poritidae:      Latest settlers, extended competency (slow-growing, persistent)")
print("• Fungiidae:      Intermediate dispersers (solitary strategy)")
print("="*70)
print("\n✓ Test completed successfully!\n")


