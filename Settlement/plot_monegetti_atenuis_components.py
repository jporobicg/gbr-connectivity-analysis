#!/usr/bin/env python3
"""
Plot Monegetti A. tenuis: Survival (decay), Competency, and their product.

Based on Moneghetti et al. (2019) piecewise models:
- Survival: Weibull-Weibull piecewise model
- Competency: Weibull-Exponential piecewise model
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.integrate import quad
import os

sns.set_style("whitegrid")

# ==========================================================================
# PARAMETERS from Moneghetti et al. (2019) — A. tenuis GBR 2012
# ==========================================================================

# Competency model parameters (Weibull-Exponential piecewise)
COMP_PARAMS = {
    'a': 1.290078,
    'b1': 0.00188621,
    'v1': 0.3652014,
    'b2': 0.3969605,
    'tc': 3.331252,
    'Tcp': 69.91391,
}

# Survival/mortality model parameters (Weibull-Weibull piecewise)
SURV_PARAMS = {
    'u1': 0.4002655,
    'v1': 2.892410,
    'u2': 0.01923417,
    'v2': 1.715599,
    'Tcp': 2.583284,
}


def survival_function(t, u1, v1, u2, v2, Tcp):
    """Piecewise Weibull-Weibull survival (probability of being alive at time t)."""
    t = np.atleast_1d(np.float64(t))
    S = np.empty_like(t)
    early = t <= Tcp
    S[early] = np.exp(-(u1 * t[early]) ** v1)
    S[~early] = (np.exp(-(u1 * Tcp) ** v1)
                 * np.exp(-(u2 * (t[~early] - Tcp)) ** v2))
    return S if len(S) > 1 else S[0]


def competency_function(age, a, b1, v1, b2, tc, Tcp):
    """Piecewise Weibull-Exponential competency (P(competent | alive) at time t)."""
    age = np.atleast_1d(np.float64(age))
    comp = np.zeros_like(age)

    for i, t in enumerate(age):
        if t < tc:
            comp[i] = 0.0
        elif t <= Tcp:
            integrand = lambda tau: (a * np.exp(-a * (tau - tc))
                                     * np.exp(-((b1 * (t - tau)) ** v1)))
            comp[i] = max(0.0, min(1.0, quad(integrand, tc, t,
                                             limit=100)[0]))
        else:
            integrand1 = lambda tau: (a * np.exp(-a * (tau - tc))
                                      * np.exp(-((b1 * (Tcp - tau)) ** v1))
                                      * np.exp(-b2 * (t - Tcp)))
            integrand2 = lambda tau: (a * np.exp(-a * (tau - tc))
                                      * np.exp(-b2 * (t - tau)))
            r1 = quad(integrand1, tc, Tcp, limit=100)[0]
            r2 = quad(integrand2, Tcp, t, limit=100)[0]
            comp[i] = max(0.0, min(1.0, r1 + r2))

    return comp if len(comp) > 1 else comp[0]


# ==========================================================================
# LOAD OBSERVED DATA
# ==========================================================================

base = os.path.dirname(os.path.abspath(__file__))

surv_df = pd.read_csv(os.path.join(base, 'monegetti',
                                   'TenuisGBR2012longtermsurvival.csv'))
surv_agg = surv_df.groupby('age (d)').agg(
    larvae=('larvae', 'first'),
    surv_mean=('surv', 'mean'),
).reset_index()
surv_agg['prop_alive'] = surv_agg['surv_mean'] / surv_agg['larvae']

meta_df = pd.read_csv(os.path.join(base, 'monegetti',
                                   'A.tenuisGBR2012metamorphosis.csv'))
meta_agg = meta_df.groupby('age (d)').agg(
    meta=('meta', 'sum'),
    larvae=('larvae', 'sum'),
).reset_index()
meta_agg['prop_meta'] = meta_agg['meta'] / meta_agg['larvae']

# ==========================================================================
# COMPUTE MODEL CURVES
# ==========================================================================

t_max = 85
t_fine = np.linspace(0, t_max, 500)

surv_curve = survival_function(t_fine, **SURV_PARAMS)
comp_curve = competency_function(t_fine, **COMP_PARAMS)
combined_curve = surv_curve * comp_curve

# ==========================================================================
# PLOT
# ==========================================================================

fig, axes = plt.subplots(3, 1, figsize=(8.4, 14), sharey=False)

colors = {
    'model': '#c0392b',
    'data': '#2c3e50',
    'fill': '#e74c3c',
    'combined': '#8e44ad',
}

# --- Panel 1: Survival (decay / mortality) ---
ax = axes[0]
ax.scatter(surv_agg['age (d)'], surv_agg['prop_alive'],
           s=50, color=colors['data'], zorder=3, label='Observed (mean)')
ax.plot(t_fine, surv_curve, lw=2.5, color=colors['model'],
        label='Weibull–Weibull fit')
ax.fill_between(t_fine, 0, surv_curve, alpha=0.12, color=colors['fill'])
ax.axvline(SURV_PARAMS['Tcp'], ls=':', color='grey', lw=1.2,
           label=f"$T_{{cp}}^{{mort}}$ = {SURV_PARAMS['Tcp']:.2f} d")
ax.set_xlabel('Larval age (days)', fontsize=12)
ax.set_ylabel('Proportion alive', fontsize=12)
ax.set_title('(a)  Survival / Decay', fontsize=13, fontweight='bold')
ax.legend(fontsize=9, loc='upper right')
ax.set_xlim(0, t_max)
ax.set_ylim(-0.02, 1.05)

param_txt = (f"$u_1$={SURV_PARAMS['u1']:.4f}  $v_1$={SURV_PARAMS['v1']:.3f}\n"
             f"$u_2$={SURV_PARAMS['u2']:.5f}  $v_2$={SURV_PARAMS['v2']:.3f}\n"
             f"$T_{{cp}}$={SURV_PARAMS['Tcp']:.3f} d")
ax.text(0.97, 0.55, param_txt, transform=ax.transAxes, fontsize=8.5,
        va='top', ha='right', family='monospace',
        bbox=dict(boxstyle='round', fc='wheat', alpha=0.85))

# --- Panel 2: Competency (given alive) ---
ax = axes[1]
ax.scatter(meta_agg['age (d)'], meta_agg['prop_meta'],
           s=50, color=colors['data'], zorder=3, label='Observed (metamorphosis)')
ax.plot(t_fine, comp_curve, lw=2.5, color=colors['model'],
        label='Weibull–Exp fit')
ax.fill_between(t_fine, 0, comp_curve, alpha=0.12, color=colors['fill'])
ax.axvline(COMP_PARAMS['tc'], ls=':', color='grey', lw=1.2,
           label=f"$t_c$ = {COMP_PARAMS['tc']:.2f} d")
ax.axvline(COMP_PARAMS['Tcp'], ls='--', color='grey', lw=1.2,
           label=f"$T_{{cp}}^{{comp}}$ = {COMP_PARAMS['Tcp']:.2f} d")
ax.set_xlabel('Larval age (days)', fontsize=12)
ax.set_ylabel('P(competent | alive)', fontsize=12)
ax.set_title('(b)  Competency', fontsize=13, fontweight='bold')
ax.legend(fontsize=9, loc='upper right')
ax.set_xlim(0, t_max)
ax.set_ylim(-0.02, 1.05)

param_txt = (f"$a$={COMP_PARAMS['a']:.4f}   $b_1$={COMP_PARAMS['b1']:.6f}\n"
             f"$v_1$={COMP_PARAMS['v1']:.4f}  $b_2$={COMP_PARAMS['b2']:.4f}\n"
             f"$t_c$={COMP_PARAMS['tc']:.3f} d  $T_{{cp}}$={COMP_PARAMS['Tcp']:.2f} d")
ax.text(0.97, 0.55, param_txt, transform=ax.transAxes, fontsize=8.5,
        va='top', ha='right', family='monospace',
        bbox=dict(boxstyle='round', fc='wheat', alpha=0.85))

# --- Panel 3: Combined = Survival × Competency ---
ax = axes[2]
ax.plot(t_fine, combined_curve, lw=2.5, color=colors['combined'],
        label='Survival × Competency')
ax.fill_between(t_fine, 0, combined_curve, alpha=0.15, color=colors['combined'])
ax.plot(t_fine, surv_curve, lw=1.2, ls='--', color=colors['model'],
        alpha=0.5, label='Survival only')
ax.plot(t_fine, comp_curve, lw=1.2, ls=':', color=colors['model'],
        alpha=0.5, label='Competency only')
ax.set_xlabel('Larval age (days)', fontsize=12)
ax.set_ylabel('P(alive & competent)', fontsize=12)
ax.set_title('(c)  Survival × Competency', fontsize=13, fontweight='bold')
ax.legend(fontsize=9, loc='upper right')
ax.set_xlim(0, t_max)
ax.set_ylim(-0.02, 1.05)

peak_idx = np.argmax(combined_curve)
peak_age = t_fine[peak_idx]
peak_val = combined_curve[peak_idx]
ax.annotate(f'Peak = {peak_val:.2f}\nat {peak_age:.1f} d',
            xy=(peak_age, peak_val),
            xytext=(peak_age + 8, peak_val + 0.12),
            fontsize=9, ha='left',
            arrowprops=dict(arrowstyle='->', color='black', lw=1.2),
            bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.9))

fig.suptitle('Moneghetti et al. (2019)  —  A. tenuis  GBR 2012',
             fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout()

os.makedirs(os.path.join(base, 'figures'), exist_ok=True)
out_path = os.path.join(base, 'figures', 'monegetti_atenuis_survival_competency.png')
plt.savefig(out_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {out_path}")
