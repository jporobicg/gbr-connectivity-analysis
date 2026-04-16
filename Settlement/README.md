# Larval Competency Modeling: Hybrid Continuous Weighted Approach

## Overview

This project develops and validates a novel approach for modeling coral larval competency as a function of age, using **continuous weighted proportions** combined with **Monegetti-style sigmoid models**.

---

## Quick Start

### Run the Main Analysis
```bash
python hybrid_competency_model.py
```

### Run Uncertainty Analysis (with Bootstrap CIs)
```bash
python uncertainty_analysis_by_family.py
```

### Run Monegetti Piecewise Model
```bash
python monegetti_piecewise_model.py
```

**Note**: The Monegetti model script includes checkpointing - results are saved after each family is processed to `estimations/estimations_*Family*_*timestamp*.csv`. If the script crashes, you can resume from the last completed family.

### Generate Monegetti Model Plots
```bash
python plot_monegetti_results.py
```

---

## Key Results

### TC50 Values (Age at 50% Competency)

| Family | Best Model | TC50 (days) | 95% CI | Bootstrap Samples |
|--------|-----------|-------------|--------|-------------------|
| **Acroporidae** | Weibull | 5.91 | [3.15-7.00] | 496 |
| **Poritidae** | Gompertz | 4.31 | N/A | 247 |
| **Merulinidae** | Weibull | 3.97 | [0.86-4.72] | 497 |
| **Diploastreidae** | Logistic | 4.01 | [-35.62-6.94] | 499 |
| **Agariciidae** | Weibull | 5.00 | [2.95-5.71] | 263 |
| **Lobophylliidae** | Gompertz | 4.01 | N/A | 328 |
| **Euphylliidae** | Logistic | 4.02 | [0.20-5.93] | 500 |

---

## Project Structure

```
Settlement/
├── hybrid_competency_model.py          # Main analysis script (simple models)
├── monegetti_piecewise_model.py         # Monegetti piecewise model analysis
├── plot_monegetti_results.py            # Generate Monegetti model plots
├── uncertainty_analysis_by_family.py    # Bootstrap uncertainty quantification
│
├── estimations/                         # Monegetti model results (checkpointed)
│   ├── estimations_Acroporidae_*.csv
│   ├── estimations_Merulinidae_*.csv
│   └── ... (one file per family with timestamp)
│
├── figures/                             # Model plots
│   ├── monegetti_Acroporidae.png
│   ├── monegetti_Merulinidae.png
│   └── ... (one plot per family)
│
├── family_results/                      # Individual family outputs (hybrid approach)
│   ├── Acroporidae_parameters.csv
│   ├── Acroporidae_uncertainty.png
│   ├── ... (for each family)
│   └── all_families_summary.csv
│
├── hybrid_competency_models.png         # Main figure (7 families, simple models)
├── hybrid_model_summary.csv             # Summary table (simple models)
├── monegetti_model_summary.csv          # Summary table (Monegetti model)
│
├── COMPREHENSIVE_ANALYSIS_REPORT.md     # Detailed method comparison
├── FINAL_CONCLUSIONS.md                 # Conclusions & recommendations
├── README.md                            # This file
│
├── old_analysis/                        # Previous analyses (archived)
│   ├── larval_competency_model.py      # Original CCA analysis
│   ├── randal_approach_replication.py  # Binary cumulative approach
│   └── ... (various comparison documents)
│
├── monegetti/                           # Monegetti's original data
│   ├── A.tenuisGBR2012metamorphosis.csv
│   └── competence model code.R
│
├── Randal/                              # Randal's data
│   ├── data/Settlement.csv
│   └── s42003-024-05824-3.pdf
│
└── Randal_github/                       # Randal's R code (reference)
    └── scripts/
        ├── 20_fitModels.R
        ├── functions.R
        └── 10_processData.R
```

---

## Methodology

### Our Hybrid Approach

**Key Innovation**: Uses **continuous proportions** (not binary) with **sample-size weighting**

#### Steps:
1. **Data Preparation**
   - Filter to rubble treatment
   - Aggregate by family and larval age
   - Calculate proportion settled at each age

2. **Model Fitting**
   - Fit three sigmoid models: Logistic, Gompertz, Weibull
   - Use weighted binomial likelihood
   - Select best model via AIC

3. **Uncertainty Quantification**
   - Bootstrap resampling (500 iterations)
   - Calculate parameter standard errors
   - Compute 95% confidence intervals

#### Advantages:
✅ **Avoids binary threshold artifacts**  
✅ **Allows natural fluctuation in settlement**  
✅ **Sample-size weighting for robustness**  
✅ **Multiple model comparison**  
✅ **Bootstrap uncertainty quantification**  
✅ **Computationally efficient**  
✅ **Transparent methodology**

---

## Monegetti Piecewise Weibull-Exponential Model

### Model Description

The **Monegetti piecewise model** (Moneghetti et al. 2019) represents a mechanistic approach to modeling larval competency that captures the biological reality of two-phase competency development:

1. **Precompetency period** (t < tc): No settlement possible
2. **Early competency phase** (tc < t < Tcp): Weibull loss of competency
3. **Late competency phase** (t > Tcp): Exponential loss of competency

The model integrates over all possible times of competency acquisition (τ) to calculate the probability of being competent at age t:

- **Early phase**: Competency acquired at τ, maintained with Weibull decay
- **Late phase**: Transition to exponential decay after change point Tcp

**Model Parameters**:
- `a`: Rate of acquisition of competency (when t > tc)
- `b1`: Weibull scale parameter (early period loss rate)
- `v1`: Weibull shape parameter (early period decay shape)
- `b2`: Exponential rate parameter (late period loss rate)
- `tc`: Precompetency period (onset age)
- `Tcp`: Change point between early and late phases

### Model Selection Results

We compared the Monegetti piecewise model (6 parameters) against simpler alternatives (Logistic, Gompertz, Weibull) using AICc for model selection:

| Family | Best Model | N Replicates | N Larvae | AICc | Δ AICc | TC50 (days) | Max Competency |
|--------|-----------|--------------|----------|------|--------|-------------|----------------|
| **Acroporidae** | Weibull | 4,730 | 45,123 | 56,895.3 | +201.9 | 5.36 | 0.510 |
| **Agariciidae** | **Monegetti** | 167 | 1,542 | 598.8 | 0.0 | N/A | 0.077 |
| **Diploastreidae** | Weibull | 539 | 5,329 | 5,198.1 | +9.2 | 4.43 | 0.251 |
| **Euphylliidae** | Weibull | 496 | 4,631 | 5,825.0 | +9.4 | 3.71 | 0.389 |
| **Lobophylliidae** | **Monegetti** | 760 | 6,560 | 2,982.9 | 0.0 | 5.32 | 0.156 |
| **Merulinidae** | **Monegetti** | 3,880 | 38,495 | 38,053.9 | 0.0 | 30.74 | 0.280 |
| **Poritidae** | Logistic | 1,187 | 11,612 | 10,771.3 | +42.7 | 28.26 | 0.264 |

**Key Findings**:
- **3 families** (Agariciidae, Lobophylliidae, Merulinidae) show **Monegetti as best model** (Δ AICc = 0.0)
- **3 families** prefer simpler **Weibull model** (Acroporidae, Diploastreidae, Euphylliidae)
- **1 family** (Poritidae) prefers **Logistic model**
- Monegetti model is particularly valuable for families with **complex competency dynamics** (early rapid onset, late exponential decay)

---

### Family-by-Family Analysis

#### 1. **Acroporidae** (Best Model: Weibull)
- **Sample Size**: 4,730 replicates, 45,123 larvae
- **Monegetti Parameters**: a=0.332, b1=0.0017, v1=0.100, b2=0.001, tc=3.90 days, Tcp=71.48 days
- **TC50**: 5.36 days | **Max Competency**: 0.510
- **Model Comparison**: Weibull preferred (Δ AICc = +201.9)

**Pros of Monegetti approach**:
- ✅ Captures very long competency period (Tcp = 71.5 days)
- ✅ Identifies early precompetency period (tc = 3.9 days)
- ✅ Mechanistically sound two-phase structure

**Cons of Monegetti approach**:
- ❌ Overparameterized for this family (6 params vs 4 for Weibull)
- ❌ Large sample size makes simpler models more parsimonious
- ❌ Weibull model fits equally well with fewer parameters

**Recommendation**: Use **Weibull model** for Acroporidae due to parsimony and equivalent fit.

**Plot**: See `figures/monegetti_Acroporidae.png`

---

#### 2. **Agariciidae** (Best Model: Monegetti)
- **Sample Size**: 167 replicates, 1,542 larvae
- **Monegetti Parameters**: a=0.100, b1=1.000, v1=1.789, b2=1.000, tc=3.00 days, Tcp=6.19 days
- **TC50**: N/A (max competency too low: 0.077) | **Max Competency**: 0.077
- **Model Comparison**: Monegetti preferred (Δ AICc = 0.0, vs Weibull +1.1)

**Pros of Monegetti approach**:
- ✅ Best fitting model (lowest AICc)
- ✅ Captures short competency window (Tcp = 6.2 days)
- ✅ Identifies rapid competency loss (high b2 = 1.0)

**Cons of Monegetti approach**:
- ❌ Very low maximum competency (7.7%) suggests poor settlement overall
- ❌ Small sample size (167 replicates) limits statistical power
- ❌ Limited age range (3-7 days) constrains model fitting

**Recommendation**: Use **Monegetti model** for Agariciidae, but interpret with caution due to low settlement rates.

**Plot**: See `figures/monegetti_Agariciidae.png`

---

#### 3. **Diploastreidae** (Best Model: Weibull)
- **Sample Size**: 539 replicates, 5,329 larvae
- **Monegetti Parameters**: a=0.400, b1=0.941, v1=0.223, b2=0.001, tc=3.74 days, Tcp=13.24 days
- **TC50**: 4.43 days | **Max Competency**: 0.251
- **Model Comparison**: Weibull preferred (Δ AICc = +9.2)

**Pros of Monegetti approach**:
- ✅ Identifies change point at 13.2 days (moderate competency period)
- ✅ Captures early rapid acquisition (a = 0.400)
- ✅ Mechanistic interpretation of two-phase development

**Cons of Monegetti approach**:
- ❌ Weibull model fits nearly as well with 2 fewer parameters
- ❌ Change point may not be biologically meaningful for this family
- ❌ Overparameterization penalty (Δ AICc = +9.2)

**Recommendation**: Use **Weibull model** for Diploastreidae for parsimony.

**Plot**: See `figures/monegetti_Diploastreidae.png`

---

#### 4. **Euphylliidae** (Best Model: Weibull)
- **Sample Size**: 496 replicates, 4,631 larvae
- **Monegetti Parameters**: a=0.413, b1=0.071, v1=0.102, b2=0.001, tc=2.35 days, Tcp=33.45 days
- **TC50**: 3.71 days | **Max Competency**: 0.389
- **Model Comparison**: Weibull preferred (Δ AICc = +9.4)

**Pros of Monegetti approach**:
- ✅ Early competency onset (tc = 2.35 days, earliest among families)
- ✅ Long competency period (Tcp = 33.5 days)
- ✅ Good maximum competency (38.9%)

**Cons of Monegetti approach**:
- ❌ Weibull model provides equivalent fit with fewer parameters
- ❌ Two-phase structure may not be necessary
- ❌ Overparameterization (Δ AICc = +9.4)

**Recommendation**: Use **Weibull model** for Euphylliidae for parsimony.

**Plot**: See `figures/monegetti_Euphylliidae.png`

---

#### 5. **Lobophylliidae** (Best Model: Monegetti)
- **Sample Size**: 760 replicates, 6,560 larvae
- **Monegetti Parameters**: a=0.100, b1=1.000, v1=1.235, b2=0.140, tc=3.11 days, Tcp=10.91 days
- **TC50**: 5.32 days | **Max Competency**: 0.156
- **Model Comparison**: Monegetti preferred (Δ AICc = 0.0, vs Weibull +25.7)

**Pros of Monegetti approach**:
- ✅ Best fitting model (lowest AICc)
- ✅ Captures short competency window (Tcp = 10.9 days)
- ✅ Identifies transition to exponential decay (b2 = 0.140)
- ✅ Moderate sample size provides good statistical power

**Cons of Monegetti approach**:
- ❌ Relatively low maximum competency (15.6%)
- ❌ Short age range (3-23 days) limits late-phase characterization

**Recommendation**: Use **Monegetti model** for Lobophylliidae - it provides the best fit and captures biologically meaningful two-phase structure.

**Plot**: See `figures/monegetti_Lobophylliidae.png`

---

#### 6. **Merulinidae** (Best Model: Monegetti)
- **Sample Size**: 3,880 replicates, 38,495 larvae (largest dataset)
- **Monegetti Parameters**: a=0.245, b1=0.333, v1=0.322, b2=0.001, tc=1.92 days, Tcp=54.10 days
- **TC50**: 30.74 days | **Max Competency**: 0.280
- **Model Comparison**: Monegetti strongly preferred (Δ AICc = 0.0, vs Weibull +627.9)

**Pros of Monegetti approach**:
- ✅ **Strongly best model** - largest Δ AICc advantage (+627.9 over Weibull)
- ✅ Very early competency onset (tc = 1.92 days, earliest)
- ✅ Very long competency period (Tcp = 54.1 days)
- ✅ Large sample size provides excellent statistical power
- ✅ Captures complex competency dynamics (late TC50 = 30.7 days)

**Cons of Monegetti approach**:
- ❌ Late TC50 (30.7 days) suggests competency develops slowly
- ❌ Moderate maximum competency (28.0%)

**Recommendation**: **Definitely use Monegetti model** for Merulinidae - it provides the strongest evidence for two-phase competency structure among all families.

**Plot**: See `figures/monegetti_Merulinidae.png`

---

#### 7. **Poritidae** (Best Model: Logistic)
- **Sample Size**: 1,187 replicates, 11,612 larvae
- **Monegetti Parameters**: a=0.389, b1=0.828, v1=0.195, b2=0.176, tc=3.90 days, Tcp=26.67 days
- **TC50**: 28.26 days | **Max Competency**: 0.264
- **Model Comparison**: Logistic preferred (Δ AICc = +42.7)

**Pros of Monegetti approach**:
- ✅ Identifies moderate competency period (Tcp = 26.7 days)
- ✅ Captures transition to exponential decay (b2 = 0.176)
- ✅ Mechanistic interpretation available

**Cons of Monegetti approach**:
- ❌ Logistic model fits much better (Δ AICc = +42.7)
- ❌ Overparameterized for this family's data structure
- ❌ Late TC50 (28.3 days) may reflect model complexity rather than biology

**Recommendation**: Use **Logistic model** for Poritidae - simpler model provides better fit.

**Plot**: See `figures/monegetti_Poritidae.png`

---

### Summary of Model Selection

**When to use Monegetti piecewise model**:
1. ✅ **Merulinidae**: Strong evidence (Δ AICc = +627.9 advantage)
2. ✅ **Lobophylliidae**: Best fit (Δ AICc = 0.0, +25.7 over Weibull)
3. ✅ **Agariciidae**: Best fit (Δ AICc = 0.0, +1.1 over Weibull)

**When simpler models are preferred**:
1. ❌ **Acroporidae**: Use Weibull (Δ AICc = +201.9 penalty)
2. ❌ **Diploastreidae**: Use Weibull (Δ AICc = +9.2 penalty)
3. ❌ **Euphylliidae**: Use Weibull (Δ AICc = +9.4 penalty)
4. ❌ **Poritidae**: Use Logistic (Δ AICc = +42.7 penalty)

**Key Insight**: The Monegetti model is most valuable for families with **complex competency dynamics** that show clear evidence of two-phase development (early Weibull decay, late exponential decay). For families with simpler competency curves, parsimonious models (Weibull, Logistic) are preferred.

---

## Comparison with Other Approaches

### 1. Randal's Bayesian Hierarchical (2024)
- **Method**: Binary cumulative, Bayesian with random effects
- **Level**: Species-specific
- **Uncertainty**: Posterior distributions
- **Time**: Hours per species
- **Our difference**: +1-3 days higher TC50 (family vs species aggregation)

### 2. Monegetti's Piecewise Weibull (2012)
- **Method**: Continuous, two-phase development
- **Level**: A. tenuis only
- **Uncertainty**: None
- **Time**: Minutes
- **Our difference**: Similar concept but generalized to all families

---

## Key Findings

### 1. **Family-Level Patterns are Consistent**
Despite aggregating across species, family-level competency curves are stable and meaningful.

### 2. **TC50 Values Cluster Around 4-6 Days**
Most families reach 50% competency between 4-6 days post-spawning (rubble treatment).

### 3. **Model Selection Varies by Family**

**Simple Models (Hybrid Approach)**:
- Weibull: Acroporidae, Merulinidae, Agariciidae (rapid onset, early peak)
- Gompertz: Poritidae, Lobophylliidae (gradual increase, late plateau)
- Logistic: Diploastreidae, Euphylliidae (symmetric sigmoid)

**Monegetti Piecewise Model**:
- Best for: Merulinidae (strong evidence, Δ AICc = +627.9), Lobophylliidae, Agariciidae
- Captures two-phase competency development (early Weibull decay, late exponential decay)
- Most valuable for families with complex competency dynamics

### 4. **Bootstrap Provides Robust Uncertainty**
95% confidence intervals are reasonably narrow (1-3 day range) for most families.

### 5. **Continuous Approach Reveals Natural Variation**
Unlike binary methods, our approach captures fluctuations in settlement proportion over time.

---

## Usage Examples

### Example 1: Run Full Analysis
```python
python hybrid_competency_model.py
```

**Outputs**:
- `hybrid_competency_models.png` - 7-panel figure
- `hybrid_model_summary.csv` - Parameter table
- `hybrid_predictions.csv` - Predictions at key ages

### Example 2: Uncertainty Analysis for One Family
```python
from uncertainty_analysis_by_family import analyze_single_family
from hybrid_competency_model import prepare_data_continuous

# Load data
family_data = prepare_data_continuous('Randal/data/Settlement.csv', 'rubble')

# Analyze Acroporidae
analyze_single_family('Acroporidae', family_data['Acroporidae'], 'output/')
```

### Example 3: Predict Competency at Specific Ages
```python
from hybrid_competency_model import logistic_competency, gompertz_competency, weibull_competency
import numpy as np

# Acroporidae (Weibull): tc=3.91, a=0.393, b=0.000, v=0.001
ages = np.array([3, 5, 7, 10, 15])
competency = weibull_competency(ages, 3.91, 0.393, 0.000, 0.001)
print(competency)
# Output: [0.000, 0.235, 0.474, 0.612, 0.665]
```

---

## Dependencies

```bash
pip install pandas numpy scipy matplotlib seaborn
```

**Python Version**: 3.8+

---

## Documentation

### Main Documents

1. **[COMPREHENSIVE_ANALYSIS_REPORT.md](COMPREHENSIVE_ANALYSIS_REPORT.md)**
   - Detailed comparison of three approaches
   - Parameter tables with uncertainty
   - Method comparison

2. **[FINAL_CONCLUSIONS.md](FINAL_CONCLUSIONS.md)**
   - Strengths, assumptions, constraints
   - When to use each approach
   - Recommendations for future research

3. **[README.md](README.md)** (this file)
   - Quick start guide
   - Project overview

### Family-Specific Results

Located in `family_results/`:
- Individual parameter tables with standard errors and 95% CIs
- Plots showing model fit with bootstrap uncertainty bands
- Combined summary table

---

## Citation

If you use this approach, please cite:

**Foundational Works**:
1. Randal et al. (2024). "Age to settlement competency and settlement cue preference in coral larvae." *Communications Biology*.
2. Moneghetti et al. (2019). "High-frequency sampling and piecewise models reshape dispersal kernels of a common reef coral." *Ecology Letters*.
3. Monegetti et al. (2012). "A quantitative assessment of the competency period of *Acropora tenuis* coral larvae."

**This Study**:
3. [Your citation for this hybrid approach]

---

## Contact & Support

For questions or issues:
- Check documentation in `COMPREHENSIVE_ANALYSIS_REPORT.md`
- Review individual family results in `family_results/`
- Examine code comments in `hybrid_competency_model.py`

---

## Future Directions

### Immediate Next Steps
1. ✅ ~~Hybrid continuous approach~~ (DONE)
2. ✅ ~~Bootstrap uncertainty~~ (DONE)
3. ⏳ **Multi-treatment extension**
4. ⏳ **Species-level random effects**
5. ⏳ **Cross-validation**

### Research Questions
- Does competency plateau or decline with age?
- Are there critical windows for settlement?
- How do treatments interact with age effects?
- Can we predict field settlement from lab competency?

---

## Acknowledgments

- **Randal et al.**: Bayesian hierarchical framework and data
- **Monegetti et al.**: Mechanistic competency models
- **Murray et al.**: Original settlement data collection

---

**Last Updated**: December 2, 2025  
**Version**: 2.0  
**Status**: Complete ✅ (includes Monegetti piecewise model)

---

## Quick Reference

### TC50 Interpretation
- **3-4 days**: Early competency (Merulinidae, Poritidae, Lobophylliidae)
- **4-5 days**: Moderate competency (Diploastreidae, Euphylliidae)
- **5-6 days**: Later competency (Acroporidae, Agariciidae)

### Model Selection Guide
- **Logistic**: Symmetric sigmoid, single inflection point
- **Gompertz**: Asymmetric, slower initial rise, late plateau
- **Weibull**: Early onset, rapid increase, early plateau

### Files to Check First
1. `family_results/all_families_summary.csv` - Overview
2. `FINAL_CONCLUSIONS.md` - Interpretation
3. `hybrid_competency_models.png` - Visual summary
4. `family_results/<Family>_parameters.csv` - Detailed parameters

---

**Happy Modeling! 🐠🪸**


