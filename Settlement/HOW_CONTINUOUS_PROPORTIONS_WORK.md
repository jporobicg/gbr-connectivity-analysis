# How Continuous Proportions Are Calculated

## Overview

The term "continuous proportions" refers to two related but distinct concepts in our analysis:

1. **Model fitting approach**: Using raw counts (binomial likelihood) rather than pre-aggregated proportions
2. **Model predictions**: Smooth, continuous curves over age ranges

## 1. Discrete Observed Proportions (Raw Data)

### At the Replicate Level

Each row in the data represents a single replicate (experimental unit):

```python
# From utility_tools.py, load_settlement_data()
df['PropSettled'] = df['NoSet'] / df['NoAlive']
```

**Example:**
- Replicate 1: Age = 5 days, NoSet = 8, NoAlive = 10 → PropSettled = 0.8 (80%)
- Replicate 2: Age = 5 days, NoSet = 7, NoAlive = 10 → PropSettled = 0.7 (70%)
- Replicate 3: Age = 5 days, NoSet = 9, NoAlive = 10 → PropSettled = 0.9 (90%)

**Key point**: These are **discrete** data points - one proportion per replicate.

## 2. Model Fitting: Using Raw Counts (Not Aggregated Proportions)

### Why "Continuous" Approach?

The term "continuous" here means we **don't pre-aggregate** the data into proportions before fitting. Instead, we use **raw counts** with binomial likelihood.

### Code Implementation

```python
# From utility_tools.py, neg_log_likelihood_binomial()
def neg_log_likelihood_binomial(
    params: np.ndarray,
    ages: np.ndarray,           # Individual age values (one per replicate)
    n_settled: np.ndarray,      # NoSet for each replicate
    n_total: np.ndarray,        # NoAlive for each replicate
    model_func,
    phi: float = 1.0
) -> float:
    # Predict competency for each observation
    pred_comp = model_func(ages, *params)  # Continuous prediction
    pred_comp = np.clip(pred_comp, 1e-10, 1 - 1e-10)
    
    # Binomial log-likelihood using RAW COUNTS
    ll = np.sum(
        n_settled * np.log(pred_comp) +           # Successes
        (n_total - n_settled) * np.log(1 - pred_comp)  # Failures
    )
    
    return -ll
```

### What This Means

**For each replicate individually:**
- Age: 5.0 days
- NoSet: 8 (observed successes)
- NoAlive: 10 (total trials)
- Model prediction: 0.75 (predicted proportion)
- Likelihood contribution: `8 * log(0.75) + 2 * log(0.25)`

**The model sees:**
- All replicates at their individual ages
- Raw counts (NoSet, NoAlive) for each replicate
- **NOT** pre-aggregated proportions

### Why This Is "Continuous"

1. **No aggregation before fitting**: We preserve all information from individual replicates
2. **Model can predict any age**: The model function is continuous (can evaluate at any age)
3. **Sample-size weighting**: Larger replicates (more larvae) automatically have more influence

### Comparison: Continuous vs. Aggregated Approach

#### ❌ Aggregated Approach (NOT used):
```python
# BAD: Aggregate first, then fit to proportions
age_summary = df.groupby('LarvalAge').agg({
    'NoSet': 'sum',
    'NoAlive': 'sum'
})
age_summary['PropSettled'] = age_summary['NoSet'] / age_summary['NoAlive']
# Fit model to age_summary['PropSettled'] - loses sample size information!
```

**Problems:**
- Loses information about sample sizes per replicate
- Treats n=10 and n=100 equally
- Can't account for replicate-level variation

#### ✅ Continuous Approach (USED):
```python
# GOOD: Fit to raw counts, model predicts continuous curve
# Model sees: (age=5.0, NoSet=8, NoAlive=10), (age=5.0, NoSet=7, NoAlive=10), ...
# Model predicts: competency(age) = continuous function
# Likelihood: uses raw counts with proper binomial weighting
```

**Benefits:**
- Preserves all information
- Properly weights by sample size
- Accounts for replicate-level variation
- Model can predict at any age (continuous)

## 3. Continuous Model Predictions (Smooth Curves)

### For Plotting

When we want to plot the fitted model, we generate a **smooth, continuous curve**:

```python
# From figures_analysis.py, create_family_figure()
age_range = np.linspace(0, ages.max() + 5, 300)  # 300 points from 0 to max_age+5
if best_model == 'Monegetti':
    pred = np.array([model_func(age, *params) for age in age_range])
else:
    pred = model_func(age_range, *params)
```

**This creates:**
- 300 evenly-spaced age points (e.g., 0.0, 0.1, 0.2, ..., 15.0 days)
- Model predictions at each point (e.g., competency = 0.0, 0.01, 0.02, ..., 0.35)
- A **smooth, continuous curve** for plotting

### Example

**Discrete observations** (what we measured):
```
Age (days)  |  Observed Proportion
------------|-------------------
5.0         |  0.80
5.0         |  0.70
5.0         |  0.90
6.0         |  0.85
6.0         |  0.90
...
```

**Continuous model prediction** (what the model predicts):
```
Age (days)  |  Predicted Competency
------------|-------------------
0.0         |  0.000
0.1         |  0.001
0.2         |  0.002
...
5.0         |  0.750  ← Model prediction at age 5
5.1         |  0.755
5.2         |  0.760
...
15.0        |  0.350
```

## 4. Aggregation for Visualization Only

### When We Aggregate

We **only aggregate** when creating plots to show observed data:

```python
# From figures_analysis.py, create_family_figure()
# Aggregate data by age FOR PLOTTING ONLY
age_summary = df_family.groupby('LarvalAge').agg({
    'NoSet': 'sum',
    'NoAlive': 'sum'
}).reset_index()
age_summary['PropSettled'] = age_summary['NoSet'] / age_summary['NoAlive']
```

**This aggregation:**
- ✅ **Only for visualization** (error bars, plotting)
- ❌ **NOT used for model fitting**
- Shows average observed proportion at each age
- Includes error bars based on sample size

### Why Aggregate for Plotting?

1. **Clarity**: Can't plot 1000+ individual replicate points
2. **Error bars**: Need to calculate standard errors from aggregated data
3. **Comparison**: Easy to compare observed (aggregated) vs. predicted (continuous curve)

## 5. Complete Workflow

### Step-by-Step

1. **Load data** (replicate level):
   ```python
   df = load_settlement_data('Settlement.csv')
   # Each row: (LarvalAge, NoSet, NoAlive, Family, ...)
   # PropSettled = NoSet / NoAlive (for each replicate)
   ```

2. **Filter by family**:
   ```python
   df_family = df[df['Family'] == 'Acroporidae']
   # Still at replicate level
   ```

3. **Fit model** (uses raw counts):
   ```python
   # Model fitting function receives:
   ages = df_family['LarvalAge'].values      # [5.0, 5.0, 5.0, 6.0, ...]
   n_settled = df_family['NoSet'].values    # [8, 7, 9, 8, ...]
   n_total = df_family['NoAlive'].values    # [10, 10, 10, 10, ...]
   
   # Model predicts continuous function: competency(age)
   # Likelihood uses raw counts: binomial(n_total, competency(age))
   ```

4. **Generate predictions** (continuous curve):
   ```python
   age_range = np.linspace(0, 15, 300)  # 300 points
   predictions = model_func(age_range, *fitted_params)  # Continuous curve
   ```

5. **Plot** (aggregate observations for visualization):
   ```python
   # Aggregate for plotting
   age_summary = df_family.groupby('LarvalAge').agg({
       'NoSet': 'sum',
       'NoAlive': 'sum'
   })
   age_summary['PropSettled'] = age_summary['NoSet'] / age_summary['NoAlive']
   
   # Plot aggregated observations + continuous model curve
   plot(age_summary['LarvalAge'], age_summary['PropSettled'])  # Discrete points
   plot(age_range, predictions)  # Continuous curve
   ```

## 6. Key Distinctions

| Aspect | Discrete Observations | Continuous Model |
|--------|----------------------|------------------|
| **What** | Individual replicate proportions | Smooth function over age |
| **When** | Data collection | Model prediction |
| **Values** | One per replicate | Can evaluate at any age |
| **Example** | PropSettled = 0.8 at age 5.0 | competency(5.0) = 0.75, competency(5.1) = 0.755, ... |
| **Used for** | Model fitting (raw counts) | Visualization, prediction |

## 7. Mathematical Formulation

### Binomial Likelihood (Continuous Approach)

For each replicate $i$:
- Age: $t_i$
- Observed: $s_i$ settled out of $n_i$ total
- Model prediction: $p(t_i; \theta)$ (continuous function of age)

**Likelihood:**
$$L(\theta) = \prod_{i=1}^{N} \binom{n_i}{s_i} [p(t_i; \theta)]^{s_i} [1-p(t_i; \theta)]^{n_i-s_i}$$

**Log-likelihood:**
$$\ell(\theta) = \sum_{i=1}^{N} \left[ s_i \log(p(t_i; \theta)) + (n_i-s_i) \log(1-p(t_i; \theta)) \right]$$

**Key points:**
- Each replicate contributes individually
- Model $p(t; \theta)$ is continuous (can evaluate at any $t$)
- Sample sizes $n_i$ automatically weight the contributions
- No pre-aggregation needed

## Summary

**"Continuous proportions" means:**

1. ✅ **Model fitting**: Uses raw counts (binomial likelihood), not pre-aggregated proportions
2. ✅ **Model predictions**: Smooth, continuous curves that can be evaluated at any age
3. ✅ **Sample-size weighting**: Larger replicates automatically have more influence
4. ✅ **No information loss**: All replicate-level information is preserved

**The "continuous" part refers to:**
- The model function being continuous (can predict at any age)
- Not aggregating data before fitting (preserving all information)
- Generating smooth prediction curves for visualization

---

**Key Takeaway**: We fit models to **raw counts** (discrete observations), but the **model itself is continuous** (can predict at any age), and we generate **continuous prediction curves** for visualization.

