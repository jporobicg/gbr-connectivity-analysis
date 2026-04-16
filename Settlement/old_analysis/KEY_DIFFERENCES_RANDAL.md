# Key Differences Between Our Approach and Randal's

After examining Randal's R code, here are the critical differences that explain why our LD50 values differ:

---

## 1. **Data Processing: SPECIES vs FAMILY**

### Randal's Approach (lines 27-42 in 20_fitModels.R):
```r
data.q1 <- data %>%
    filter(Q1 == 'Y') %>%
    # ... filtering ...
    mutate(FractionSettled = as.numeric(NoSet/(NoSet + NoNotSet))) %>%
    mutate(Settle = as.numeric(NoSet/(NoSet + NoNotSet) > thresholdProp)) %>% 
    group_by(Species, SpecificTreatment) %>%  # ← GROUPS BY SPECIES
    arrange(LarvalAge) %>%
    mutate(cumsumSettle = cumsum(Settle),
           Settlement = ifelse(cumsumSettle>0, 1, 0)) %>%
    ungroup() %>%
    group_by(Species) %>%  # ← NESTS BY SPECIES
    nest()
```

### Our Approach:
```python
# We group by FAMILY, not species
for family, group in df_treatment.groupby(['Family']):
    group_processed = calculate_cumulative_settlement(group, threshold)
```

**Impact:** 
- Randal keeps species separate throughout the entire analysis
- We aggregate all species within a family immediately
- This loses species-level structure and plate-level replication

---

## 2. **Statistical Model: BAYESIAN HIERARCHICAL vs SIMPLE GLM**

### Randal's Model (lines 82-108 in functions.R):
```r
fitModel <- function(dat) {
    # Bayesian hierarchical model with:
    form <- bf(Settlement|trials(1) ~ SpecificTreatment * LarvalAge + (1|Plate),
               family = 'binomial')
    mod <- brm(form,
               prior = priors,
               data = dat,
               iter = 6000,    # 6000 MCMC iterations
               warmup = 2000,  # 2000 burn-in
               chains = 3,     # 3 MCMC chains
               cores = 3,
               thin = 10,
               control = list(adapt_delta = 0.99, max_treedepth = 20),
               backend = 'cmdstanr')  # Uses Stan (full Bayesian inference)
}
```

**Key features:**
- Treatment × Age interaction
- **Random effects for Plate**: `(1|Plate)`
- Bayesian priors on all parameters
- MCMC sampling (18,000 total draws, thinned to 1,200)
- Full uncertainty propagation

### Our Model:
```python
def fit_logistic_model(df_family):
    X = df_family[['LarvalAge']].values
    y = df_family['Settlement'].values
    
    model = LogisticRegression(max_iter=1000, random_state=42)
    model.fit(X_scaled, y)  # Simple frequentist GLM
```

**Key features:**
- No treatment effects (we only use rubble)
- **No random effects**
- No hierarchical structure
- Maximum likelihood point estimates
- No uncertainty propagation

---

## 3. **LD50 Calculation: FROM POSTERIOR vs INVERSE LOGIT**

### Randal's LD50 Calculation (lines 481-520 in functions.R):
```r
LD50 <- function(mod) {
    # Extract intercept and slope for each treatment from posterior
    e1 <- mod %>% emmeans(~SpecificTreatment|LarvalAge, at = list(LarvalAge=0)) %>%
        tidy_draws()  # Get posterior draws of intercept
    
    e2 <- mod %>% emtrends(specs = 'SpecificTreatment', var = 'LarvalAge') %>%
        tidy_draws()  # Get posterior draws of slope
    
    # Join and calculate LD50 for each posterior draw
    ee %>% 
        pivot_wider() %>%
        mutate(LD50 = -(intercept / slope) %>% inv.logit())  # LD50 calculation
    
    # Returns full posterior distribution of LD50
}
```

**Process:**
1. Extract intercept (α) and slope (β) from Bayesian posterior
2. Calculate for EACH posterior draw: LD50 = -α/β transformed through inverse logit
3. Returns distribution of LD50 with uncertainty
4. Takes median of posterior as TC50 estimate

### Our LD50 Calculation:
```python
def calculate_ld50(model_info):
    model = model_info['model']
    scaler = model_info['scaler']
    
    def prob_at_age(age):
        age_scaled = scaler.transform([[age]])
        prob = model.predict_proba(age_scaled)[0, 1]
        return abs(prob - 0.5)  # Find where probability = 0.5
    
    # Numerical optimization to find age where P(Settlement=1) = 0.5
    result = minimize_scalar(prob_at_age, bounds=(min_age, max_age))
    return result.x  # Returns single point estimate
```

**Process:**
1. Get single point estimate of model parameters
2. Numerically search for age where predicted probability = 0.5
3. Returns single value with no uncertainty

---

## 4. **Random Effects / Plate Structure**

### Randal:
```r
Settlement ~ SpecificTreatment * LarvalAge + (1|Plate)
```
- Each experimental plate gets its own random intercept
- Accounts for plate-to-plate variation
- Partial pooling across plates
- More robust to outliers

### Ours:
- No plate structure retained
- All observations pooled
- Single fixed effects model
- More sensitive to outliers/early settlers

---

## 5. **Treatment Handling**

### Randal:
```r
# Models ALL treatments simultaneously
Settlement ~ SpecificTreatment * LarvalAge

# Then extracts LD50 for EACH treatment
LD50_CCA, LD50_rubble, LD50_disc, etc.

# Reports the EARLIEST TC50 across treatments
```

### Ours:
```python
# Filter to rubble ONLY
df_treatment = df[df['SpecificTreatment'] == treatment]

# Model only rubble
Settlement ~ LarvalAge  # No treatment effects
```

**Impact:** Randal's approach allows comparison across treatments and uses the best one for TC50.

---

## 6. **Why Our Values Are Lower (Especially the 0s)**

### Problem: Immediate Competency Artifact

When we aggregate to family level and use cumulative binary:

**Scenario:**
- Merulinidae has 6 species
- Species A: First settles at age 2
- Species B: First settles at age 3  
- Species C: First settles at age 0 (maybe spurious/outlier)
- Species D-F: First settle at age 2-3

**Randal's handling:**
- Fits separate model for each species
- Random effect for each plate
- Species C's age-0 settlement is moderated by:
  - Plate random effect (maybe just that plate was weird)
  - Bayesian priors (pulls extreme values toward reasonable)
  - Separate species model (doesn't affect others)
- **Result:** Species C gets TC50 ~ 2 days, others get their own values
- **Family average:** ~2.5 days

**Our handling:**
- All species pooled immediately
- No plate structure
- No hierarchical borrowing of strength
- Species C's age-0 settlement counts equally with others
- Cumsum triggers at age 0 for the pooled data
- Logistic regression fits: "Settlement = 1 starting at age ~0"
- **Result:** Family LD50 = 0.0 days

---

## 7. **Data Filtering Difference**

### Randal (line 28-29):
```r
filter(Q1 == 'Y') %>%  # Only include designated observations
filter(!(NoSet ==0 & NoNotSet ==0))  # Remove empty observations
```

They have a quality control flag (Q1) that we don't use.

### Ours:
```python
# We use all data, no Q1 filter
# May include lower-quality observations
```

---

## 8. **Summary of Why Randal's Approach Is Better**

| Feature | Randal | Us | Impact |
|---------|--------|-----|--------|
| **Taxonomic Level** | Species | Family | Loses species variation |
| **Statistical Framework** | Bayesian MCMC | Frequentist MLE | No uncertainty, no pooling |
| **Random Effects** | Yes (Plate) | No | Can't handle hierarchical data |
| **Treatment Modeling** | All simultaneously | One at a time | Can't compare, can't use best |
| **Uncertainty** | Full posterior | Point estimate | No confidence intervals |
| **Outlier Robustness** | High (priors + RE) | Low | Sensitive to early settlers |
| **Computational Cost** | High (hours) | Low (seconds) | Trade-off |

---

## 9. **The Real Issue with Our 0-Day Values**

Looking at Randal's code, the cumulative binary calculation is **identical**:

```r
# Randal (line 37-38)
mutate(cumsumSettle = cumsum(Settle),
       Settlement = ifelse(cumsumSettle>0, 1, 0))
```

```python
# Ours
df['cumsumSettle'] = df['Settle'].cumsum()
df['Settlement'] = (df['cumsumSettle'] > 0).astype(int)
```

**So the data processing is the same.**

The difference is what happens AFTER:

### Randal:
1. Creates binary Settlement variable (cumulative)
2. **Fits Bayesian hierarchical model to binary outcome**
3. **Random effects absorb plate-level variation**
4. **Treatment effects model different cues**
5. Extracts LD50 from smooth posterior probability curves
6. Takes median across posterior draws

### Ours:
1. Creates binary Settlement variable (cumulative) 
2. **Fits simple logistic regression to binary outcome**
3. **No random effects** - all variation attributed to age
4. **No treatment effects** - single cue only
5. Finds age where predicted probability = 0.5
6. Single point estimate

---

## 10. **Specific Code Causing Our Problem**

In our family-level aggregation:

```python
# This is the problem:
for family, group in df_treatment.groupby(['Family']):
    group_processed = calculate_cumulative_settlement(group, threshold)
    # ↑ Cumsum is calculated across ALL observations in family
    # If ANY observation hits threshold early, cumsum = 1 for that age onward
    # Loses track of which species/plate it came from
```

Should be:

```python
# Better approach (like Randal):
for species in df_treatment['Species'].unique():
    df_species = df_treatment[df_treatment['Species'] == species]
    for plate in df_species['Plate'].unique():
        df_plate = df_species[df_species['Plate'] == plate]
        # Calculate cumsum at plate level
        # Then model with random effects
```

---

## 11. **Conclusion**

Our simplified approach fails for families with:
- High inter-species variation
- Early/sporadic settlers
- Small sample sizes per species/plate

It works reasonably for:
- Acroporidae (large N, consistent behavior)
- Families with uniform competency patterns

**To match Randal's results, we would need:**
1. ✅ Species-level analysis (not family)
2. ✅ Bayesian hierarchical models with random effects
3. ✅ Multiple treatment modeling
4. ✅ Posterior distribution of LD50
5. ✅ Quality control filtering (Q1 flag)

**Or better yet:** Just use Randal's published values! They did it right.

---

## 12. **Why Acroporidae Works Better**

Acroporidae shows better agreement (3.36 vs 4.33, 22% diff) because:
- Large sample size (12 species, 934 observations)
- More consistent competency timing across species
- Law of large numbers smooths out outliers
- Family-level aggregation less problematic

But even for Acroporidae, we're still underestimating by ~1 day (22%) because:
- No treatment comparison (some species prefer CCA/disc)
- No random effects (plate variation)
- No hierarchical structure (species variation)

---

**Final Answer:** The differences aren't about the cue (rubble is correct), but about **statistical methodology**. Randal's Bayesian hierarchical approach properly handles the complex data structure, while our simple frequentist approach cannot.

