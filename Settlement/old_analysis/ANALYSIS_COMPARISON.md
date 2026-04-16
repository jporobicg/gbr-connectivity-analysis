# Comparison: Randal's Approach vs. Current Analysis

## Randal's Methodology (from his GitHub repository)

### Data Processing Approach:

1. **Age Zero Initialization**:
   ```R
   # Adds age=0 observations with NoSet=0, NoNotSet=10, NoAlive=10
   # Forces settlement to be 0 at birth (biological constraint)
   ```

2. **Binary Settlement Scoring with Thresholds**:
   ```R
   # Tests multiple thresholds (0.1 to 0.9)
   mutate(Settle = as.numeric(NoSet/(NoSet + NoNotSet) > thresholdProp))
   ```

3. **Cumulative Settlement Tracking**:
   ```R
   # Once settlement exceeds threshold, it stays "competent"
   arrange(LarvalAge) %>%
   mutate(cumsumSettle = cumsum(Settle),
          Settlement = ifelse(cumsumSettle>0, 1, 0))
   ```

4. **Bayesian Logistic Regression**:
   ```R
   # Models binary Settlement outcome
   Settlement|trials(1) ~ SpecificTreatment * LarvalAge + (1|Plate)
   family = 'binomial'
   ```

5. **Treatment Focus**:
   - Excludes: control, peptide, extract treatments
   - Analyzes: CCA, rubble, disc (settlement cue treatments)
   - Models treatment effects explicitly

6. **Analysis by Species** (not Family):
   - Each species analyzed separately
   - Random effects for experimental Plate

---

## Current Analysis Methodology

### Key Features:

1. **Cumulative Maximum Proportion**:
   ```python
   # Direct calculation of cumulative competency
   agg_data['CumulativeCompetency'] = agg_data['PropSettled'].cummax()
   ```

2. **Continuous Competency Models**:
   - Logistic: `L / (1 + exp(-k * (age - x0)))`
   - Gompertz: `a * exp(-b * exp(-c * age))`
   - Weibull: Integral-based (following Monegetti)

3. **Family-Level Analysis**:
   - Groups multiple species by taxonomic family
   - Broader taxonomic patterns vs. species-specific

4. **CCA Treatment Only**:
   - Uses only CCA treatment data
   - Doesn't model treatment effects
   - Focuses on competency development with settlement cues

5. **Maximum Likelihood Estimation**:
   - Binomial likelihood: `NoSet * log(pred) + NoNotSet * log(1-pred)`
   - AIC for model selection

---

## Key Similarities

✓ **Cumulative Competency Concept**: Both track irreversible competency acquisition  
✓ **CCA Treatment**: Both recognize CCA as indicator of competency  
✓ **Age-Dependent Models**: Both model competency as function of larval age  
✓ **Binomial Framework**: Both account for binomial nature of settlement data

---

## Key Differences

| Aspect | Randal | Current Analysis |
|--------|--------|------------------|
| **Binary vs. Continuous** | Binary Settlement (0/1) | Continuous proportion (0-1) |
| **Taxonomic Level** | Species-specific | Family-level |
| **Statistical Framework** | Bayesian (brms/Stan) | Frequentist (MLE) |
| **Treatment Effects** | Models explicitly | Uses CCA only |
| **Random Effects** | Plate-level | None (aggregated) |
| **Threshold Sensitivity** | Tests 0.1-0.9 | Uses cummax (no threshold) |
| **Model Complexity** | Logistic regression + RE | Multiple functional forms |

---

## Strengths of Each Approach

### Randal's Strengths:
1. **Species-specific resolution** - captures fine-scale variation
2. **Treatment comparisons** - can compare CCA vs. rubble vs. disc
3. **Uncertainty quantification** - Bayesian credible intervals
4. **Random effects** - accounts for plate-to-plate variation
5. **Threshold sensitivity** - robust to scoring criteria

### Current Analysis Strengths:
1. **Taxonomic generality** - applicable across families
2. **Model flexibility** - tests multiple functional forms
3. **Comparison with Monegetti** - directly comparable curves
4. **Computational efficiency** - fast, simple models
5. **Visual clarity** - easy to interpret competency curves

---

## Recommendations for Integration

### Option 1: Hybrid Approach
- Use Randal's data processing (age=0, cumulative binary)
- Fit family-level models (like current analysis)
- Compare with both Randal's species results and Monegetti

### Option 2: Replicate Randal Exactly
- Fit species-specific Bayesian models
- Extract LD50 (age at 50% competency) for each species
- Aggregate to family level for comparison

### Option 3: Enhance Current Analysis
- Add threshold sensitivity analysis (like Randal)
- Include random effects if data structure permits
- Generate prediction intervals

---

## What the Current Analysis Provides That Randal's Doesn't

1. **Direct comparison with Monegetti model**:
   - Same functional forms (Weibull, exponential)
   - Overlaid competency curves
   - Explicit parameter comparisons

2. **Family-level patterns**:
   - Useful for taxa where species ID uncertain
   - Reveals phylogenetic patterns
   - More generalizable to connectivity models

3. **Multiple model forms**:
   - Logistic, Gompertz, Weibull all tested
   - AIC-based selection
   - Reveals which functional form fits best

4. **Simple implementation**:
   - Can be easily integrated into connectivity code
   - No Bayesian infrastructure needed
   - Fast computation

---

## Suggested Updates to Current Analysis

Based on Randal's approach, consider:

### 1. Add Age Zero Constraint
```python
# Force competency = 0 at age = 0
# Currently models may predict small positive values at age=0
```

### 2. Add Threshold Analysis
```python
# Test cumulative competency at different settlement thresholds
thresholds = [0.1, 0.3, 0.5, 0.7, 0.9]
for thresh in thresholds:
    # Calculate LD50 (age when 50% of cases exceed threshold)
```

### 3. Calculate LD50 Values
```python
# Age at which 50% of plates/replicates show competency
# Directly comparable to Randal's metrics
```

### 4. Add Species-Level Option
```python
# Allow analysis at species OR family level
# User can choose taxonomic resolution
```

---

## Implementation Priority

**High Priority:**
1. Add age=0 constraint to models ✓ (already somewhat addressed by tc parameter)
2. Calculate LD50/x0 values for each family ✓ (already in logistic model)
3. Document differences from Randal clearly ✓ (this document)

**Medium Priority:**
4. Add threshold sensitivity analysis
5. Compare LD50 values with Randal's results
6. Create species-level option

**Low Priority:**
7. Implement Bayesian versions
8. Add random effects structure
9. Full replication of Randal's analysis

---

## Conclusion

**Current Analysis Status**: ✓ **VALID AND COMPLEMENTARY**

The current analysis provides:
- Family-level competency models suitable for connectivity applications
- Direct comparison with Monegetti's framework
- Multiple functional forms with model selection
- Clear, interpretable competency-by-age curves

It differs from Randal's approach in:
- Taxonomic resolution (family vs. species)
- Statistical framework (frequentist vs. Bayesian)
- Not modeling treatment effects explicitly

**Both approaches are scientifically valid** and address slightly different questions:
- **Randal**: "How do different settlement cues affect species-specific competency timing?"
- **Current**: "What are family-level competency curves comparable to Monegetti's model?"

**For connectivity modeling**, the current analysis provides **ready-to-use competency functions** at the family level, which is often the appropriate scale for large-scale dispersal models where species-level detail may not be available or necessary.

