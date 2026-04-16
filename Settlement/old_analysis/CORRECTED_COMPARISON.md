# CORRECTED COMPARISON: Our Results vs. Randal et al. (2024)

## ❌ My Previous Error

I incorrectly stated that Merulinidae, Diploastreidae, and Euphylliidae prefer CCA over rubble.

## ✅ CORRECT Information from Table 2

**Settlement Cue Rankings** (from Randal's Table 2):

### Merulinidae (ALL prefer rubble most):
- *Dipsastraea matthaii*: **Rubble** > Disc = CCA
- *Dipsastraea pallida*: **Rubble** > CCA = Disc
- *Goniastrea retiformis*: **Rubble** = CCA > Disc
- *Mycedium elephantotus*: **Rubble** > CCA = Disc
- *Oulophyllia crispa*: **Rubble** > CCA
- *Platygyra daedalea*: **Rubble** > CCA > Disc

### Diploastraeidae:
- *Diploastrea heliopora*: **Rubble** > CCA > Disc

### Euphylliidae:
- *Galaxea fascicularis*: **Rubble** > CCA > Extract

### Conclusion from Table 2:
**RUBBLE is the most effective cue for essentially ALL families tested**, including those where our LD50 values are problematic.

---

## 🤔 So Why Are Our LD50 Values So Different?

### Comparison (Threshold = 0.3):

| Family | Our LD50 | Randal Mean TC50 | Difference | Issue |
|--------|----------|------------------|------------|-------|
| Acroporidae | 3.36 days | 4.33 days | -22% | ✅ Reasonable |
| Lobophylliidae | 2.54 days | 3.48 days | -27% | ✅ Reasonable |
| Poritidae | 3.01 days | 4.38 days | -31% | ⚠️ Moderate |
| **Diploastreidae** | **0.00 days** | **3.84 days** | **-100%** | ❌ **Major problem** |
| **Euphylliidae** | **0.00 days** | **2.92 days** | **-100%** | ❌ **Major problem** |
| **Merulinidae** | **0.00 days** | **2.53 days** | **-100%** | ❌ **Major problem** |

---

## 🔍 Possible Explanations for ~0 Day LD50 Values

Since rubble IS the right cue, the problem must be elsewhere:

### 1. **Binary Threshold Artifact**
   - Our method: If ANY observation exceeds 30% settlement, cumsum becomes 1
   - These families may have a few individuals settling very early (even at age 0)
   - This pulls the entire family LD50 to near-zero
   
### 2. **Data Aggregation Issue**
   - Pooling ALL species in a family together
   - If even ONE species/replicate shows early settlement, it affects the whole family
   - Merulinidae has 6 species - if any show early settlement, family average crashes

### 3. **Age=0 Initialization**
   - We added age=0 observations with NoSet=0
   - But cumulative calculation may be incorrectly triggering at age=0

### 4. **Model Fitting Problem**
   - Logistic regression can't handle cases where settlement happens immediately
   - Returns age~0 as the LD50
   - Need to check if model is converging properly

### 5. **Statistical vs Biological Interpretation**
   - Randal's TC50: "When does 50% of replicates exceed 30% settlement?"
   - Our LD50: "At what age is there 50% probability of the cumulative binary being 1?"
   - These may not be exactly equivalent

---

## 🔬 Let's Investigate Our Data

Looking at our output for these problematic families:

```
Diploastreidae
  Observations: 91
  Age range: 0-38 days
  Species in family: 1
  LD50: 0.0 days

Euphylliidae
  Observations: 122
  Age range: 0-36 days
  Species in family: 1
  LD50: 0.0 days

Merulinidae
  Observations: 719
  Age range: 0-77 days
  Species in family: 6
  LD50: 0.0 days
```

All show LD50 = 0.0 days, which means the model predicts immediate competency.

---

## 💡 Most Likely Explanation

**The cumulative binary approach is too sensitive to early settlers.**

### How Randal Handles This:
- Uses Bayesian hierarchical models with random effects
- Models individual replicates separately
- TC50 is calculated AFTER modeling settlement probability curves
- Considers plate-level variation

### How We Handle This:
- Aggregate all data to family level immediately
- Calculate cumulative binary (once ANY observation >30%, stays 1 forever)
- Fit simple logistic regression to binary outcome
- If a few individuals settle early, the entire family appears to have immediate competency

### Example Scenario:
```
Age 0: 0/100 settled (0%)     → Settlement = 0
Age 2: 5/100 settled (5%)     → Settlement = 0 (below 30%)
Age 3: 35/100 settled (35%)   → Settlement = 1 (exceeded 30%!)
Age 4: 40/100 settled (40%)   → Settlement = 1 (stays 1)
Age 5: 80/100 settled (80%)   → Settlement = 1
...

Cumsum kicks in at Age 3, so model fits:
  Ages 0-2: Settlement = 0
  Ages 3+: Settlement = 1
  
LD50 ≈ Age 1-2 (interpolation between 0 and 3)
```

But if there's ANY variability and a replicate hits 30% at age 0 or 1, the LD50 crashes to 0.

---

## 🛠️ How to Fix This

### Option 1: Use Higher Threshold
- Try threshold = 0.5 or 0.7 instead of 0.3
- More stringent criterion might delay cumsum trigger
- Check if LD50 values improve

### Option 2: Replicate-Level Analysis
- Don't aggregate immediately to family
- Calculate cumulative binary for each replicate separately
- Then model: "At what age do 50% of replicates show competency?"
- This is closer to Randal's approach

### Option 3: Use Continuous Proportion
- Instead of binary 0/1, use actual proportion settled
- Model proportion directly (like our first analysis)
- Calculate LD50 as age when proportion reaches 0.5

### Option 4: Filter Out Age 0-2
- Exclude very early ages where settlement shouldn't be possible
- Forces model to find LD50 in biologically reasonable range

---

## 📊 Quick Test: Check Other Thresholds

From our output:

**Threshold = 0.5:**
- Diploastreidae: 1.72 days (better!)
- Euphylliidae: 0.00 days (still bad)
- Merulinidae: 0.00 days (still bad)

**Threshold = 0.7:**
- Diploastreidae: 2.27 days (even better! Close to Randal's 3.84)
- Euphylliidae: 0.00 days (still bad)
- Merulinidae: 0.00 days (still bad)

### Conclusion:
Higher thresholds help Diploastreidae but not the other two families.

---

## 🎯 Final Interpretation

### What Went Right:
✅ Acroporidae (3.36 vs 4.33) - Good agreement  
✅ Lobophylliidae (2.54 vs 3.48) - Good agreement  
✅ Poritidae (3.01 vs 4.38) - Acceptable

### What Went Wrong:
❌ Diploastreidae, Euphylliidae, Merulinidae show immediate competency

### Why:
1. **Cumulative binary method is too sensitive** to early/sporadic settlement
2. **Family-level aggregation** loses replicate-level structure
3. **Need Randal's full Bayesian approach** with random effects to handle this properly

### The Real Lesson:
**Randal's Bayesian hierarchical approach is necessary for accurate LD50 estimation.**

Our simplified frequentist approach works reasonably well for:
- Families with consistent, gradual competency development (Acroporidae)
- Families with little early settlement noise

But fails for:
- Families with high variability
- Families with early settlers
- Situations requiring proper handling of hierarchical data structure

---

## ✅ Corrected Recommendation

**For connectivity modeling:**

1. **Use Randal's species-specific TC50 values directly** - they're the gold standard
2. **Our family-level approach has limited utility**:
   - OK for Acroporidae (within ~25%)
   - Problematic for most other families
3. **Don't use our LD50 values for Merulinidae, Diploastreidae, Euphylliidae** - they're artifactually low

**The comparison confirms that Randal's comprehensive Bayesian approach is necessary for robust estimates across diverse coral families.**

---

## 📚 Key Takeaway

I was wrong about the cue preference, but that wasn't the problem. The real issue is that:

1. ✅ Rubble IS the best cue for all families (you were right!)
2. ❌ Our simplified statistical approach can't handle the data structure properly
3. ✅ Bayesian hierarchical models (like Randal's) are needed for this type of data
4. ✅ Our Acroporidae estimate is decent, others are not reliable

**Thank you for the correction!** This highlights the importance of carefully reading the source material.

