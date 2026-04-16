# Why Do Settlement Values Reach 100% in Randal's Data?

## The Question

**Observation**: In Randal's data, some rows show `PropSettled = 1.0` (100% settlement).  
**Question**: Why does this happen, and is it a problem?

## The Answer: Replicate-Level vs. Aggregated Data

### Key Finding

**100% settlement occurs at the REPLICATE LEVEL, not when aggregated.**

### Evidence

1. **Replicate-Level Data** (raw observations):
   - **247 rows** with 100% settlement in Merulinidae
   - **527 rows** with 100% settlement in Acroporidae
   - Typical sample size: **~10 larvae per replicate**
   - When all 10 larvae in a replicate settle → 100% for that replicate

2. **Aggregated Data** (by age + treatment):
   - **0 rows** with 100% settlement
   - When aggregating across multiple replicates, values are lower
   - Averaging across replicates eliminates perfect settlement

### Why This Happens

#### 1. **Small Sample Sizes Per Replicate**

Each replicate typically has **~10 larvae**:
- If all 10 settle → 100% for that replicate
- This is **statistically valid** but represents a **single experimental unit**
- With small N, 100% is possible even if true probability < 100%

#### 2. **Experimental Design**

Randal's experimental design:
- Multiple replicates per age/treatment combination
- Each replicate is an independent experimental unit
- Replicates are pooled for analysis

**Example**:
- Replicate 1 (age 5, rubble): 10/10 settled = 100%
- Replicate 2 (age 5, rubble): 8/10 settled = 80%
- Replicate 3 (age 5, rubble): 9/10 settled = 90%
- **Aggregated**: 27/30 = 90% (not 100%)

#### 3. **Treatment Effects**

100% settlement is more common in:
- **Rubble treatment**: 176 cases (best treatment)
- **CCA treatment**: 63 cases
- **Disc, peptide**: Fewer cases

This makes biological sense: better treatments → higher settlement → more cases of 100%

### Is This a Problem?

**No, this is normal and expected:**

1. ✅ **Replicate-level variation**: Individual replicates can show 100%
2. ✅ **Aggregation removes extremes**: When aggregated, values are more moderate
3. ✅ **Model fits to aggregated data**: The model doesn't see these 100% values directly
4. ✅ **Statistically valid**: With small N per replicate, 100% is possible

### What the Model Sees

The model fits to **aggregated data** (by age):
- Data is grouped by `LarvalAge`
- Multiple replicates are summed: `NoSet` and `NoAlive` are summed across replicates
- Result: **No 100% values** in aggregated data
- Model predictions: **0.2-0.4** (20-40%), which is correct

### Comparison: Replicate vs. Aggregated

**Merulinidae Example**:

| Level | 100% Cases | Max Settlement | Interpretation |
|-------|------------|---------------|---------------|
| **Replicate** | 247 rows | 100% | Individual replicates can show perfect settlement |
| **Aggregated** | 0 rows | ~30% | Aggregated across replicates shows realistic values |

### Biological Interpretation

**Why 100% at replicate level is biologically meaningful:**

1. **Optimal conditions**: Some replicates have perfect conditions
   - All larvae are competent
   - Cues are optimal
   - Environmental conditions are ideal
   - Result: 100% settlement in that replicate

2. **Small sample sizes**: With N=10, 100% is possible even if true probability is 80-90%
   - Statistical variation allows perfect outcomes
   - Not necessarily indicating 100% true competency

3. **Treatment quality**: Better treatments (rubble, CCA) show more 100% cases
   - Reflects treatment effectiveness
   - Biologically meaningful pattern

### Why the Model Doesn't Predict 100%

The model fits to **aggregated data**, which:
- Averages across replicates
- Includes variation between replicates
- Shows realistic settlement rates (20-40%)
- Doesn't capture replicate-level extremes

**This is correct behavior:**
- Models should fit to aggregated data for robustness
- Replicate-level extremes are averaged out
- Predictions represent expected values, not extremes

### Data Structure

```
Raw Data (Replicate Level):
- Age 5, Rubble, Replicate 1: 10/10 = 100% ← Individual replicate
- Age 5, Rubble, Replicate 2: 8/10 = 80%
- Age 5, Rubble, Replicate 3: 9/10 = 90%

Aggregated Data (Model Input):
- Age 5, Rubble: 27/30 = 90% ← Summed across replicates
```

### Conclusion

**100% settlement values are:**
1. ✅ **Normal**: Occur at replicate level with small sample sizes
2. ✅ **Biologically meaningful**: Reflect optimal conditions in some replicates
3. ✅ **Not a problem**: Aggregation removes them before model fitting
4. ✅ **Expected**: Better treatments show more 100% cases

**The model correctly:**
- Fits to aggregated data (no 100% values)
- Predicts realistic settlement rates (20-40%)
- Accounts for variation between replicates
- Produces biologically meaningful results

---

**Key Takeaway**: 100% values at the replicate level are a feature of the experimental design (small N per replicate), not a data quality issue. Aggregation before model fitting correctly handles this, producing realistic predictions.

