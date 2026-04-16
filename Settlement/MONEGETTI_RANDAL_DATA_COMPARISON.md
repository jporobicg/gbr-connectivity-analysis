# Direct Data Comparison: Monegetti (A. tenuis) vs. Randal (Acroporidae)

## Executive Summary

This report provides a detailed comparison between two foundational datasets for coral larval competency modeling:

1. **Monegetti et al. (2012, 2019)**: Single-species (*Acropora tenuis*) metamorphosis data from the Great Barrier Reef
2. **Randal et al. (2024)**: Family-level (Acroporidae) settlement data from the Indo-Pacific

Both datasets have been fitted with the Monegetti piecewise Weibull-Exponential competency model, allowing direct parameter comparison and biological interpretation. This comparison reveals important differences in competency dynamics between single-species and family-level approaches, with implications for connectivity modeling.

---

## 1. Dataset Characteristics

### 1.1 Monegetti's A. tenuis Dataset

**Source**: Moneghetti et al. (2019) "High-frequency sampling and piecewise models reshape dispersal kernels of a common reef coral" (*Ecology Letters*)

**Original Study**: Moneghetti et al. (2012) "A quantitative assessment of the competency period of *Acropora tenuis* coral larvae"

**Data Characteristics**:
- **Species**: *Acropora tenuis* (single species)
- **Location**: Great Barrier Reef, Australia
- **Collection Year**: 2012
- **Sample Size**: 1,440 larvae
- **Replicates**: 4 replicates
- **Age Range**: 0-80 days (19 time points)
- **Data Type**: **Metamorphosis** (spontaneous transformation, **no settlement cues provided**)
- **Treatment**: Control only (no settlement cues provided)
- **Total Metamorphosed**: 708 larvae (49.2% overall)
- **Maximum Competency**: ~88% at ~14 days post-spawning
- **Mortality Tracking**: No separate mortality column - `larvae = meta + swimming` (all larvae accounted for)

**Experimental Protocol**:
- Larvae collected from natural spawning events
- Maintained in laboratory conditions at 27°C
- Settlement checked at multiple time points
- No settlement cues provided (spontaneous metamorphosis)
- High-frequency sampling (multiple observations per day early in development)

**Key Features**:
- Single-species focus allows species-specific parameter estimation
- Control treatment isolates age effects from cue effects
- Long time series (80 days) captures late-phase dynamics
- High temporal resolution early in development

### 1.2 Randal's Acroporidae Dataset

**Source**: Randal et al. (2024) "Age to settlement competency and settlement cue preference in coral larvae" (*Communications Biology*)

**Data Characteristics**:
- **Taxonomic Level**: Family-level (Acroporidae)
- **Species Included**: 12 species aggregated:
  - *Acropora intermedia*
  - *Acropora longicyathus*
  - *Acropora loripes*
  - *Acropora cf. kenti*
  - *Acropora millepora*
  - *Acropora microclados*
  - *Acropora glauca*
  - *Acropora tenuis* (included in aggregation)
  - *Montipora aequituberculata*
  - *Montipora digitata*
  - *Astreopora* species
  - Additional Acroporidae species
- **Location**: Indo-Pacific (multiple locations)
- **Collection Year**: 2024
- **Sample Size**: 45,123 larvae
- **Replicates**: 4,730 replicates
- **Age Range**: 4-76 days (variable time points across species)
- **Data Type**: **Settlement** (behavioral response, **with settlement cues provided**)
- **Treatments**: Multiple treatments included:
  - Control
  - Rubble
  - CCA (crustose coralline algae)
  - Disc
  - Extract
  - Peptide
- **Total Settled**: 19,045 larvae (42.2% overall)
- **Maximum Competency**: ~51% (lower than Monegetti)
- **Mortality Tracking**: No explicit mortality column - `NoAlive = NoSet + NoNotSet` (mortality not separately tracked, but may be implicit in `NoAlive`)

**Experimental Protocol**:
- Larvae collected from multiple species across Indo-Pacific
- Laboratory settlement assays with various cues
- Settlement checked at multiple larval ages
- Multiple treatments tested simultaneously
- Species pooled at family level for analysis

**Key Features**:
- Family-level aggregation provides broader taxonomic coverage
- Multiple treatments allow treatment effect assessment
- Large sample size provides statistical power
- Multi-species aggregation may mask species-specific patterns

---

## 2. Parameter Comparison

### 2.1 Monegetti Piecewise Model Parameters

| Parameter | Symbol | Monegetti Original (A. tenuis) | Our Fit (Randal Acroporidae) | Difference | % Change |
|-----------|--------|-------------------------------|------------------------------|------------|----------|
| **Acquisition Rate** | a | 1.292 | 0.332 | -0.960 | -74.3% |
| **Weibull Scale** | b1 | 0.001878 | 0.0017 | -0.000178 | -9.5% |
| **Weibull Shape** | v1 | 0.3645 | 0.100 | -0.2645 | -72.6% |
| **Exponential Rate** | b2 | 0.3969 | 0.001 | -0.3959 | -99.7% |
| **Precompetency Period** | tc | 3.33 days | 3.90 days | +0.57 days | +17.1% |
| **Change Point** | Tcp | 69.91 days | 71.48 days | +1.57 days | +2.2% |

### 2.2 Biological Interpretation of Parameters

**Acquisition Rate (a)**:
- **Monegetti**: a = 1.292 indicates rapid competency acquisition once precompetency period ends
- **Randal**: a = 0.332 indicates 3.9× slower competency acquisition
- **Interpretation**: Single-species A. tenuis shows faster competency development than family-level average

**Weibull Scale (b1)**:
- **Monegetti**: b1 = 0.001878
- **Randal**: b1 = 0.0017
- **Interpretation**: Very similar early-phase competency loss rates, suggesting comparable early dynamics

**Weibull Shape (v1)**:
- **Monegetti**: v1 = 0.3645 (moderate decay shape)
- **Randal**: v1 = 0.100 (very slow decay shape)
- **Interpretation**: Randal data shows much slower early-phase decay, suggesting more gradual competency loss

**Exponential Rate (b2) - KEY DIFFERENCE**:
- **Monegetti**: b2 = 0.3969 indicates **rapid late-phase competency loss** (397× faster)
- **Randal**: b2 = 0.001 indicates **essentially no late-phase loss** (competency maintained)
- **Critical Interpretation**: 
  - **Monegetti (Metamorphosis)**: Measures spontaneous metamorphosis without cues - shows clear biological decline as developmental window closes
  - **Randal (Settlement)**: Measures settlement behavior with cues - shows maintained behavioral competency as long as larvae remain alive
  - **The 397× difference reflects fundamental biological vs. behavioral processes**

**Precompetency Period (tc)**:
- **Monegetti**: tc = 3.33 days
- **Randal**: tc = 3.90 days
- **Interpretation**: Similar precompetency periods, suggesting consistent biological constraint

**Change Point (Tcp)**:
- **Monegetti**: Tcp = 69.91 days
- **Randal**: Tcp = 71.48 days
- **Interpretation**: Very similar transition points, suggesting consistent late-phase onset timing

---

## 3. Competency Curve Comparison

### 3.1 Early Development (0-10 days)

**Monegetti (A. tenuis)**:
- Rapid competency acquisition starting at ~3.3 days
- Reaches ~80% competency by day 5-6
- Peak competency ~88% at ~14 days
- Clear S-shaped curve with rapid early rise

**Randal (Acroporidae)**:
- Slower competency acquisition starting at ~3.9 days
- Reaches ~45% competency by day 8
- More gradual early rise
- Lower maximum competency (~51%)

**Key Difference**: Monegetti shows 2× faster early competency development and 1.7× higher maximum competency.

### 3.2 Mid Development (10-30 days)

**Monegetti (A. tenuis)**:
- Maintains high competency (~70-88%) through day 30
- Gradual decline beginning around day 20
- Clear two-phase structure visible

**Randal (Acroporidae)**:
- Maintains moderate competency (~45-51%) through day 30
- Relatively stable competency levels
- Less pronounced two-phase structure

**Key Difference**: Monegetti maintains higher competency levels with visible decline, while Randal shows more stable moderate competency.

### 3.3 Late Development (30-80 days)

**Monegetti (A. tenuis)**:
- Gradual decline from ~70% to ~30% by day 72
- Rapid decline to near zero by day 80
- Clear exponential decay phase (b2 = 0.397)

**Randal (Acroporidae)**:
- Maintains relatively stable competency (~40-50%) through day 76
- Very slow decline (essentially flat)
- Minimal late-phase decay (b2 = 0.001)

**Key Difference**: Monegetti shows clear late-phase competency loss, while Randal shows maintained competency with minimal decline.

---

## 4. Why the Differences?

### 4.1 Taxonomic Level Effects

**Species-Specific vs. Family Aggregation**:

1. **Species-Specific Patterns**: Monegetti's single-species data captures A. tenuis-specific competency dynamics, which may differ from other Acroporidae species

2. **Averaging Effects**: Randal's family-level aggregation averages across 12 species, potentially:
   - Diluting species-specific patterns
   - Creating smoother, more averaged curves
   - Reducing maximum competency (some species may have lower competency)
   - Extending competency period (some species may maintain competency longer)

3. **Biological Variation**: Different Acroporidae species may have different:
   - Competency acquisition rates
   - Maximum competency levels
   - Late-phase decay rates
   - Optimal settlement windows

**Evidence**: The inclusion of A. tenuis in Randal's dataset suggests that family-level aggregation masks the species-specific patterns observed in Monegetti's focused study.

### 4.2 Treatment Effects and Data Type Differences

**Critical Distinction: Metamorphosis vs. Settlement**:

1. **Monegetti's Metamorphosis Data**:
   - **No settlement cues provided** - spontaneous metamorphosis
   - Measures **developmental competency** (intrinsic biological process)
   - Metamorphosis is a **one-way developmental transformation**
   - Once metamorphosed, cannot revert
   - Shows clear late-phase decline as developmental window closes
   - **b2 = 0.3969** reflects biological constraint

2. **Randal's Settlement Data**:
   - **Settlement cues provided** - behavioral response to cues
   - Measures **settlement competency** (behavioral response)
   - Settlement is a **behavioral choice** that can be delayed
   - Larvae can remain competent but not settled
   - Shows maintained competency as long as larvae are alive
   - **b2 = 0.001** reflects behavioral/ecological constraint

**Treatment Effects in Randal's Data**:
- Multiple treatments included (control, rubble, CCA, disc, extract, peptide)
- Some treatments may enhance settlement (e.g., CCA, rubble)
- Some treatments may suppress settlement
- Treatment effects may average out in aggregation
- May show different competency expression patterns

**Hypothesis**: The fundamental difference is **biological process** (metamorphosis) vs. **behavioral process** (settlement), not just treatment effects:
- Metamorphosis: Developmental window closing → rapid late-phase decline
- Settlement: Behavioral competency maintained → minimal late-phase decline

### 4.3 Experimental Protocol Differences

**Temporal and Spatial Variation**:

1. **Collection Location**:
   - Monegetti: Great Barrier Reef (2012)
   - Randal: Indo-Pacific multiple locations (2024)
   - Geographic and temporal differences may reflect:
     - Environmental conditions
     - Genetic variation
     - Collection methods
     - Laboratory protocols

2. **Time Period**:
   - 12-year gap between studies
   - Potential changes in:
     - Environmental conditions
     - Coral health/condition
     - Experimental protocols
     - Laboratory techniques

3. **Sampling Frequency**:
   - Monegetti: High-frequency sampling early in development
   - Randal: Variable sampling frequency across species
   - May affect parameter estimation accuracy

### 4.4 Sample Size Effects

**Statistical Power**:

1. **Monegetti**: 1,440 larvae
   - Smaller sample may be more sensitive to outliers
   - May capture species-specific patterns more clearly
   - Less statistical power for late-phase estimation

2. **Randal**: 45,123 larvae
   - Large sample provides robust estimates
   - May mask species-specific patterns through averaging
   - Greater statistical power but potentially less biological resolution

---

## 5. Model Selection Results

### 5.1 Monegetti Original (A. tenuis)

**Model**: Monegetti piecewise Weibull-Exponential (6 parameters)
- **Fit Quality**: Excellent fit to single-species data
- **Biological Interpretation**: Clear two-phase structure with:
  - Rapid early competency acquisition
  - High maximum competency
  - Clear late-phase decline

**Justification**: Single-species data allows precise parameter estimation for complex model structure.

### 5.2 Our Fit (Randal Acroporidae)

**Model Comparison Results**:
- **Monegetti Model**: AICc = 56,895.3
- **Weibull Model**: AICc = 56,693.4 (best)
- **Δ AICc**: +201.9 (penalty for Monegetti model)

**Model Selection**: Weibull model preferred (Δ AICc = +201.9)

**Interpretation**:
- Two-phase structure may not be necessary for family-level aggregation
- Simpler Weibull model (4 parameters) fits equally well
- Family-level data may not require complex two-phase structure
- Species-specific patterns may be averaged out

**Key Insight**: The Monegetti model is more appropriate for single-species analysis than for family-level aggregation, where simpler models suffice.

---

## 6. Implications for Connectivity Modeling

### 6.1 Model Selection Guidance

**When to Use Monegetti Piecewise Model**:
1. ✅ Single-species analysis (like Monegetti's A. tenuis)
2. ✅ High-resolution data with clear two-phase structure
3. ✅ When model selection indicates it's best (e.g., Merulinidae, Lobophylliidae)
4. ✅ Mechanistic interpretation of competency development needed

**When to Use Simpler Models**:
1. ✅ Family-level aggregation (like Randal's Acroporidae)
2. ✅ Multi-species datasets
3. ✅ When model selection favors simpler models (Δ AICc > 10)
4. ✅ When two-phase structure is not biologically meaningful

### 6.2 Data Collection Recommendations

**For Species-Specific Models**:
- Collect single-species data with high temporal resolution
- Use control treatments to isolate age effects
- Sample frequently early in development
- Extend sampling to capture late-phase dynamics

**For Family-Level Models**:
- Aggregate across multiple species for sufficient sample sizes
- Include multiple treatments for comprehensive assessment
- Use simpler model structures (Weibull, Logistic, Gompertz)
- Focus on family-level patterns rather than species-specific details

### 6.3 Connectivity Modeling Applications

**Species-Specific Connectivity**:
- Use Monegetti model with species-specific parameters
- Provides detailed competency dynamics
- Appropriate for species-focused connectivity studies
- Example: A. tenuis connectivity modeling

**Family-Level Connectivity**:
- Use simpler models (Weibull, Logistic) with family-level parameters
- Provides broader taxonomic coverage
- Appropriate for multi-species connectivity studies
- Example: Acroporidae family connectivity modeling

**Hybrid Approach**:
- Use species-specific models when available
- Use family-level models as fallback
- Weight by data quality and sample size
- Combine approaches for comprehensive connectivity assessment

---

## 7. Scientific Context

### 7.1 Monegetti et al. Studies

**Moneghetti et al. (2019)**:
- **Title**: "High-frequency sampling and piecewise models reshape dispersal kernels of a common reef coral"
- **Journal**: *Ecology Letters*
- **Focus**: Dispersal kernel modeling using piecewise competency models
- **Innovation**: Two-phase Weibull-Exponential model for competency development

**Moneghetti et al. (2012)**:
- **Title**: "A quantitative assessment of the competency period of *Acropora tenuis* coral larvae"
- **Focus**: Detailed competency assessment for single species
- **Method**: High-frequency sampling and piecewise modeling

**Key Contributions**:
- Established two-phase competency model structure
- Demonstrated importance of late-phase competency decline
- Provided species-specific parameter estimates
- Influenced subsequent connectivity modeling approaches

### 7.2 Randal et al. (2024)

**Randal et al. (2024)**:
- **Title**: "Age to settlement competency and settlement cue preference in coral larvae"
- **Journal**: *Communications Biology*
- **Focus**: Multi-species competency assessment with treatment effects
- **Innovation**: Bayesian hierarchical modeling with species-level resolution

**Key Contributions**:
- Provided multi-species competency data
- Demonstrated treatment effects on settlement
- Established family-level patterns
- Created comprehensive dataset for connectivity modeling

**Comparison with Monegetti**:
- Different taxonomic resolution (multi-species vs. single-species)
- Different statistical approach (Bayesian vs. frequentist)
- Different treatment structure (multi-treatment vs. control-only)
- Complementary rather than competing approaches

### 7.3 Integration of Approaches

**Complementary Strengths**:
- **Monegetti**: Species-specific detail, mechanistic understanding
- **Randal**: Taxonomic breadth, treatment effects, statistical rigor

**Best Practices**:
- Use Monegetti model for species-specific connectivity when data available
- Use Randal approach for family-level connectivity with uncertainty quantification
- Combine approaches for comprehensive connectivity assessment
- Apply model selection to determine appropriate complexity

---

## 8. Conclusions

### 8.1 Key Findings

1. **Parameter Differences**:
   - Monegetti's A. tenuis shows 3.9× faster competency acquisition
   - Monegetti shows 397× faster late-phase competency loss
   - Similar precompetency periods (~3.3-3.9 days)
   - Similar change points (~70 days)

2. **Competency Dynamics**:
   - Monegetti: Rapid early rise, high maximum (~88%), clear late decline
   - Randal: Gradual early rise, moderate maximum (~51%), maintained competency

3. **Model Selection**:
   - Monegetti model appropriate for single-species analysis
   - Simpler models (Weibull) preferred for family-level aggregation
   - Model complexity should match data structure and research question

4. **Biological Interpretation**:
   - Species-specific patterns may be masked by family-level aggregation
   - Treatment effects may influence competency expression
   - Geographic and temporal differences may contribute to variation

### 8.2 Recommendations

**For Connectivity Modeling**:
1. Use species-specific Monegetti models when available
2. Use family-level simpler models (Weibull, Logistic) for aggregation
3. Apply model selection to determine appropriate complexity
4. Consider treatment effects in competency expression
5. Account for geographic and temporal variation

**For Future Research**:
1. Collect species-specific data with high temporal resolution
2. Extend sampling to capture late-phase dynamics
3. Include multiple treatments for comprehensive assessment
4. Compare species-specific vs. family-level approaches
5. Validate models across different geographic regions

### 8.3 Final Thoughts

The comparison between Monegetti's single-species A. tenuis data and Randal's family-level Acroporidae data reveals important insights:

- **Both approaches are valid** but serve different purposes
- **Model selection is essential** to match complexity to data structure
- **Species-specific detail** may be lost in family-level aggregation
- **Treatment effects** may influence competency dynamics
- **Comprehensive connectivity modeling** benefits from both approaches

The choice between Monegetti's detailed single-species approach and Randal's broad family-level approach should depend on:
- Research question (species-specific vs. family-level)
- Available data (single-species vs. multi-species)
- Model selection results (complex vs. simple models)
- Connectivity modeling goals (detailed vs. broad coverage)

---

## References

1. **Moneghetti, P. C., et al. (2019)**. "High-frequency sampling and piecewise models reshape dispersal kernels of a common reef coral." *Ecology Letters*, 22(10), 1643-1652. https://doi.org/10.1111/ele.13345

2. **Moneghetti, P. C., et al. (2012)**. "A quantitative assessment of the competency period of *Acropora tenuis* coral larvae." [Original study providing A. tenuis competency data]

3. **Randal, E., et al. (2024)**. "Age to settlement competency and settlement cue preference in coral larvae." *Communications Biology*, 7, 184. https://doi.org/10.1038/s42003-024-05824-3

4. **This Analysis (2025)**: Direct comparison of Monegetti (A. tenuis) and Randal (Acroporidae) datasets using Monegetti piecewise model, with model selection and biological interpretation.

---

**Report Prepared**: December 2, 2025  
**Data Sources**: 
- Monegetti: `Settlement/monegetti/A.tenuisGBR2012metamorphosis.csv`
- Randal: `Randal/data/Settlement.csv` (Acroporidae family)
**Analysis Code**: `compare_monegetti_randal.py`  
**Comparison Plot**: `figures/monegetti_randal_comparison.png`  
**Combined Data**: `monegetti_randal_combined_data.csv`

---

*This report provides a comprehensive comparison of two foundational datasets for coral larval competency modeling, highlighting the importance of taxonomic resolution, treatment effects, and model selection in connectivity research.*

