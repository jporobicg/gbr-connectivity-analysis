# Full Analysis Results Summary

## Analysis Completed Successfully ✅

**Date**: December 17, 2025  
**Total Families Analyzed**: 7  
**Total Observations**: 11,759  
**Status**: Complete

## Model Selection Results

### Summary by Family

| Family | Best Model | N Replicates | N Larvae | AICc | TC50 (days) |
|--------|------------|--------------|----------|------|-------------|
| Acroporidae | **Weibull** | 4,730 | 45,123 | 56,693.4 | N/A |
| Agariciidae | **Weibull** | 167 | 1,542 | 599.9 | N/A |
| Diploastreidae | **Weibull** | 539 | 5,329 | 5,188.8 | N/A |
| Euphylliidae | **Weibull** | 496 | 4,631 | 5,815.6 | N/A |
| Lobophylliidae | **Monegetti** | 760 | 6,560 | 2,984.0 | 8.79 |
| Merulinidae | **Monegetti** | 3,880 | 38,495 | 38,033.0 | 34.46 |
| Poritidae | **Monegetti** | 1,187 | 11,612 | 10,726.7 | 4.77 |

### Key Findings

1. **Weibull Model Selected** (4 families):
   - Acroporidae (largest dataset: 45K larvae)
   - Agariciidae (smallest dataset: 1.5K larvae)
   - Diploastreidae
   - Euphylliidae
   
   These families show simpler competency dynamics that are well-captured by the 4-parameter Weibull model.

2. **Monegetti Model Selected** (3 families):
   - Lobophylliidae (TC50 = 8.79 days)
   - Merulinidae (TC50 = 34.46 days) - Strongest evidence (Δ AICc = 648.7)
   - Poritidae (TC50 = 4.77 days)
   
   These families show complex two-phase competency dynamics requiring the 6-parameter Monegetti model.

3. **Model Selection Evidence**:
   - **Merulinidae**: Strongest preference for Monegetti (Δ AICc = 648.7 vs Weibull)
   - **Lobophylliidae**: Moderate preference (Δ AICc = 24.5 vs Weibull)
   - **Poritidae**: Weak preference (Δ AICc = 1.9 vs Logistic)
   - Other families: Clear preference for Weibull over Monegetti

## Generated Outputs

### Results Files
- ✅ `results/model_comparison_summary.csv` - Summary table
- ✅ `results/*_timestamp.csv` - Individual family results (7 files)
- ✅ `results/tables/model_comparison_table.tex` - Full model comparison LaTeX table
- ✅ `results/tables/final_model_selection.tex` - Final selection LaTeX table

### Figures
- ✅ ✅ `figures_analysis/theoretical_model_comparison.png` - Theoretical model shapes
- ✅ `figures_analysis/family_Acroporidae_fit.png` - Acroporidae fit with diagnostics
- ✅ `figures_analysis/family_Agariciidae_fit.png` - Agariciidae fit with diagnostics
- ✅ `figures_analysis/family_Diploastreidae_fit.png` - Diploastreidae fit with diagnostics
- ✅ `figures_analysis/family_Euphylliidae_fit.png` - Euphylliidae fit with diagnostics
- ✅ `figures_analysis/family_Lobophylliidae_fit.png` - Lobophylliidae fit with diagnostics
- ✅ `figures_analysis/family_Merulinidae_fit.png` - Merulinidae fit with diagnostics
- ✅ `figures_analysis/family_Poritidae_fit.png` - Poritidae fit with diagnostics
- ✅ `figures_analysis/all_families_comparison.png` - All families comparison

**Total**: 9 publication-ready figures

## Statistical Summary

### Overdispersion (φ)
- All families show overdispersion (φ > 1)
- Range: 3.95 (Lobophylliidae) to 6.97 (Diploastreidae)
- Indicates extra-binomial variation, which is expected for biological data

### Sample Sizes
- **Largest**: Acroporidae (45,123 larvae, 4,730 replicates)
- **Smallest**: Agariciidae (1,542 larvae, 167 replicates)
- All families have sufficient data for robust model fitting

## Biological Interpretation

1. **Simple Competency Dynamics** (Weibull families):
   - Four families show straightforward competency development
   - Single-phase model sufficient
   - No evidence for two-phase decay structure

2. **Complex Competency Dynamics** (Monegetti families):
   - Three families show two-phase competency loss
   - Early phase: Weibull decay
   - Late phase: Exponential decay
   - Change point (Tcp) marks transition between phases

3. **TC50 Values** (Monegetti families only):
   - **Poritidae**: 4.77 days (fastest)
   - **Lobophylliidae**: 8.79 days
   - **Merulinidae**: 34.46 days (slowest, longest competency window)

## Next Steps

1. ✅ **Analysis Complete** - All families analyzed
2. ⏳ **Review Figures** - Check all family fits for quality
3. ⏳ **Compile LaTeX** - Generate PDF supplementary material
4. ⏳ **Final Review** - Verify all results are publication-ready

## Files Generated

All outputs are in:
- **Results**: `results/` directory
- **Figures**: `figures_analysis/` directory
- **Log**: `full_analysis.log` (if saved)

---

**Analysis Status**: ✅ **COMPLETE**  
**Ready for**: Publication review and LaTeX compilation

