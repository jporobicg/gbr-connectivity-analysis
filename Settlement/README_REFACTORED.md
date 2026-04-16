# Refactored Analysis Pipeline - Usage Guide

## Overview

This is a complete refactoring of the larval competency modeling analysis into a clean, reproducible, publication-ready pipeline.

## New Structure

```
Settlement/
├── monegetti_piecewise_model_refactored.py  # Main orchestrator (entry point)
├── utility_tools.py                         # All helper functions
├── figures_analysis.py                      # All plotting code
├── generate_latex_tables.py                 # LaTeX table generation
├── analysis_documents.tex                   # LaTeX supplementary material
├── figures_analysis/                        # All publication figures
│   ├── theoretical_model_comparison.png
│   ├── family_*_fit.png (one per family)
│   └── all_families_comparison.png
├── results/                                 # Analysis results
│   ├── model_comparison_summary.csv
│   ├── *_timestamp.csv (one per family)
│   └── tables/
│       ├── model_comparison_table.tex
│       └── final_model_selection.tex
└── old_settlement_analysis/                 # Legacy code (archived)
```

## Quick Start

### Run Full Analysis

```bash
cd Settlement
python monegetti_piecewise_model_refactored.py
```

### Run in Test Mode (First Family Only)

```bash
python monegetti_piecewise_model_refactored.py --test
```

## File Responsibilities

### 1. `monegetti_piecewise_model_refactored.py` (Main Orchestrator)

**Purpose**: Main entry point, orchestrates full workflow

**Responsibilities**:
- Defines all model formulations (Logistic, Gompertz, Weibull, Monegetti)
- Controls model fitting workflow
- Orchestrates: data → fit → diagnostics → outputs
- Calls utility and plotting functions
- NO duplicated logic

**Key Functions**:
- `fit_all_models()`: Fit all candidate models
- `select_best_model()`: Select best model by AICc
- `analyze_family()`: Analyze single family
- `run_full_analysis()`: Complete pipeline

### 2. `utility_tools.py` (Helper Functions)

**Purpose**: All shared utility functions

**Responsibilities**:
- Data loading and preprocessing
- Likelihood calculations
- AICc computation
- Statistical diagnostics (overdispersion, residuals)
- File utilities
- Parameter constraints
- TC50 calculations
- Model functions (Logistic, Gompertz, Weibull, Monegetti)

**Key Functions**:
- `load_settlement_data()`: Load Randal data
- `load_monegetti_data()`: Load Monegetti data
- `fit_*_binomial()`: Fit each model type
- `calculate_aicc()`: Compute AICc
- `estimate_overdispersion()`: Calculate phi
- `calculate_residuals()`: Compute residuals

### 3. `figures_analysis.py` (All Plotting)

**Purpose**: Generate all publication-ready figures

**Responsibilities**:
- Theoretical model comparison figure
- Family-specific fits (data + model + diagnostics)
- All families comparison
- All figures saved to `figures_analysis/` folder
- Reproducible from code only

**Key Functions**:
- `create_theoretical_comparison_figure()`: Theoretical model shapes
- `create_family_figure()`: Single family fit with diagnostics
- `create_all_families_figure()`: All families comparison

### 4. `analysis_documents.tex` (Supplementary Material)

**Purpose**: LaTeX document for supplementary material

**Structure**:
1. Data comparison (why re-analysis needed)
2. Modeling decay (model descriptions, theoretical comparison)
3. Fitting the data (model comparison, diagnostics, final selection)
4. References to all figures and tables

## Outputs

### Results Directory (`results/`)

- `model_comparison_summary.csv`: Summary table with best model for each family
- `*_timestamp.csv`: Individual family results (one per family)
- `tables/`: LaTeX tables for document

### Figures Directory (`figures_analysis/`)

- `theoretical_model_comparison.png`: Theoretical shapes of all models
- `family_*_fit.png`: Individual family fits with diagnostics
- `all_families_comparison.png`: All families comparison

## Model Selection

The analysis fits all four models to each family:
1. **Logistic** (3 parameters)
2. **Gompertz** (3 parameters)
3. **Weibull** (4 parameters)
4. **Monegetti** (6 parameters)

The best model is selected based on **AICc** (Akaike Information Criterion corrected for small samples). Only statistically defensible (best) models are reported in the final document.

## Reproducibility

The entire analysis is reproducible from raw data:

1. Load data: `load_settlement_data(data_path)`
2. Fit models: `fit_all_models(df)`
3. Select best: `select_best_model(results)`
4. Generate figures: `figures_analysis.create_*_figure()`
5. Generate tables: `generate_latex_tables.generate_*_table()`

All random seeds are fixed (where applicable) to ensure reproducibility.

## Code Standards

- **PEP8 compliant**: All code follows PEP8 style guide
- **Type hints**: Functions include type hints where appropriate
- **Docstrings**: All functions have clear docstrings
- **Explicit imports**: No wildcard imports
- **Clear structure**: Logical organization and separation of concerns

## Legacy Code

All legacy scripts and outputs have been moved to `old_settlement_analysis/` for reference. The new structure uses only:
- `monegetti_piecewise_model_refactored.py`
- `utility_tools.py`
- `figures_analysis.py`
- `generate_latex_tables.py`

## Next Steps

1. Run full analysis: `python monegetti_piecewise_model_refactored.py`
2. Review results in `results/` directory
3. Review figures in `figures_analysis/` directory
4. Compile LaTeX document: `pdflatex analysis_documents.tex`
5. Verify all outputs match regenerated results

## Troubleshooting

### Import Errors
- Ensure all modules are in the same directory
- Check Python path includes Settlement directory

### Figure Generation Fails
- Check that `figures_analysis.py` is in the same directory
- Verify matplotlib and seaborn are installed

### Model Fitting Fails
- Check data requirements (minimum replicates, larvae)
- Review parameter bounds in `utility_tools.get_monegetti_bounds()`

## Contact

For questions about the refactored analysis, refer to:
- Code comments in each module
- This README
- `REFACTORING_PLAN.md` for detailed migration notes

