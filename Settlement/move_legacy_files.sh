#!/bin/bash
# Move legacy files to old_settlement_analysis/

# Legacy Python scripts (except the new refactored ones)
mv improved_competency_model.py old_settlement_analysis/ 2>/dev/null
mv compare_monegetti_randal.py old_settlement_analysis/ 2>/dev/null
mv plot_monegetti_results.py old_settlement_analysis/ 2>/dev/null
mv plot_all_families_with_uncertainty.py old_settlement_analysis/ 2>/dev/null
mv analyze_settlement_data.py old_settlement_analysis/ 2>/dev/null
mv run_monegetti_all_families.py old_settlement_analysis/ 2>/dev/null
mv parameter_uncertainty_analysis.py old_settlement_analysis/ 2>/dev/null
mv uncertainty_analysis_by_family.py old_settlement_analysis/ 2>/dev/null
mv hybrid_competency_model.py old_settlement_analysis/ 2>/dev/null
mv plot_improved_vs_original.py old_settlement_analysis/ 2>/dev/null
mv plot_settlement_by_family.py old_settlement_analysis/ 2>/dev/null
mv test_family_competence.py old_settlement_analysis/ 2>/dev/null
mv connectivity.py old_settlement_analysis/ 2>/dev/null

# Legacy output files
mv *.png old_settlement_analysis/ 2>/dev/null
mv *.csv old_settlement_analysis/ 2>/dev/null
mv *.pkl old_settlement_analysis/ 2>/dev/null
mv *.log old_settlement_analysis/ 2>/dev/null

# Keep only new structure files
echo "Legacy files moved to old_settlement_analysis/"
