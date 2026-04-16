#!/bin/bash
# Run full analysis pipeline for all families

echo "=========================================="
echo "Running Full Analysis Pipeline"
echo "=========================================="
echo ""

# Run the analysis
python monegetti_piecewise_model_refactored.py

echo ""
echo "=========================================="
echo "Analysis Complete!"
echo "=========================================="
echo ""
echo "Results saved to: results/"
echo "Figures saved to: figures_analysis/"
echo ""
echo "To compile LaTeX document:"
echo "  pdflatex analysis_documents.tex"
echo ""

