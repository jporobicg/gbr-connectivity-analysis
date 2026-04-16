#!/usr/bin/env python3
"""
Run Monegetti model for all families with progress tracking.
"""

import sys
import time
from monegetti_piecewise_model import main

if __name__ == "__main__":
    print("Starting Monegetti piecewise model analysis for ALL families...")
    print("This will analyze: Acroporidae, Poritidae, Merulinidae, Diploastreidae, Euphylliidae, Agariciidae, Lobophylliidae")
    print("Each family may take 2-5 minutes due to numerical integration.\n")
    
    start_time = time.time()
    
    # Run without test mode (all families)
    results, comparisons = main(test_mode=False)
    
    elapsed = time.time() - start_time
    print(f"\n{'='*70}")
    print(f"ANALYSIS COMPLETE!")
    print(f"Total time: {elapsed/60:.1f} minutes")
    print(f"{'='*70}")


