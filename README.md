# GBR connectivity analysis

<p align="center">
  <img src="assets/repo-icon.svg" alt="GBR connectivity analysis icon" width="96" height="96"/>
</p>

Code and workflows for analysing **Great Barrier Reef (GBR) larval connectivity**: comparing single connectivity estimates with bootstrap ensembles, summarising network structure, uncertainty, and ecological interpretation metrics.

This repository is intended for **reproducible analysis scripts** only. Large **datasets**, **raster outputs**, and **generated figures** are excluded via `.gitignore` and should be kept locally or distributed separately.

## Purpose

- **Connectivity comparison** (`connectivity_comparison/`): end-to-end pipeline (`run_analysis.py`) that reads single-matrix and bootstrap NetCDF products, computes matrix and node-level metrics, rank and profile stability, concentration and backbone diagnostics, coverage-style uncertainty, ordination, and writes tables plus optional plots and markdown reports.
- **Broader scripts** (`Codes/`, `Connectivity/`, `Settlement/`, `Spawning/`): notebooks and utilities for mapping, clustering reefs, plotting connectivity matrices, and related analyses used in the wider GBR modelling context.

## Repository structure

```text
.
├── README.md                 # This file
├── requirements.txt          # Core Python dependencies (notebooks & scripts)
├── assets/
│   └── repo-icon.svg         # Repository icon (network / reef motif)
├── Codes/                    # Standalone scripts, notebooks, plotting helpers
├── Connectivity/             # Connectivity-related utilities
├── connectivity_comparison/  # Main reproducible single vs bootstrap comparison workflow
│   ├── README.md             # Detailed workflow documentation
│   ├── config.py
│   ├── run_analysis.py       # CLI entry point
│   ├── requirements.txt      # Dependencies for the comparison pipeline
│   ├── src/                  # Metrics, I/O, plotting, reporting
│   └── notebooks/            # Optional exploration
├── Documents/                # Notes and presentation materials (figures ignored)
├── Settlement/               # Settlement-related analysis code
├── Spawning/                 # Spawning-related scripts
└── Tex_report/               # LaTeX catalog of analyses (figure copies ignored)
```

Ignored paths (see `.gitignore`) include: `datasets/`, `**/figures/`, `connectivity_comparison/outputs/`, binary geospatial and NetCDF data, and LaTeX build products under `Tex_report/`.

## Quick start (connectivity comparison)

```bash
cd connectivity_comparison
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
# Configure paths in config or pass CLI args, then:
python run_analysis.py --bootstrap-path /path/to/bootstrap.nc --single-path /path/to/single.nc --output-dir ./outputs
```

See `connectivity_comparison/README.md` for full options, outputs, and interpretation notes.

## Requirements

- Python 3.10+ recommended.
- Root `requirements.txt` lists common scientific stack packages; `connectivity_comparison/requirements.txt` targets the comparison workflow specifically.

## Citation

If you use this code in a publication, please cite the relevant GBR modelling or connectivity study and this repository as appropriate.

## License

Add a `LICENSE` file if you want to specify terms (e.g. MIT, CC-BY). Until then, all rights reserved unless you state otherwise.
