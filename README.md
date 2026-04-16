# GBR connectivity analysis

**Repository:** [github.com/jporobicg/gbr-connectivity-analysis](https://github.com/jporobicg/gbr-connectivity-analysis)

<p align="center">
  <img src="assets/repo-icon.svg" alt="GBR connectivity analysis icon" width="96" height="96"/>
</p>

Analysis code for **Great Barrier Reef (GBR) larval connectivity**: comparing single connectivity estimates with bootstrap ensembles, mapping and plotting helpers, settlement and spawning workflows, and reproducible reports.

Large **datasets**, **generated figures**, and **tabular exports** stay local (see `.gitignore`).

---

## 1. Connectivity

**Path:** [`Connectivity/`](Connectivity/)

**Purpose:** Main hub for connectivity science code. It merges the former `connectivity_comparison/` pipeline with the older `Connectivity/` plotting helpers.

| What | Role |
|------|------|
| `run_analysis.py`, `config.py`, `src/` | Reproducible **single vs bootstrap** NetCDF comparison: matrix metrics, ranks, profiles, concentration, stable edges, coverage, ordination, regional and family extensions, markdown reports and plots. |
| `requirements.txt` | Python dependencies for that pipeline. |
| `notebooks/` | Optional exploration (e.g. `optional_exploration.ipynb`). |
| `make_plots.py`, `plot_images_functions.py`, `*.ipynb` | Legacy or ad hoc **plotting and parcel** workflows. |
| `outputs/` | **Generated** tables, figures, intermediates (not tracked in git). |

**Quick start:** see [`Connectivity/README.md`](Connectivity/README.md).

```bash
cd Connectivity
pip install -r requirements.txt
python run_analysis.py --help
```

---

## 2. Codes

**Path:** [`Codes/`](Codes/)

**Purpose:** Standalone notebooks and scripts used across GBR modelling (clustering reefs, connectivity-by-area summaries, kernel helpers, shared plotting utilities).

---

## 3. Settlement

**Path:** [`Settlement/`](Settlement/)

**Purpose:** Larval settlement and competency modelling, comparisons with literature parametrisations, R reproduction material (`Randal_github/`), and analysis notes (Markdown / Org).

---

## 4. Spawning

**Path:** [`Spawning/`](Spawning/)

**Purpose:** Spawning timing and related **R** analyses for coral reproductive schedules.

---

## 5. Documents

**Path:** [`Documents/`](Documents/)

**Purpose:** Manuscript sources, prompts, and internal notes (presentation figures are not committed).

---

## 6. Tex_report

**Path:** [`Tex_report/`](Tex_report/)

**Purpose:** LaTeX version of the analysis catalog (`analysis_catalog.tex`). Local `figures/` copies and PDF build artifacts are gitignored.

---

## 7. Data layout (local only)

**Path:** [`datasets/`](datasets/) (not in git)

Place NetCDF connectivity matrices, reef shapefiles / CSV metadata, and kernels here as expected by `Connectivity/config.py` defaults or your own paths.

---

## Root `requirements.txt`

Lightweight stack for notebooks and scripts used outside the pipeline. For the **comparison workflow**, prefer:

`Connectivity/requirements.txt`

---

## License

This project is released under the [MIT License](LICENSE).

---

## Citation

If you use this code, cite the relevant GBR or connectivity publication and link to this repository.

---

## GitHub (maintainers)

1. **Repository description:** e.g. *“GBR larval connectivity analysis: single vs bootstrap ensembles, settlement, spawning.”*  
   Settings → General → **Repository name** (already set) and **Description**.

2. **Social preview:** Settings → General → **Social preview** → upload  
   [`assets/social-preview.png`](assets/social-preview.png)  
   (1200×630 recommended for Open Graph / link cards.)

3. **Topics:** e.g. `great-barrier-reef`, `connectivity`, `larval-dispersal`, `python`, `netcdf`, `bootstrap`.

4. **Default branch:** `main` (already used for pushes).
