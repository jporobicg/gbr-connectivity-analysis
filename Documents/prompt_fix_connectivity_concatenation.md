# Prompt: Fix connectivity NetCDF concatenation bug

Use this prompt in another AI session (with access to the codebase that **produces** the `*_single.nc` files from the per-run NetCDFs) to locate and fix the concatenation bug.

---

## Task

Find the code that **concatenates** per-run connectivity NetCDF files (e.g. `connectivity_2015-10-29_merulinidae_run0.nc`, `connectivity_2015-10-29_acroporidae_run0.nc`, possibly multiple runs and dates) into **single** combined files (e.g. `connectivity_merulinidae_single.nc`, `connectivity_acroporidae_single.nc`). There is a bug in that concatenation: the output files have far fewer sources with valid connectivity than the input files.

Fix the bug so that the combined files preserve the same number of sources with valid connections as the originals.

---

## Observed behaviour (evidence of the bug)

**Definition of “valid connection” (for counting):**  
A matrix cell is valid if: value is finite, not the NetCDF fill value (e.g. `_FillValue` ~ 9.97e36), and **0 < value < 1e16**.

**“Sources with valid connections”** = number of rows (sources) that have at least one valid cell in that row.

### Original (per-run) files — correct

- **connectivity_2015-10-29_merulinidae_run0.nc**  
  - Dimensions: `(source=3806, sink=3806, treatment=1, sample=10)` (no time dim).  
  - **Sources with valid connections:** 3,806 (all sources have at least one valid value in the first slice).  
  - **Sources with at least one value > 0:** 3,618.

- **connectivity_2015-10-29_acroporidae_run0.nc**  
  - Similar layout.  
  - **Sources with valid connections:** 3,806 (or similar, e.g. 3,605 with value > 0).

So the **originals have on the order of 3,500+ sources with valid connections**.

### Concatenated (“single”) files — wrong

- **connectivity_merulinidae_single.nc**  
  - Dimensions: e.g. `(source=3806, sink=3806, time=27, treatment=3, sample=100)`.  
  - **Only time=0** has valid data.  
  - **Sources with valid connections at time=0:** only **27** (or 26 in one slice).  
  - For **time=1, 2, …, 26**: **0** sources with valid connections.

- **connectivity_acroporidae_single.nc**  
  - Similar structure.  
  - **Sources with valid connections (e.g. at time=0):** only **27**.

So the **concatenated results have only ~27 sources with valid connections** instead of 3,500+.

---

## What to look for in the concatenation code

1. **Dimension / axis used for concatenation**  
   - Are runs (or times) concatenated along the correct dimension (e.g. `time`)?  
   - Is data from runs/times being written into the wrong dimension or the wrong slice (e.g. only the first time step ever filled)?

2. **Source/sink indexing**  
   - When merging multiple runs or files, is the **source** (and sink) index preserved correctly?  
   - Could the code be writing each run’s full matrix into a **subset** of sources (e.g. only the first 27) instead of all 3,806?

3. **Fill values and missing data**  
   - Is the concatenation initialising a large array with fill/missing and then only filling a small part (e.g. one time or one block of sources)?  
   - Are slices that should contain data from other runs/times left as fill?

4. **Time vs run vs sample**  
   - How are “run”, “time”, “treatment”, and “sample” defined in the original files vs the single file?  
   - Is “time” in the single file intended to correspond to different release dates or runs? If so, only time=0 being non-fill suggests only the first run/date is being copied and the rest are never written.

5. **xarray concat/merge**  
   - If using `xr.concat()` or `xr.merge()`, check the `dim=` argument and that all input datasets are actually included and aligned on `source` and `sink`.

---

## Files and paths (for reference)

- **Original (per-run) files:**  
  `datasets/connectivity_matrices/connectivity_2015-10-29_merulinidae_run0.nc`,  
  `datasets/connectivity_matrices/connectivity_2015-10-29_acroporidae_run0.nc`,  
  (and possibly other runs/dates in the same or sibling directory.)

- **Incorrect output (concatenated):**  
  `datasets/connectivity_matrices/connectivity_merulinidae_single.nc`,  
  `datasets/connectivity_matrices/connectivity_acroporidae_single.nc`.

The concatenation script may live in a **different repository** (e.g. a preprocessing or GBR modelling pipeline), not in the Connectivity_analysis repo that only **reads** these NetCDFs for plotting and analysis.

