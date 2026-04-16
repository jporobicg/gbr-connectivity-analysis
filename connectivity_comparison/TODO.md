# TODO

## Data Gaps And Assumptions

- Confirm the ancillary distance fields in the connectivity NetCDF files. The current files behave as if `direction` is continuous distance in km and `distance` is a 0-35 distance-bin index, despite the variable attributes saying the opposite.
- Confirm that reef ordering in `datasets/reefs/Reefs2024.csv` matches the source/sink ordering in every family NetCDF file. The workflow currently assumes row order alignment because counts match at 3806 reefs.
- Replace the default region mapping if a management-specific reef-to-region table is available. The current regional exchange summaries default to `AREA_DESCR`, and community anchor profiles default to `AIMS_sector`.
- Add more bootstrap family NetCDF files. The workflow now compares all discovered single-family products, but bootstrap-aware family uncertainty remains limited by the files currently present in the repository.

## Method Refinements

- Evaluate whether a dedicated directed community method such as Infomap or Leiden-on-derived-flow-graphs is worth adding once a scalable dependency path is acceptable. The current implementation uses directed anchor-profile clustering for tractability.
- Revisit the stepping-stone metric if exact or sampled weighted betweenness becomes computationally practical for the full GBR network. The current bridge score is a brokerage-style approximation based on inter-anchor participation and throughflow.
- Consider caching time-mean family summaries so cross-family comparison does not have to recompute them from raw NetCDF files every run.
- Consider caching the selected-time ecological threshold robustness pass. It is intentionally limited to representative times, but it is still a second blockwise pass over the bootstrap data.

## Reporting And Interpretation

- Review whether the default long-distance classes (`0-100`, `100-300`, `>300 km`) are the right ecological bins for every family.
- If a paper uses only a subset of figures, keep the core set listed in the README and analysis catalog, then demote the rest to supporting or supplement material.
