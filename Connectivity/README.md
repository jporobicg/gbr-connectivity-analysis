# Connectivity (workflow and tools)

This folder holds the **single vs bootstrap connectivity comparison** pipeline (`run_analysis.py`, `config.py`, `src/`) together with **older plotting helpers** (`make_plots.py`, `plot_images_functions.py`, notebooks).

It compares a single connectivity estimate against a bootstrap connectivity ensemble stored in NetCDF files.

The central question is:

> Does the single matrix lead to the same ecological interpretation as the bootstrap ensemble, and how much extra uncertainty information does the bootstrap approach provide?

## Layout of this folder

```text
Connectivity/
├── README.md                 # This file
├── requirements.txt          # Pipeline dependencies
├── config.py                 # Defaults and CLI merge
├── run_analysis.py           # Main CLI entry point
├── TODO.md
├── codex_connectivity_extended_prompt.md
├── make_plots.py             # Helper plotting scripts
├── plot_images_functions.py
├── Plot_parcels.ipynb
├── Untitled-1.ipynb
├── src/
│   ├── io_utils.py
│   ├── preprocessing.py
│   ├── metrics.py
│   ├── network_metrics.py
│   ├── extended_metrics.py
│   ├── metadata.py
│   ├── high_level_metrics.py
│   ├── family_comparison.py
│   ├── plotting.py
│   └── reporting.py
├── notebooks/
│   └── optional_exploration.ipynb
└── outputs/                  # Created when you run the pipeline (gitignored)
    ├── tables/
    ├── figures/
    ├── intermediate/
    ├── analysis_catalog.md
    └── report_summary.md
```

## What The Workflow Covers

The workflow now combines the original matrix/network comparison with extra analyses focused on:

- system structure
- strength and concentration of connectivity
- overrepresentation and stable backbones
- uncertainty and robustness
- community-structure uncertainty
- stepping-stone / bridge reefs
- temporal vs bootstrap variance decomposition
- family-level comparison
- distance-dependent uncertainty
- regional exchange
- mean-matrix bias diagnostics
- spatial hotspots of disagreement
- ecologically interpretable threshold robustness

For each time step it computes:

- the single matrix
- the bootstrap mean matrix
- bootstrap SD, CV, and optional quantiles
- matrix similarity metrics
- node strength and rank stability
- destination/origin profile divergence
- simple higher-level grouping stability
- concentration and upper-tail metrics
- stable-edge frequency and top-k churn
- local entropy, dominance, and effective partner counts
- coverage-style uncertainty summaries
- discrepancy vs within-bootstrap variability
- bootstrap ordination with the single and bootstrap mean overlaid
- consensus connectivity regions and co-assignment summaries
- bridge-reef scores and removal-impact summaries
- distance-decay envelopes and long-distance fractions
- region-to-region exchange and regional source/sink roles
- interpretation-focused mean-bias diagnostics
- spatial clustering summaries for disagreement hotspots
- family summaries across all detected single-matrix products
- variance-component tables for interpretation metrics

## Default Input Files

By default the project uses:

- `datasets/connectivity_matrices/connectivity_acroporidae.nc`
- `datasets/connectivity_matrices/connectivity_acroporidae_single.nc`
- `datasets/reefs/Reefs2024.csv`

You can change these in [config.py](config.py) or pass them on the command line.

## Suggested order to run

Run these from the `Connectivity/` directory (so `config.py` and `src/` resolve correctly).

1. Install dependencies:

```bash
cd Connectivity
pip install -r requirements.txt
```

2. Run a quick pilot on one or two time steps:

```bash
python run_analysis.py --time-indices 1,2 --skip-plots
```

3. Run the full analysis:

```bash
python run_analysis.py
```

### SLURM (one job per family)

From the **repository root**, set the directory that contains `connectivity_*.nc` and submit the array job (see `slurm/run_family_array.slurm` for memory/time defaults you can tune):

```bash
bash Connectivity/slurm/submit_family_array.sh /path/to/your/netcdf_directory
```

Or manually:

```bash
export CONNECTIVITY_MATRICES_DIR=/path/to/your/netcdf_directory
mkdir -p Connectivity/slurm/logs
sbatch Connectivity/slurm/run_family_array.slurm
```

Outputs go to `Connectivity/outputs_<family>/` so families do not overwrite each other.

4. Open the main outputs:

- [outputs/analysis_catalog.md](outputs/analysis_catalog.md)
- [outputs/report_summary.md](outputs/report_summary.md)
- [TODO.md](TODO.md)
- [outputs/tables](outputs/tables)
- [outputs/figures](outputs/figures)

## Main Output Tables

Original core outputs:

- `per_time_matrix_comparison.csv`
- `per_time_network_metrics.csv`
- `top_link_overlap.csv`
- `rank_stability_summary.csv`
- `top_changing_sources.csv`
- `top_changing_sinks.csv`
- `top_changing_links.csv`
- `edge_ci_summary.csv`
- `node_ci_summary.csv`
- `sample_similarity_summary.csv`
- `ordination_coordinates.csv`

New extended outputs:

- `source_rank_comparison.csv`
- `sink_rank_comparison.csv`
- `source_rank_uncertainty.csv`
- `sink_rank_uncertainty.csv`
- `profile_divergence_summary.csv`
- `source_profile_changes.csv`
- `sink_profile_changes.csv`
- `grouping_similarity_summary.csv`
- `grouping_size_summary.csv`
- `concentration_metrics_summary.csv`
- `lorenz_curve_points.csv`
- `tail_metrics_summary.csv`
- `bootstrap_sample_tail_summary.csv`
- `node_strength_uncertainty.csv`
- `top_uncertain_nodes.csv`
- `stable_edge_summary.csv`
- `stable_edge_detail.csv`
- `local_node_metrics.csv`
- `backbone_summary.csv`
- `coverage_summary.csv`
- `discrepancy_vs_variability.csv`

Higher-level ecological outputs:

- `community_structure_summary.csv`
- `community_assignments.csv`
- `community_coassignment_summary.csv`
- `community_bootstrap_similarity.csv`
- `bridge_reef_summary.csv`
- `bridge_reef_time_summary.csv`
- `distance_uncertainty_curve.csv`
- `distance_uncertainty_classes.csv`
- `regional_exchange_summary.csv`
- `regional_role_summary.csv`
- `regional_exchange_time_summary.csv`
- `mean_bias_summary.csv`
- `spatial_clustering_summary.csv`
- `spatial_hotspots.csv`
- `variance_decomposition_summary.csv`
- `variance_decomposition_detail.csv`
- `family_summary.csv`
- `family_similarity_summary.csv`
- `family_topk_overlap.csv`
- `family_ordination.csv`
- `family_uncertainty_summary.csv`
- `ecological_threshold_summary.csv`

## Main Figures

Core figures:

- `time_series_matrix_metrics.png`
- `threshold_sensitivity.png`
- `uncertainty_time_series.png`
- `heatmap_panel_time_XX.png`
- `scatter_time_XX.png`
- `distributions_time_XX.png`
- `rank_panel_time_XX.png`
- `ordination_time_XX.png`

New figures:

- `rank_correlation_time_series.png`
- `rank_change_bars_time_XX.png`
- `profile_divergence_hist_time_XX.png`
- `profile_top_changes_time_XX.png`
- `grouping_similarity_time_series.png`
- `group_size_summary_time_XX.png`
- `lorenz_curves_time_XX.png`
- `concentration_time_series.png`
- `tail_metric_time_series.png`
- `top_edge_strength_distribution_time_XX.png`
- `node_strength_uncertainty_time_XX.png`
- `node_strength_cv_hist_time_XX.png`
- `topk_overlap_time_series.png`
- `edge_churn_time_XX.png`
- `stable_edge_counts.png`
- `stable_edge_overlap.png`
- `entropy_strength_time_XX.png`
- `dominance_ratio_time_XX.png`
- `effective_links_distribution_time_XX.png`
- `backbone_overlap_time_series.png`
- `coverage_rates_time_series.png`
- `coverage_examples_time_XX.png`
- `single_vs_bootstrap_variability.png`
- `topk_probability_time_XX.png`
- `rank_uncertainty_intervals_time_XX.png`
- `community_similarity_time_series.png`
- `community_consensus_time_XX.png`
- `community_stability_time_XX.png`
- `community_coassignment_time_XX.png`
- `bridge_map_time_XX.png`
- `bridge_uncertainty_time_XX.png`
- `bridge_strength_scatter_time_XX.png`
- `bridge_removal_impact_time_XX.png`
- `distance_decay_time_XX.png`
- `distance_uncertainty_time_XX.png`
- `long_distance_fraction_time_series.png`
- `regional_exchange_time_XX.png`
- `regional_roles_time_XX.png`
- `regional_pair_intervals_time_XX.png`
- `mean_bias_summary.png`
- `spatial_clustering_summary.png`
- `spatial_hotspots_time_XX.png`
- `variance_decomposition.png`
- `family_similarity_heatmap.png`
- `family_overlap_summary.png`
- `family_ordination.png`
- `family_uncertainty_panel.png`
- `ecological_threshold_robustness.png`

## What The New Analyses Mean

`Source and sink rank uncertainty`
: Tests whether the same reefs stay important once the full bootstrap distribution is considered.

`Profile divergence`
: Tests whether reefs connect to the same places even if total connectivity stays similar.

`Grouping similarity`
: Tests whether the system shows similar large-scale organization under both approaches.

`Concentration and tail metrics`
: Tests whether the single estimate is sharper or more extreme than the bootstrap mean.

`Stable-edge summaries`
: Tests which pathways stay important across many bootstrap samples.

`Local entropy and dominance`
: Tests whether nodes have diffuse or highly focused connectivity patterns.

`Coverage and discrepancy analyses`
: Tests whether the single estimate behaves like a plausible bootstrap realization or falls outside the ensemble.

## Practical Notes

- The bootstrap dataset is too large to hold as one full `(time, source, sink, sample)` array.
- The script processes one time step at a time and one source/sink block at a time.
- The default `block_size` is `512`.
- Quantiles are useful but slower than means and SDs. Use `--skip-quantiles` for a faster pilot run.
- Some analyses use sampled edges for tractability. Those are intended for uncertainty scaling and plotting, not to replace the exact node-level summaries.
- The example Acroporidae files contain a metadata mismatch at time index `0`: the bootstrap file has `NaT`, while the single file has `2015-10-29`. The workflow keeps going and uses the valid time label.
- The current NetCDF ancillaries appear to store the continuous reef-to-reef distance matrix in `direction`, while `distance` behaves like a distance-bin index. The loader detects this and records the mismatch in `TODO.md`.
- Community structure is implemented as directed anchor-profile clustering using reef flow profiles to coarse reef sectors. This is a practical large-network approximation, not an exact Infomap or Leiden partition.
- Bridge-reef importance is implemented with inter-anchor participation and brokerage-style throughflow rather than exact all-node weighted betweenness. This keeps the workflow tractable at GBR scale.
- Family comparison uses all detected single-matrix family files and adds bootstrap-aware uncertainty summaries when bootstrap files are available. At the moment the repository only contains one full bootstrap family product.

## Core, Supporting, Optional

Core paper analyses:

- matrix comparison and uncertainty time series
- source/sink rank stability
- community similarity and consensus maps
- stepping-stone maps and uncertainty
- distance-decay envelopes and long-distance fractions
- regional exchange summaries
- mean-bias diagnostics
- variance decomposition

Supporting analyses:

- stable edges and backbone overlap
- profile divergence
- concentration and tail metrics
- community co-assignment heatmaps
- regional pair intervals
- family similarity heatmaps
- ecological threshold robustness

Optional analyses:

- block heatmaps
- edge scatter and raw distributions
- entropy, dominance, and effective-link views
- spatial hotspot maps
- family ordination

## Checklist Before Using With Real Data

- Confirm the connectivity variable is still `connectivity`.
- Confirm `treatment` is length 1 or change the treatment selection.
- Decide whether the diagonal should be interpreted as self-recruitment.
- Check whether the default thresholds are biologically meaningful for your species and units.
- Check whether tiny non-zero values should count as active links.
- Confirm the source and sink ordering matches your reef metadata.
- Inspect one pilot run before running all time steps.

## Adapting To Other Species

- Change the file paths.
- Adjust thresholds if the connectivity scale differs.
- Adjust `top_k_edges` and `rank_top_k_values` if the number of biologically important links or reefs differs.
- Revisit the ecological meaning of self-recruitment and local dominance if life histories differ.

## Extending To More Than Two Products

- Generalize the script to loop over named products instead of assuming only `single` and `bootstrap`.
- Reuse the same matrix, rank, concentration, and uncertainty summaries for each pairwise comparison.
- Add a product label column and write a comparison matrix across all products.

## Short Interpretation Guide

Start with the matrix comparison, rank correlation, top-link overlap, and discrepancy-vs-variability figures. If those stay stable, the ecological story is probably robust. Then use the uncertainty, stable-edge, concentration, and profile figures to see what extra information the bootstrap ensemble adds. The main goal is not perfect edge-by-edge equality. It is whether the same ecological conclusions hold once uncertainty is represented explicitly.
