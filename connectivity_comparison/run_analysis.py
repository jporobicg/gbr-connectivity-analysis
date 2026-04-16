#!/usr/bin/env python3
"""Run a reproducible comparison between single and bootstrap connectivity datasets."""

from __future__ import annotations

import argparse
import gc
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from config import build_runtime_config
from src.io_utils import (
    dataset_summary,
    ensure_output_directories,
    open_connectivity_dataset,
    prepare_connectivity_array,
    resolve_time_values,
    validate_datasets,
    write_dataframe,
    write_markdown,
)
from src.family_comparison import build_family_comparison, discover_family_datasets
from src.extended_metrics import (
    compute_backbone_summary,
    compute_concentration_and_tail,
    compute_discrepancy_vs_variability,
    compute_edge_churn,
    compute_grouping_summary,
    compute_local_node_metrics,
    compute_node_strength_uncertainty,
    compute_profile_divergence,
    compute_rank_uncertainty,
    compute_stable_edge_summary,
    compute_system_coverage_summary,
    sample_topk_tracker_update,
    summarize_local_node_metrics,
)
from src.metrics import (
    bootstrap_interval_summary,
    compute_matrix_comparison_metrics,
    compute_rank_stability,
    pca_ordination,
    summarize_sample_similarity,
    summarize_single_vs_bootstrap_ci,
    top_changing_links,
)
from src.high_level_metrics import (
    compute_bridge_reefs,
    compute_community_structure,
    compute_distance_uncertainty,
    compute_ecological_thresholds,
    compute_mean_bias,
    compute_regional_exchange,
    compute_spatial_hotspots,
    compute_variance_decomposition,
)
from src.metadata import infer_distance_matrix, load_reef_metadata
from src.network_metrics import gini, network_metric_rows, top_nodes_table
from src.preprocessing import (
    assign_sampled_edges_from_block,
    block_average,
    coefficient_of_variation,
    format_time_label,
    iter_matrix_blocks,
    npz_path,
    sampled_edge_indices,
    clean_array,
)
from src.reporting import build_analysis_catalog, build_report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bootstrap-path", help="Path to the bootstrap NetCDF file.")
    parser.add_argument("--single-path", help="Path to the single-matrix NetCDF file.")
    parser.add_argument("--output-dir", help="Directory for all tables, figures, and intermediate outputs.")
    parser.add_argument("--time-indices", help="Comma-separated list of time indices to analyse. Default: all.")
    parser.add_argument("--block-size", type=int, help="Source/sink block size used while reading the bootstrap file.")
    parser.add_argument("--skip-quantiles", action="store_true", help="Skip bootstrap quantiles and CI-based analyses.")
    parser.add_argument("--skip-plots", action="store_true", help="Skip plotting and only write tables/report.")
    return parser.parse_args()


def print_summary(summary: dict[str, object]) -> None:
    print(f"[summary] {summary['label']}")
    print(f"  path: {summary['path']}")
    print(f"  dims: {summary['dims']}")
    print(f"  preview time: {summary['time_preview_label']}")
    print(f"  finite fraction (preview): {summary['finite_fraction_preview']:.4f}")
    print(f"  zero fraction (preview): {summary['zero_fraction_preview']:.4f}")
    print(f"  preview min/max: {summary['min_preview']:.4g} / {summary['max_preview']:.4g}")
    print(f"  preview median non-zero: {summary['median_nonzero_preview']:.4g}")


def _diagonal_overlap(source_slice: slice, sink_slice: slice) -> tuple[np.ndarray, np.ndarray] | None:
    start = max(int(source_slice.start), int(sink_slice.start))
    stop = min(int(source_slice.stop), int(sink_slice.stop))
    if start >= stop:
        return None
    source_idx = np.arange(start, stop) - int(source_slice.start)
    sink_idx = np.arange(start, stop) - int(sink_slice.start)
    return source_idx.astype(np.int32), sink_idx.astype(np.int32)


def _family_name_from_single_path(path: Path) -> str:
    stem = path.stem
    if stem.startswith("connectivity_"):
        stem = stem[len("connectivity_") :]
    if stem.endswith("_single"):
        stem = stem[: -len("_single")]
    return stem


def _time_components(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty or "time_label" not in frame:
        return frame
    result = frame.copy()
    timestamps = pd.to_datetime(result["time_label"], errors="coerce")
    result["year"] = timestamps.dt.year
    result["spawning_period"] = timestamps.dt.strftime("%m")
    result.loc[timestamps.isna(), "year"] = np.nan
    result.loc[timestamps.isna(), "spawning_period"] = "unknown"
    return result


def analyse_time_step(
    bootstrap_da,
    single_da,
    time_index: int,
    time_label: str,
    config: dict,
    metadata_context: dict,
    distance_matrix: np.ndarray | None,
) -> dict[str, pd.DataFrame | Path | str]:
    """Analyse one time slice and write compact intermediate arrays."""
    n_source = int(single_da.sizes["source"])
    n_sink = int(single_da.sizes["sink"])
    n_samples = int(bootstrap_da.sizes["sample"])
    thresholds = [float(value) for value in config["thresholds"]]
    stable_thresholds = [float(value) for value in config["stable_edge_thresholds"]]
    edge_top_k_values = sorted(set(int(value) for value in config["top_k_edges"]))
    rank_top_k_values = sorted(set(int(value) for value in config["rank_top_k_values"]))
    max_edge_top_k = max(edge_top_k_values)
    community_anchor_onehot = np.asarray(metadata_context["community_anchor_onehot"], dtype=np.float32)
    region_onehot = np.asarray(metadata_context["region_onehot"], dtype=np.float32)
    n_anchor = int(community_anchor_onehot.shape[1])
    n_region = int(region_onehot.shape[1])

    distance_curve_edges = np.asarray(config["distance_curve_edges_km"], dtype=np.float32)
    distance_class_edges = np.asarray(config["distance_class_edges_km"], dtype=np.float32)
    distance_class_labels = list(config["distance_class_labels"])
    n_distance_curve = max(distance_curve_edges.size - 1, 0)
    n_distance_class = len(distance_class_labels)

    single_matrix = np.full((n_source, n_sink), np.nan, dtype=np.float32)
    bootstrap_mean = np.full_like(single_matrix, np.nan)
    bootstrap_std = np.full_like(single_matrix, np.nan)
    q_low = q_med = q_high = None
    if config["compute_quantiles"]:
        q_low = np.full_like(single_matrix, np.nan)
        q_med = np.full_like(single_matrix, np.nan)
        q_high = np.full_like(single_matrix, np.nan)

    row_sums_boot = np.zeros((n_samples, n_source), dtype=np.float64)
    col_sums_boot = np.zeros((n_samples, n_sink), dtype=np.float64)
    diagonal_boot = np.zeros(n_samples, dtype=np.float64)
    connectance_counts = np.zeros((n_samples, len(thresholds)), dtype=np.int64)
    threshold_presence_counts = {
        threshold: np.zeros((n_source, n_sink), dtype=np.uint8) for threshold in stable_thresholds
    }
    row_wlogw_boot = np.zeros((n_samples, n_source), dtype=np.float64)
    col_wlogw_boot = np.zeros((n_samples, n_sink), dtype=np.float64)
    row_max_boot = np.zeros((n_samples, n_source), dtype=np.float32)
    col_max_boot = np.zeros((n_samples, n_sink), dtype=np.float32)
    sample_top_values = np.full((n_samples, max_edge_top_k), -np.inf, dtype=np.float32)
    sample_top_indices = np.full((n_samples, max_edge_top_k), -1, dtype=np.int64)
    single_out_anchor = np.zeros((n_source, n_anchor), dtype=np.float64)
    single_in_anchor = np.zeros((n_source, n_anchor), dtype=np.float64)
    mean_out_anchor = np.zeros((n_source, n_anchor), dtype=np.float64)
    mean_in_anchor = np.zeros((n_source, n_anchor), dtype=np.float64)
    sample_out_anchor = np.zeros((n_samples, n_source, n_anchor), dtype=np.float32)
    sample_in_anchor = np.zeros((n_samples, n_source, n_anchor), dtype=np.float32)
    single_region_matrix = np.zeros((n_region, n_region), dtype=np.float64)
    mean_region_matrix = np.zeros((n_region, n_region), dtype=np.float64)
    sample_region_matrices = np.zeros((n_samples, n_region, n_region), dtype=np.float32)
    single_curve_sums = np.zeros(n_distance_curve, dtype=np.float64)
    mean_curve_sums = np.zeros(n_distance_curve, dtype=np.float64)
    sample_curve_sums = np.zeros((n_samples, n_distance_curve), dtype=np.float64)
    curve_counts = np.zeros(n_distance_curve, dtype=np.int64)
    single_curve_positive = np.zeros(n_distance_curve, dtype=np.int64)
    mean_curve_positive = np.zeros(n_distance_curve, dtype=np.int64)
    sample_curve_positive = np.zeros((n_samples, n_distance_curve), dtype=np.int64)
    single_class_sums = np.zeros(n_distance_class, dtype=np.float64)
    mean_class_sums = np.zeros(n_distance_class, dtype=np.float64)
    sample_class_sums = np.zeros((n_samples, n_distance_class), dtype=np.float64)

    max_sampled_edges = max(int(config["scatter_sample_size"]), int(config["pairwise_sample_size"]))
    rng = np.random.default_rng(int(config["random_seed"]) + time_index)
    sample_source_idx, sample_sink_idx = sampled_edge_indices(
        n_source=n_source,
        n_sink=n_sink,
        n_edges=max_sampled_edges,
        rng=rng,
        include_diagonal=bool(config["include_diagonal"]),
    )
    sampled_single = np.full(max_sampled_edges, np.nan, dtype=np.float32)
    sampled_bootstrap = np.full((n_samples, max_sampled_edges), np.nan, dtype=np.float32)

    block_count = 0
    for source_slice, sink_slice in iter_matrix_blocks(n_source, n_sink, int(config["block_size"])):
        block_count += 1
        if block_count % 10 == 0:
            print(f"    processed {block_count} blocks")

        single_block = clean_array(
            single_da.isel(time=time_index, source=source_slice, sink=sink_slice, sample=0).values
        )
        bootstrap_block = clean_array(
            bootstrap_da.isel(time=time_index, source=source_slice, sink=sink_slice)
            .transpose("sample", "source", "sink")
            .values
        )
        positive_bootstrap = np.where(bootstrap_block > 0, bootstrap_block, 0.0)

        single_matrix[source_slice, sink_slice] = single_block

        mean_block = np.nanmean(bootstrap_block, axis=0, dtype=np.float64).astype(np.float32)
        std_block = np.nanstd(bootstrap_block, axis=0, ddof=1).astype(np.float32)
        bootstrap_mean[source_slice, sink_slice] = mean_block
        bootstrap_std[source_slice, sink_slice] = std_block

        if config["compute_quantiles"]:
            quantiles = np.nanquantile(bootstrap_block, config["quantiles"], axis=0).astype(np.float32)
            q_low[source_slice, sink_slice] = quantiles[0]
            q_med[source_slice, sink_slice] = quantiles[1]
            q_high[source_slice, sink_slice] = quantiles[2]

        row_sums_boot[:, source_slice] += np.nansum(bootstrap_block, axis=2, dtype=np.float64)
        col_sums_boot[:, sink_slice] += np.nansum(bootstrap_block, axis=1, dtype=np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            wlogw = np.where(positive_bootstrap > 0, positive_bootstrap * np.log(positive_bootstrap), 0.0)
        row_wlogw_boot[:, source_slice] += np.sum(wlogw, axis=2, dtype=np.float64)
        col_wlogw_boot[:, sink_slice] += np.sum(wlogw, axis=1, dtype=np.float64)
        row_max_boot[:, source_slice] = np.maximum(row_max_boot[:, source_slice], np.max(positive_bootstrap, axis=2))
        col_max_boot[:, sink_slice] = np.maximum(col_max_boot[:, sink_slice], np.max(positive_bootstrap, axis=1))
        for idx, threshold in enumerate(thresholds):
            connectance_counts[:, idx] += np.sum(bootstrap_block > threshold, axis=(1, 2))
        for threshold, counts in threshold_presence_counts.items():
            counts[source_slice, sink_slice] += np.sum(bootstrap_block > threshold, axis=0).astype(np.uint8)

        source_anchor_block = community_anchor_onehot[source_slice]
        sink_anchor_block = community_anchor_onehot[sink_slice]
        source_region_block = region_onehot[source_slice]
        sink_region_block = region_onehot[sink_slice]

        single_out_anchor[source_slice] += single_block @ sink_anchor_block
        mean_out_anchor[source_slice] += mean_block @ sink_anchor_block
        sample_out_anchor[:, source_slice, :] += np.einsum("aij,jk->aik", bootstrap_block, sink_anchor_block, optimize=True)
        single_in_anchor[sink_slice] += single_block.T @ source_anchor_block
        mean_in_anchor[sink_slice] += mean_block.T @ source_anchor_block
        sample_in_anchor[:, sink_slice, :] += np.einsum("aij,ik->ajk", bootstrap_block, source_anchor_block, optimize=True)

        single_region_matrix += source_region_block.T @ single_block @ sink_region_block
        mean_region_matrix += source_region_block.T @ mean_block @ sink_region_block
        sample_region_matrices += np.einsum(
            "aij,ik,jl->akl",
            bootstrap_block,
            source_region_block,
            sink_region_block,
            optimize=True,
        )

        if distance_matrix is not None and n_distance_curve > 0:
            distance_block = distance_matrix[source_slice, sink_slice]
            for bin_index in range(n_distance_curve):
                low = distance_curve_edges[bin_index]
                high = distance_curve_edges[bin_index + 1]
                if bin_index == n_distance_curve - 1:
                    mask = (distance_block >= low) & (distance_block <= high)
                else:
                    mask = (distance_block >= low) & (distance_block < high)
                if not np.any(mask):
                    continue
                curve_counts[bin_index] += int(np.sum(mask))
                single_curve_sums[bin_index] += float(np.nansum(single_block[mask], dtype=np.float64))
                mean_curve_sums[bin_index] += float(np.nansum(mean_block[mask], dtype=np.float64))
                sample_curve_sums[:, bin_index] += np.nansum(bootstrap_block[:, mask], axis=1, dtype=np.float64)
                single_curve_positive[bin_index] += int(np.sum(single_block[mask] > 0))
                mean_curve_positive[bin_index] += int(np.sum(mean_block[mask] > 0))
                sample_curve_positive[:, bin_index] += np.sum(bootstrap_block[:, mask] > 0, axis=1, dtype=np.int64)

            for class_index in range(n_distance_class):
                low = distance_class_edges[class_index]
                high = distance_class_edges[class_index + 1]
                if class_index == n_distance_class - 1:
                    class_mask = (distance_block >= low) & (distance_block <= high)
                else:
                    class_mask = (distance_block >= low) & (distance_block < high)
                if not np.any(class_mask):
                    continue
                single_class_sums[class_index] += float(np.nansum(single_block[class_mask], dtype=np.float64))
                mean_class_sums[class_index] += float(np.nansum(mean_block[class_mask], dtype=np.float64))
                sample_class_sums[:, class_index] += np.nansum(bootstrap_block[:, class_mask], axis=1, dtype=np.float64)

        diagonal_overlap = _diagonal_overlap(source_slice, sink_slice)
        if diagonal_overlap is not None:
            source_idx, sink_idx = diagonal_overlap
            diagonal_boot += np.nansum(bootstrap_block[:, source_idx, sink_idx], axis=1, dtype=np.float64)

        global_block_indices = (
            (np.arange(int(source_slice.start), int(source_slice.stop))[:, None] * n_sink)
            + np.arange(int(sink_slice.start), int(sink_slice.stop))[None, :]
        ).astype(np.int64).ravel()
        sample_top_values, sample_top_indices = sample_topk_tracker_update(
            current_values=sample_top_values,
            current_indices=sample_top_indices,
            block_values=bootstrap_block,
            global_block_indices=global_block_indices,
            keep_k=max_edge_top_k,
        )

        assign_sampled_edges_from_block(
            source_idx=sample_source_idx,
            sink_idx=sample_sink_idx,
            source_slice=source_slice,
            sink_slice=sink_slice,
            single_block=single_block,
            bootstrap_block=bootstrap_block,
            sampled_single=sampled_single,
            sampled_bootstrap=sampled_bootstrap,
        )

    bootstrap_cv = coefficient_of_variation(
        bootstrap_std,
        bootstrap_mean,
        epsilon=float(config["cv_epsilon"]),
    )

    row_sums_single = np.nansum(single_matrix, axis=1, dtype=np.float64)
    row_sums_mean = np.nansum(bootstrap_mean, axis=1, dtype=np.float64)
    col_sums_single = np.nansum(single_matrix, axis=0, dtype=np.float64)
    col_sums_mean = np.nansum(bootstrap_mean, axis=0, dtype=np.float64)
    sample_total_weight = np.sum(row_sums_boot, axis=1)

    row_entropy_single, row_effective_single, row_dominance_single, col_entropy_single, col_effective_single, col_dominance_single = compute_local_node_metrics(
        single_matrix,
        row_sums_single,
        col_sums_single,
    )
    row_entropy_mean, row_effective_mean, row_dominance_mean, col_entropy_mean, col_effective_mean, col_dominance_mean = compute_local_node_metrics(
        bootstrap_mean,
        row_sums_mean,
        col_sums_mean,
    )

    sample_row_entropy = np.full(row_sums_boot.shape, np.nan, dtype=np.float64)
    sample_col_entropy = np.full(col_sums_boot.shape, np.nan, dtype=np.float64)
    sample_row_effective = np.full(row_sums_boot.shape, np.nan, dtype=np.float64)
    sample_col_effective = np.full(col_sums_boot.shape, np.nan, dtype=np.float64)
    sample_row_dominance = np.full(row_sums_boot.shape, np.nan, dtype=np.float64)
    sample_col_dominance = np.full(col_sums_boot.shape, np.nan, dtype=np.float64)
    valid_row = row_sums_boot > 0
    valid_col = col_sums_boot > 0
    sample_row_entropy[valid_row] = np.log(row_sums_boot[valid_row]) - (row_wlogw_boot[valid_row] / row_sums_boot[valid_row])
    sample_col_entropy[valid_col] = np.log(col_sums_boot[valid_col]) - (col_wlogw_boot[valid_col] / col_sums_boot[valid_col])
    sample_row_effective[valid_row] = np.exp(sample_row_entropy[valid_row])
    sample_col_effective[valid_col] = np.exp(sample_col_entropy[valid_col])
    sample_row_dominance[valid_row] = row_max_boot[valid_row] / row_sums_boot[valid_row]
    sample_col_dominance[valid_col] = col_max_boot[valid_col] / col_sums_boot[valid_col]

    matrix_metrics = compute_matrix_comparison_metrics(
        single=single_matrix,
        bootstrap_mean=bootstrap_mean,
        thresholds=thresholds,
        time_index=time_index,
        time_label=time_label,
        row_sums_single=row_sums_single,
        row_sums_mean=row_sums_mean,
        col_sums_single=col_sums_single,
        col_sums_mean=col_sums_mean,
        default_change_threshold=float(config["default_change_threshold"]),
        relative_diff_epsilon=float(config["relative_diff_epsilon"]),
        spearman_max_points=int(config["spearman_max_points"]),
        random_seed=int(config["random_seed"]),
    )

    network_single = network_metric_rows(
        matrix=single_matrix,
        row_sums=row_sums_single,
        col_sums=col_sums_single,
        thresholds=thresholds,
        component_thresholds=config["component_thresholds"],
        time_index=time_index,
        time_label=time_label,
        dataset_label="single",
        max_component_edges=int(config["max_component_edges"]),
    )
    network_bootstrap = network_metric_rows(
        matrix=bootstrap_mean,
        row_sums=row_sums_mean,
        col_sums=col_sums_mean,
        thresholds=thresholds,
        component_thresholds=config["component_thresholds"],
        time_index=time_index,
        time_label=time_label,
        dataset_label="bootstrap_mean",
        max_component_edges=int(config["max_component_edges"]),
    )
    network_metrics_df = pd.concat([network_single, network_bootstrap], ignore_index=True)

    edge_overlap_df, edge_churn_detail_df = compute_edge_churn(
        single_matrix=single_matrix,
        bootstrap_mean=bootstrap_mean,
        top_k_values=edge_top_k_values,
        time_index=time_index,
        time_label=time_label,
        top_n_per_status=int(config["top_n_links"]),
    )

    outgoing_rank_df, changing_sources_df = compute_rank_stability(
        single_values=row_sums_single,
        bootstrap_mean_values=row_sums_mean,
        top_k_values=rank_top_k_values,
        entity_label="source",
        time_index=time_index,
        time_label=time_label,
        top_n_changes=int(config["top_n_rank_changes"]),
    )
    incoming_rank_df, changing_sinks_df = compute_rank_stability(
        single_values=col_sums_single,
        bootstrap_mean_values=col_sums_mean,
        top_k_values=rank_top_k_values,
        entity_label="sink",
        time_index=time_index,
        time_label=time_label,
        top_n_changes=int(config["top_n_rank_changes"]),
    )
    rank_summary_df = pd.concat([outgoing_rank_df, incoming_rank_df], ignore_index=True)
    source_rank_uncertainty_summary_df, source_rank_uncertainty_df = compute_rank_uncertainty(
        single_values=row_sums_single,
        bootstrap_mean_values=row_sums_mean,
        bootstrap_samples=row_sums_boot,
        entity_type="source",
        top_k_values=rank_top_k_values,
        time_index=time_index,
        time_label=time_label,
    )
    sink_rank_uncertainty_summary_df, sink_rank_uncertainty_df = compute_rank_uncertainty(
        single_values=col_sums_single,
        bootstrap_mean_values=col_sums_mean,
        bootstrap_samples=col_sums_boot,
        entity_type="sink",
        top_k_values=rank_top_k_values,
        time_index=time_index,
        time_label=time_label,
    )

    top_sources_df = top_nodes_table(
        single_values=row_sums_single,
        bootstrap_mean_values=row_sums_mean,
        entity_type="source",
        time_index=time_index,
        time_label=time_label,
        top_n=int(config["top_k_nodes"]),
    )
    top_sinks_df = top_nodes_table(
        single_values=col_sums_single,
        bootstrap_mean_values=col_sums_mean,
        entity_type="sink",
        time_index=time_index,
        time_label=time_label,
        top_n=int(config["top_k_nodes"]),
    )

    pairwise_sample_n = int(config["pairwise_sample_size"])
    scatter_sample_n = int(config["scatter_sample_size"])
    sampled_mean = np.nanmean(sampled_bootstrap, axis=0).astype(np.float32)
    similarity_df = summarize_sample_similarity(
        sampled_bootstrap=sampled_bootstrap[:, :pairwise_sample_n],
        sampled_single=sampled_single[:pairwise_sample_n],
        time_index=time_index,
        time_label=time_label,
    )

    edge_ci_df = pd.DataFrame()
    outgoing_ci_summary_df = pd.DataFrame()
    outgoing_ci_detail_df = pd.DataFrame()
    incoming_ci_summary_df = pd.DataFrame()
    incoming_ci_detail_df = pd.DataFrame()
    if config["compute_quantiles"]:
        edge_ci_df = summarize_single_vs_bootstrap_ci(
            single=single_matrix,
            q_low=q_low,
            q_high=q_high,
            thresholds=thresholds,
            time_index=time_index,
            time_label=time_label,
        )
        outgoing_ci_summary_df, outgoing_ci_detail_df = bootstrap_interval_summary(
            single_values=row_sums_single,
            bootstrap_samples=row_sums_boot,
            entity_type="outgoing",
            time_index=time_index,
            time_label=time_label,
            top_n=int(config["top_n_rank_changes"]),
        )
        incoming_ci_summary_df, incoming_ci_detail_df = bootstrap_interval_summary(
            single_values=col_sums_single,
            bootstrap_samples=col_sums_boot,
            entity_type="incoming",
            time_index=time_index,
            time_label=time_label,
            top_n=int(config["top_n_rank_changes"]),
        )

    changing_links_df = top_changing_links(
        single=single_matrix,
        bootstrap_mean=bootstrap_mean,
        bootstrap_std=bootstrap_std,
        bootstrap_cv=bootstrap_cv,
        q_low=q_low,
        q_high=q_high,
        time_index=time_index,
        time_label=time_label,
        top_n=int(config["top_n_links"]),
        relative_diff_epsilon=float(config["relative_diff_epsilon"]),
    )
    profile_summary_df, profile_detail_df = compute_profile_divergence(
        single_matrix=single_matrix,
        bootstrap_mean=bootstrap_mean,
        row_sums_single=row_sums_single,
        row_sums_mean=row_sums_mean,
        col_sums_single=col_sums_single,
        col_sums_mean=col_sums_mean,
        time_index=time_index,
        time_label=time_label,
    )

    concentration_summary_base_df, lorenz_curve_df, bootstrap_sample_tail_df = compute_concentration_and_tail(
        single_matrix=single_matrix,
        bootstrap_mean=bootstrap_mean,
        sampled_bootstrap=sampled_bootstrap[:, :pairwise_sample_n],
        sample_top_values=sample_top_values,
        time_index=time_index,
        time_label=time_label,
        top_share_fracs=config["top_share_fractions"],
        top_k_values=edge_top_k_values,
    )
    bootstrap_concentration_row = {
        "time_index": time_index,
        "time_label": time_label,
        "dataset": "bootstrap_samples_sampled",
        "gini": float(np.nanmedian(bootstrap_sample_tail_df["sampled_edge_gini"])),
        "max_edge": float(np.nanmedian(bootstrap_sample_tail_df["max_edge"])),
        "overall_mean": float(np.nanmean(sample_total_weight / (n_source * n_sink))),
        "overall_median": np.nan,
        "positive_edge_fraction": np.nan,
        "positive_median": float(np.nanmedian(bootstrap_sample_tail_df["sampled_positive_median"])),
        "positive_p90": np.nan,
        "positive_p95": float(np.nanmedian(bootstrap_sample_tail_df["sampled_positive_p95"])),
        "positive_p99": np.nan,
    }
    for fraction in config["top_share_fractions"]:
        bootstrap_concentration_row[f"share_top_{int(fraction * 100)}pct"] = np.nan
    for top_k in edge_top_k_values:
        bootstrap_concentration_row[f"top_{top_k}_mean"] = float(np.nanmedian(bootstrap_sample_tail_df[f"top_{top_k}_mean"]))
        bootstrap_concentration_row[f"top_{top_k}_sum"] = float(np.nanmedian(bootstrap_sample_tail_df[f"top_{top_k}_mean"])) * top_k
    bootstrap_concentration_row["strongest_to_positive_median"] = (
        bootstrap_concentration_row["max_edge"] / bootstrap_concentration_row["positive_median"]
        if bootstrap_concentration_row["positive_median"] and np.isfinite(bootstrap_concentration_row["positive_median"])
        else np.nan
    )
    bootstrap_concentration_row["top_1pct_mean_to_overall_mean"] = np.nan
    concentration_summary_df = pd.concat(
        [concentration_summary_base_df, pd.DataFrame([bootstrap_concentration_row])],
        ignore_index=True,
    )
    tail_summary_df = concentration_summary_df.copy()

    outgoing_strength_uncertainty_df = compute_node_strength_uncertainty(
        single_values=row_sums_single,
        bootstrap_mean_values=row_sums_mean,
        bootstrap_samples=row_sums_boot,
        entity_type="source",
        time_index=time_index,
        time_label=time_label,
    )
    incoming_strength_uncertainty_df = compute_node_strength_uncertainty(
        single_values=col_sums_single,
        bootstrap_mean_values=col_sums_mean,
        bootstrap_samples=col_sums_boot,
        entity_type="sink",
        time_index=time_index,
        time_label=time_label,
    )
    node_strength_uncertainty_df = pd.concat(
        [outgoing_strength_uncertainty_df, incoming_strength_uncertainty_df],
        ignore_index=True,
    )

    local_node_metrics_df = summarize_local_node_metrics(
        single_row_metrics=(row_entropy_single, row_effective_single, row_dominance_single),
        mean_row_metrics=(row_entropy_mean, row_effective_mean, row_dominance_mean),
        bootstrap_row_metrics=(sample_row_entropy, sample_row_effective, sample_row_dominance),
        single_col_metrics=(col_entropy_single, col_effective_single, col_dominance_single),
        mean_col_metrics=(col_entropy_mean, col_effective_mean, col_dominance_mean),
        bootstrap_col_metrics=(sample_col_entropy, sample_col_effective, sample_col_dominance),
        time_index=time_index,
        time_label=time_label,
    )

    single_group_features = np.column_stack(
        [
            np.log1p(row_sums_single),
            np.log1p(col_sums_single),
            row_entropy_single,
            col_entropy_single,
            row_dominance_single,
            col_dominance_single,
        ]
    )
    mean_group_features = np.column_stack(
        [
            np.log1p(row_sums_mean),
            np.log1p(col_sums_mean),
            row_entropy_mean,
            col_entropy_mean,
            row_dominance_mean,
            col_dominance_mean,
        ]
    )
    sample_group_features = np.stack(
        [
            np.column_stack(
                [
                    np.log1p(row_sums_boot[sample_idx]),
                    np.log1p(col_sums_boot[sample_idx]),
                    sample_row_entropy[sample_idx],
                    sample_col_entropy[sample_idx],
                    sample_row_dominance[sample_idx],
                    sample_col_dominance[sample_idx],
                ]
            )
            for sample_idx in range(n_samples)
        ]
    )
    grouping_summary_df, grouping_size_df = compute_grouping_summary(
        single_features=single_group_features,
        mean_features=mean_group_features,
        sample_features=sample_group_features,
        n_groups=int(config["grouping_n_clusters"]),
        time_index=time_index,
        time_label=time_label,
        random_seed=int(config["random_seed"]),
    )

    community_summary_df, community_assignments_df, community_coassignment_df, community_sample_similarity_df, community_diagnostics = compute_community_structure(
        single_out_anchor=single_out_anchor,
        single_in_anchor=single_in_anchor,
        mean_out_anchor=mean_out_anchor,
        mean_in_anchor=mean_in_anchor,
        sample_out_anchor=sample_out_anchor,
        sample_in_anchor=sample_in_anchor,
        row_sums_single=row_sums_single,
        row_sums_mean=row_sums_mean,
        col_sums_single=col_sums_single,
        col_sums_mean=col_sums_mean,
        time_index=time_index,
        time_label=time_label,
        reef_metadata=metadata_context.get("reef_metadata"),
        config=config,
    )
    bridge_summary_df, bridge_time_summary_df, bridge_diagnostics = compute_bridge_reefs(
        single_out_anchor=single_out_anchor,
        single_in_anchor=single_in_anchor,
        mean_out_anchor=mean_out_anchor,
        mean_in_anchor=mean_in_anchor,
        sample_out_anchor=sample_out_anchor,
        sample_in_anchor=sample_in_anchor,
        row_sums_single=row_sums_single,
        row_sums_mean=row_sums_mean,
        col_sums_single=col_sums_single,
        col_sums_mean=col_sums_mean,
        community_assignments_df=community_assignments_df,
        time_index=time_index,
        time_label=time_label,
        reef_metadata=metadata_context.get("reef_metadata"),
        config=config,
    )
    distance_curve_df, distance_class_df, distance_diagnostics = compute_distance_uncertainty(
        curve_edges=distance_curve_edges,
        curve_counts=curve_counts,
        single_curve_sums=single_curve_sums,
        mean_curve_sums=mean_curve_sums,
        sample_curve_sums=sample_curve_sums,
        single_curve_positive=single_curve_positive,
        mean_curve_positive=mean_curve_positive,
        sample_curve_positive=sample_curve_positive,
        class_labels=distance_class_labels,
        single_class_sums=single_class_sums,
        mean_class_sums=mean_class_sums,
        sample_class_sums=sample_class_sums,
        time_index=time_index,
        time_label=time_label,
    )
    regional_exchange_df, regional_role_df, regional_summary_df = compute_regional_exchange(
        single_region_matrix=single_region_matrix,
        mean_region_matrix=mean_region_matrix,
        sample_region_matrices=sample_region_matrices,
        region_labels=list(metadata_context["region_labels"]),
        time_index=time_index,
        time_label=time_label,
    )

    stable_edge_summary_df, stable_edge_detail_df = compute_stable_edge_summary(
        single_matrix=single_matrix,
        bootstrap_mean=bootstrap_mean,
        sample_top_indices=sample_top_indices,
        threshold_presence_counts=threshold_presence_counts,
        top_k_values=edge_top_k_values,
        stable_frequency_rules=config["stable_frequency_rules"],
        time_index=time_index,
        time_label=time_label,
        top_n_details=int(config["top_n_links"]),
    )
    backbone_summary_df = compute_backbone_summary(
        single_matrix=single_matrix,
        bootstrap_mean=bootstrap_mean,
        row_sums_single=row_sums_single,
        row_sums_mean=row_sums_mean,
        thresholds=stable_thresholds,
        local_share_thresholds=config["local_share_thresholds"],
        time_index=time_index,
        time_label=time_label,
    )
    mean_bias_df = compute_mean_bias(
        row_sums_single=row_sums_single,
        row_sums_mean=row_sums_mean,
        row_sums_boot=row_sums_boot,
        col_sums_single=col_sums_single,
        col_sums_mean=col_sums_mean,
        col_sums_boot=col_sums_boot,
        bridge_diagnostics=bridge_diagnostics,
        community_diagnostics=community_diagnostics,
        sample_top_indices=sample_top_indices,
        single_matrix=single_matrix,
        bootstrap_mean=bootstrap_mean,
        single_region_matrix=single_region_matrix,
        mean_region_matrix=mean_region_matrix,
        sample_region_matrices=sample_region_matrices,
        distance_diagnostics=distance_diagnostics,
        time_index=time_index,
        time_label=time_label,
        config=config,
    )
    spatial_summary_df, spatial_hotspots_df = compute_spatial_hotspots(
        source_rank_df=source_rank_uncertainty_df,
        sink_rank_df=sink_rank_uncertainty_df,
        bridge_df=bridge_summary_df,
        community_assignments_df=community_assignments_df,
        reef_metadata=metadata_context.get("reef_metadata"),
        time_index=time_index,
        time_label=time_label,
        config=config,
    )
    discrepancy_summary_df = compute_discrepancy_vs_variability(
        sampled_single=sampled_single[:pairwise_sample_n],
        sampled_mean=sampled_mean[:pairwise_sample_n],
        sampled_bootstrap=sampled_bootstrap[:, :pairwise_sample_n],
        row_sums_single=row_sums_single,
        row_sums_mean=row_sums_mean,
        row_sums_boot=row_sums_boot,
        col_sums_single=col_sums_single,
        col_sums_mean=col_sums_mean,
        col_sums_boot=col_sums_boot,
        sample_top_indices=sample_top_indices,
        single_matrix=single_matrix,
        bootstrap_mean=bootstrap_mean,
        top_k_values=edge_top_k_values,
        time_index=time_index,
        time_label=time_label,
    )

    threshold_rows = []
    n_valid_edges = int(np.isfinite(single_matrix).sum())
    for idx, threshold in enumerate(thresholds):
        bootstrap_connectance = connectance_counts[:, idx] / n_valid_edges if n_valid_edges else np.nan
        single_connectance = float(np.mean(single_matrix[np.isfinite(single_matrix)] > threshold)) if n_valid_edges else np.nan
        mean_connectance = float(np.mean(bootstrap_mean[np.isfinite(bootstrap_mean)] > threshold)) if n_valid_edges else np.nan
        threshold_rows.append(
            {
                "time_index": time_index,
                "time_label": time_label,
                "threshold": threshold,
                "bootstrap_sample_connectance_mean": float(np.nanmean(bootstrap_connectance)),
                "bootstrap_sample_connectance_std": float(np.nanstd(bootstrap_connectance, ddof=1)),
                "bootstrap_sample_connectance_q025": float(np.nanquantile(bootstrap_connectance, 0.025)),
                "bootstrap_sample_connectance_q975": float(np.nanquantile(bootstrap_connectance, 0.975)),
                "single_connectance": single_connectance,
                "bootstrap_mean_connectance": mean_connectance,
                "single_self_recruitment": float(np.nansum(np.diag(single_matrix))),
                "bootstrap_mean_self_recruitment": float(np.nansum(np.diag(bootstrap_mean))),
                "bootstrap_sample_self_recruitment_mean": float(np.nanmean(diagonal_boot)),
                "bootstrap_sample_self_recruitment_q025": float(np.nanquantile(diagonal_boot, 0.025)),
                "bootstrap_sample_self_recruitment_q975": float(np.nanquantile(diagonal_boot, 0.975)),
            }
        )
    threshold_summary_df = pd.DataFrame(threshold_rows)

    mean_positive = bootstrap_mean[np.isfinite(bootstrap_mean) & (bootstrap_mean > 0)]
    std_positive = bootstrap_std[np.isfinite(bootstrap_std) & (bootstrap_mean > 0)]
    cv_positive = bootstrap_cv[np.isfinite(bootstrap_cv) & (bootstrap_mean > 0)]
    uncertainty_summary = pd.DataFrame(
        [
            {
                "time_index": time_index,
                "time_label": time_label,
                "mean_edge_sd_nonzero": float(np.nanmean(std_positive)) if std_positive.size else np.nan,
                "median_edge_sd_nonzero": float(np.nanmedian(std_positive)) if std_positive.size else np.nan,
                "median_edge_cv_nonzero": float(np.nanmedian(cv_positive)) if cv_positive.size else np.nan,
                "p95_edge_cv_nonzero": float(np.nanquantile(cv_positive, 0.95)) if cv_positive.size else np.nan,
                "sampled_single_vs_bootstrap_mean_r": float(similarity_df["single_vs_bootstrap_mean_r"].iloc[0])
                if not similarity_df.empty and "single_vs_bootstrap_mean_r" in similarity_df
                else np.nan,
                "mean_pairwise_bootstrap_r": float(similarity_df["mean_pairwise_bootstrap_r"].iloc[0])
                if not similarity_df.empty and "mean_pairwise_bootstrap_r" in similarity_df
                else np.nan,
                "fraction_edge_values_outside_ci": float(edge_ci_df.loc[edge_ci_df["subset"] == "all", "fraction_outside_ci"].iloc[0])
                if not edge_ci_df.empty and (edge_ci_df["subset"] == "all").any()
                else np.nan,
                "fraction_outgoing_outside_ci": float(outgoing_ci_summary_df["fraction_outside_ci"].iloc[0])
                if not outgoing_ci_summary_df.empty
                else np.nan,
                "fraction_incoming_outside_ci": float(incoming_ci_summary_df["fraction_outside_ci"].iloc[0])
                if not incoming_ci_summary_df.empty
                else np.nan,
            }
        ]
    )

    ordination_coords_df = pd.DataFrame()
    ordination_meta_df = pd.DataFrame()
    if config["ordination_enabled"]:
        bootstrap_features = np.hstack([row_sums_boot, col_sums_boot])
        single_features = np.concatenate([row_sums_single, col_sums_single])
        ordination_coords_df, ordination_meta_df = pca_ordination(
            bootstrap_features=bootstrap_features,
            single_features=single_features,
            time_index=time_index,
            time_label=time_label,
        )

    single_metrics_for_coverage = {
        "self_recruitment": float(np.nansum(np.diag(single_matrix))),
        "outgoing_gini": float(gini(row_sums_single)),
        "incoming_gini": float(gini(col_sums_single)),
        "max_edge": float(np.nanmax(single_matrix)),
    }
    bootstrap_mean_metrics_for_coverage = {
        "self_recruitment": float(np.nansum(np.diag(bootstrap_mean))),
        "outgoing_gini": float(gini(row_sums_mean)),
        "incoming_gini": float(gini(col_sums_mean)),
        "max_edge": float(np.nanmax(bootstrap_mean)),
    }
    bootstrap_sample_metrics_for_coverage = {
        "self_recruitment": diagonal_boot,
        "outgoing_gini": np.array([gini(values) for values in row_sums_boot], dtype=np.float64),
        "incoming_gini": np.array([gini(values) for values in col_sums_boot], dtype=np.float64),
        "max_edge": sample_top_values[:, 0],
    }
    for idx, threshold in enumerate(thresholds):
        metric_name = f"connectance_{threshold:g}"
        single_metrics_for_coverage[metric_name] = float(np.mean(single_matrix[np.isfinite(single_matrix)] > threshold))
        bootstrap_mean_metrics_for_coverage[metric_name] = float(
            np.mean(bootstrap_mean[np.isfinite(bootstrap_mean)] > threshold)
        )
        bootstrap_sample_metrics_for_coverage[metric_name] = connectance_counts[:, idx] / n_valid_edges
    for top_k in edge_top_k_values:
        metric_name = f"top_{top_k}_mean"
        single_metrics_for_coverage[metric_name] = float(
            np.mean(np.sort(single_matrix.ravel())[-top_k:])
        )
        bootstrap_mean_metrics_for_coverage[metric_name] = float(
            np.mean(np.sort(bootstrap_mean.ravel())[-top_k:])
        )
        bootstrap_sample_metrics_for_coverage[metric_name] = np.mean(sample_top_values[:, :top_k], axis=1)

    coverage_summary_df = compute_system_coverage_summary(
        single_metrics=single_metrics_for_coverage,
        bootstrap_mean_metrics=bootstrap_mean_metrics_for_coverage,
        bootstrap_sample_metrics=bootstrap_sample_metrics_for_coverage,
        time_index=time_index,
        time_label=time_label,
    )

    source_top_n = min(int(config["top_k_nodes"]), row_sums_boot.shape[1])
    sink_top_n = min(int(config["top_k_nodes"]), col_sums_boot.shape[1])
    bridge_top_n = min(int(config["bridge_top_k"]), bridge_diagnostics["sample_bridge_score"].shape[1])
    sample_balance_matrix = np.sum(sample_region_matrices, axis=2) - np.sum(sample_region_matrices, axis=1)
    variance_input_df = pd.DataFrame(
        [
            {
                "time_index": time_index,
                "time_label": time_label,
                "metric_name": "source_top_mean_strength",
                "mean_value": float(np.mean(np.sort(row_sums_mean)[-source_top_n:])),
                "within_bootstrap_var": float(
                    np.nanvar(np.mean(np.sort(row_sums_boot, axis=1)[:, -source_top_n:], axis=1), ddof=1)
                ),
            },
            {
                "time_index": time_index,
                "time_label": time_label,
                "metric_name": "sink_top_mean_strength",
                "mean_value": float(np.mean(np.sort(col_sums_mean)[-sink_top_n:])),
                "within_bootstrap_var": float(
                    np.nanvar(np.mean(np.sort(col_sums_boot, axis=1)[:, -sink_top_n:], axis=1), ddof=1)
                ),
            },
            {
                "time_index": time_index,
                "time_label": time_label,
                "metric_name": "source_topk_membership",
                "mean_value": float(np.nanmean(source_rank_uncertainty_df[f"prob_top_{rank_top_k_values[0]}"])),
                "within_bootstrap_var": float(
                    np.nanmean(
                        source_rank_uncertainty_df[f"prob_top_{rank_top_k_values[0]}"]
                        * (1.0 - source_rank_uncertainty_df[f"prob_top_{rank_top_k_values[0]}"])
                    )
                ),
            },
            {
                "time_index": time_index,
                "time_label": time_label,
                "metric_name": "bridge_top_mean_score",
                "mean_value": float(np.mean(np.sort(bridge_diagnostics["mean_bridge_score"])[-bridge_top_n:])),
                "within_bootstrap_var": float(
                    np.nanvar(
                        np.mean(np.sort(bridge_diagnostics["sample_bridge_score"], axis=1)[:, -bridge_top_n:], axis=1),
                        ddof=1,
                    )
                ),
            },
            {
                "time_index": time_index,
                "time_label": time_label,
                "metric_name": "community_stability",
                "mean_value": float(np.nanmedian(community_diagnostics["sample_vs_consensus_ari"])),
                "within_bootstrap_var": float(np.nanvar(community_diagnostics["sample_vs_consensus_ari"], ddof=1)),
            },
            {
                "time_index": time_index,
                "time_label": time_label,
                "metric_name": "long_distance_fraction",
                "mean_value": float(distance_diagnostics["mean_long_fraction"]),
                "within_bootstrap_var": float(np.nanvar(distance_diagnostics["sample_long_fraction"], ddof=1)),
            },
            {
                "time_index": time_index,
                "time_label": time_label,
                "metric_name": "regional_balance_sd",
                "mean_value": float(np.nanstd(np.sum(mean_region_matrix, axis=1) - np.sum(mean_region_matrix, axis=0), ddof=1)),
                "within_bootstrap_var": float(np.nanvar(np.nanstd(sample_balance_matrix, axis=1, ddof=1), ddof=1)),
            },
        ]
    )

    source_rank_comparison_df = rank_summary_df[rank_summary_df["entity_type"] == "source"].reset_index(drop=True)
    sink_rank_comparison_df = rank_summary_df[rank_summary_df["entity_type"] == "sink"].reset_index(drop=True)
    source_profile_changes_df = profile_detail_df[profile_detail_df["entity_type"] == "source"].reset_index(drop=True)
    sink_profile_changes_df = profile_detail_df[profile_detail_df["entity_type"] == "sink"].reset_index(drop=True)
    top_uncertain_nodes_df = node_strength_uncertainty_df.sort_values(
        ["single_outside_ci", "bootstrap_sample_cv", "single_z_score"],
        ascending=[False, False, False],
    ).head(int(config["top_n_rank_changes"])).reset_index(drop=True)

    intermediate_path = npz_path(Path(config["intermediate_dir"]), time_index, "time_summary")
    if config["save_intermediate_npz"]:
        single_flat = single_matrix.ravel()
        mean_flat = bootstrap_mean.ravel()
        single_top_weights = np.sort(single_flat[np.isfinite(single_flat)])[-max_edge_top_k:]
        mean_top_weights = np.sort(mean_flat[np.isfinite(mean_flat)])[-max_edge_top_k:]
        np.savez_compressed(
            intermediate_path,
            time_index=np.array([time_index], dtype=np.int32),
            time_label=np.array([time_label]),
            single=block_average(single_matrix, output_size=int(config["heatmap_size"])),
            bootstrap_mean=block_average(bootstrap_mean, output_size=int(config["heatmap_size"])),
            difference=block_average(single_matrix - bootstrap_mean, output_size=int(config["heatmap_size"])),
            bootstrap_sd=block_average(bootstrap_std, output_size=int(config["heatmap_size"])),
            bootstrap_cv=block_average(bootstrap_cv, output_size=int(config["heatmap_size"])),
            scatter_single=sampled_single[:scatter_sample_n],
            scatter_mean=sampled_mean[:scatter_sample_n],
            row_sums_single=row_sums_single.astype(np.float32),
            row_sums_mean=row_sums_mean.astype(np.float32),
            col_sums_single=col_sums_single.astype(np.float32),
            col_sums_mean=col_sums_mean.astype(np.float32),
            single_top_edge_weights=single_top_weights.astype(np.float32),
            bootstrap_mean_top_edge_weights=mean_top_weights.astype(np.float32),
            community_coassignment=community_diagnostics["coassignment_block"].astype(np.float32),
        )

    del single_matrix, bootstrap_mean, bootstrap_std, bootstrap_cv, row_sums_boot, col_sums_boot
    del sampled_bootstrap, sampled_single, threshold_presence_counts, sample_top_values, sample_top_indices
    del sample_out_anchor, sample_in_anchor, sample_region_matrices, sample_curve_sums, sample_curve_positive, sample_class_sums
    gc.collect()

    node_ci_summary_df = pd.concat([outgoing_ci_summary_df, incoming_ci_summary_df], ignore_index=True)
    node_ci_detail_df = pd.concat([outgoing_ci_detail_df, incoming_ci_detail_df], ignore_index=True)
    top_nodes_df = pd.concat([top_sources_df, top_sinks_df], ignore_index=True)

    return {
        "matrix_metrics": matrix_metrics,
        "network_metrics": network_metrics_df,
        "edge_overlap": edge_overlap_df,
        "edge_churn": edge_churn_detail_df,
        "rank_summary": rank_summary_df,
        "source_rank_comparison": source_rank_comparison_df,
        "sink_rank_comparison": sink_rank_comparison_df,
        "source_rank_uncertainty_summary": source_rank_uncertainty_summary_df,
        "sink_rank_uncertainty_summary": sink_rank_uncertainty_summary_df,
        "source_rank_uncertainty": source_rank_uncertainty_df,
        "sink_rank_uncertainty": sink_rank_uncertainty_df,
        "changing_sources": changing_sources_df,
        "changing_sinks": changing_sinks_df,
        "changing_links": changing_links_df,
        "profile_divergence_summary": profile_summary_df,
        "source_profile_changes": source_profile_changes_df,
        "sink_profile_changes": sink_profile_changes_df,
        "grouping_summary": grouping_summary_df,
        "grouping_sizes": grouping_size_df,
        "community_summary": community_summary_df,
        "community_assignments": community_assignments_df,
        "community_coassignment": community_coassignment_df,
        "community_sample_similarity": community_sample_similarity_df,
        "bridge_summary": bridge_summary_df,
        "bridge_time_summary": bridge_time_summary_df,
        "distance_curve_summary": distance_curve_df,
        "distance_class_summary": distance_class_df,
        "regional_exchange": regional_exchange_df,
        "regional_roles": regional_role_df,
        "regional_summary": regional_summary_df,
        "mean_bias": mean_bias_df,
        "spatial_summary": spatial_summary_df,
        "spatial_hotspots": spatial_hotspots_df,
        "variance_input": variance_input_df,
        "concentration_summary": concentration_summary_df,
        "lorenz_curve_points": lorenz_curve_df,
        "tail_summary": tail_summary_df,
        "bootstrap_sample_tail": bootstrap_sample_tail_df,
        "node_strength_uncertainty": node_strength_uncertainty_df,
        "top_uncertain_nodes": top_uncertain_nodes_df,
        "stable_edge_summary": stable_edge_summary_df,
        "stable_edge_detail": stable_edge_detail_df,
        "local_node_metrics": local_node_metrics_df,
        "backbone_summary": backbone_summary_df,
        "coverage_summary": coverage_summary_df,
        "discrepancy_summary": discrepancy_summary_df,
        "top_nodes": top_nodes_df,
        "edge_ci": edge_ci_df,
        "node_ci_summary": node_ci_summary_df,
        "node_ci_detail": node_ci_detail_df,
        "threshold_summary": threshold_summary_df,
        "uncertainty_summary": uncertainty_summary,
        "similarity": similarity_df,
        "ordination_coords": ordination_coords_df,
        "ordination_meta": ordination_meta_df,
        "intermediate_path": intermediate_path,
        "time_label": time_label,
    }


def select_times_for_figures(comparison_df: pd.DataFrame, max_times: int = 3) -> list[int]:
    """Choose representative times spanning low, medium, and high difference."""
    comparison_all = comparison_df[comparison_df["subset"] == "all"].sort_values("frobenius_norm")
    if comparison_all.empty:
        return []
    candidate_positions = np.linspace(0, len(comparison_all) - 1, num=min(max_times, len(comparison_all)), dtype=int)
    selected = comparison_all.iloc[candidate_positions]["time_index"].tolist()
    return list(dict.fromkeys(int(value) for value in selected))


def concat_nonempty_frames(*frames: pd.DataFrame) -> pd.DataFrame:
    """Concatenate only non-empty dataframes."""
    valid_frames = [frame for frame in frames if isinstance(frame, pd.DataFrame) and not frame.empty]
    if not valid_frames:
        return pd.DataFrame()
    return pd.concat(valid_frames, ignore_index=True)


def compute_ecological_threshold_suite(
    bootstrap_da,
    single_da,
    time_indices: list[int],
    all_time_values: np.ndarray,
    metadata_context: dict[str, Any],
    config: dict,
) -> pd.DataFrame:
    """Recompute selected time steps for structural threshold robustness diagnostics."""
    if not time_indices:
        return pd.DataFrame()

    n_source = int(single_da.sizes["source"])
    n_sink = int(single_da.sizes["sink"])
    thresholds = [0.0, 1e-8]
    frames = []
    for time_index in time_indices:
        time_label = format_time_label(all_time_values[time_index], time_index)
        print(f"[threshold] analysing {time_label}")
        bootstrap_mean = np.full((n_source, n_sink), np.nan, dtype=np.float32)
        threshold_presence_counts = {
            threshold: np.zeros((n_source, n_sink), dtype=np.uint8) for threshold in thresholds
        }
        for source_slice, sink_slice in iter_matrix_blocks(n_source, n_sink, int(config["block_size"])):
            bootstrap_block = clean_array(
                bootstrap_da.isel(time=time_index, source=source_slice, sink=sink_slice)
                .transpose("sample", "source", "sink")
                .values
            )
            mean_block = np.nanmean(bootstrap_block, axis=0, dtype=np.float64).astype(np.float32)
            bootstrap_mean[source_slice, sink_slice] = mean_block
            for threshold, counts in threshold_presence_counts.items():
                counts[source_slice, sink_slice] += np.sum(bootstrap_block > threshold, axis=0).astype(np.uint8)
        frames.append(
            compute_ecological_thresholds(
                bootstrap_mean=bootstrap_mean,
                threshold_presence_counts=threshold_presence_counts,
                anchor_onehot=np.asarray(metadata_context["community_anchor_onehot"], dtype=np.float32),
                time_index=time_index,
                time_label=time_label,
                config=config,
            )
        )
        del bootstrap_mean, threshold_presence_counts
        gc.collect()
    return concat_nonempty_frames(*frames)


def generate_figures(
    comparison_df: pd.DataFrame,
    edge_ci_df: pd.DataFrame,
    node_ci_df: pd.DataFrame,
    similarity_df: pd.DataFrame,
    network_df: pd.DataFrame,
    rank_uncertainty_summary_df: pd.DataFrame,
    source_changes_df: pd.DataFrame,
    sink_changes_df: pd.DataFrame,
    source_profile_df: pd.DataFrame,
    sink_profile_df: pd.DataFrame,
    grouping_summary_df: pd.DataFrame,
    grouping_size_df: pd.DataFrame,
    community_summary_df: pd.DataFrame,
    community_assignments_df: pd.DataFrame,
    bridge_summary_df: pd.DataFrame,
    bridge_time_summary_df: pd.DataFrame,
    distance_curve_df: pd.DataFrame,
    distance_class_df: pd.DataFrame,
    regional_exchange_df: pd.DataFrame,
    regional_role_df: pd.DataFrame,
    mean_bias_df: pd.DataFrame,
    spatial_summary_df: pd.DataFrame,
    spatial_hotspots_df: pd.DataFrame,
    variance_summary_df: pd.DataFrame,
    family_similarity_df: pd.DataFrame,
    family_overlap_df: pd.DataFrame,
    family_ordination_df: pd.DataFrame,
    family_uncertainty_df: pd.DataFrame,
    ecological_threshold_df: pd.DataFrame,
    concentration_df: pd.DataFrame,
    lorenz_df: pd.DataFrame,
    tail_df: pd.DataFrame,
    node_strength_uncertainty_df: pd.DataFrame,
    edge_overlap_df: pd.DataFrame,
    edge_churn_df: pd.DataFrame,
    stable_edge_summary_df: pd.DataFrame,
    local_node_metrics_df: pd.DataFrame,
    backbone_summary_df: pd.DataFrame,
    coverage_summary_df: pd.DataFrame,
    discrepancy_summary_df: pd.DataFrame,
    source_rank_uncertainty_df: pd.DataFrame,
    sink_rank_uncertainty_df: pd.DataFrame,
    ordination_coords_df: pd.DataFrame,
    ordination_meta_df: pd.DataFrame,
    config: dict,
) -> list[int]:
    """Write all summary figures to disk."""
    cache_dir = Path(config["output_dir"]) / ".cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))

    from src.plotting import (
        plot_bridge_map,
        plot_bridge_removal_impact,
        plot_bridge_strength_scatter,
        plot_bridge_uncertainty,
        plot_coassignment_heatmap,
        plot_community_consensus_map,
        plot_community_similarity_timeseries,
        plot_community_stability_map,
        plot_distribution_panel,
        plot_distance_decay,
        plot_distance_uncertainty,
        plot_discrepancy_vs_variability,
        plot_dominance_ratio,
        plot_ecological_threshold_robustness,
        plot_edge_churn,
        plot_edge_overlap_timeseries,
        plot_effective_links_distribution,
        plot_entropy_strength,
        plot_family_ordination,
        plot_family_overlap,
        plot_family_similarity_heatmap,
        plot_family_uncertainty,
        plot_group_sizes,
        plot_grouping_similarity,
        plot_heatmap_panel,
        plot_concentration_timeseries,
        plot_coverage_examples,
        plot_coverage_rates,
        plot_long_distance_panel,
        plot_lorenz_curves,
        plot_mean_bias_summary,
        plot_metric_timeseries,
        plot_backbone_overlap,
        plot_node_cv_hist,
        plot_node_strength_uncertainty,
        plot_ordination,
        plot_profile_divergence_hist,
        plot_profile_top_changes,
        plot_rank_change_bars,
        plot_rank_correlation_timeseries,
        plot_rank_panel,
        plot_scatter,
        plot_stable_edge_counts,
        plot_stable_edge_overlap,
        plot_spatial_clustering,
        plot_spatial_hotspots,
        plot_tail_timeseries,
        plot_top_edge_strength_distribution,
        plot_topk_probability,
        plot_rank_uncertainty_intervals,
        plot_regional_exchange_heatmap,
        plot_regional_pair_intervals,
        plot_regional_roles,
        plot_threshold_sensitivity,
        plot_uncertainty_timeseries,
        plot_variance_decomposition,
    )

    figures_dir = Path(config["figures_dir"])
    plot_metric_timeseries(comparison_df, figures_dir / "time_series_matrix_metrics.png")
    plot_threshold_sensitivity(
        comparison_df=comparison_df,
        ci_df=edge_ci_df,
        network_df=network_df,
        path=figures_dir / "threshold_sensitivity.png",
    )
    plot_uncertainty_timeseries(
        edge_ci_df=edge_ci_df,
        node_ci_df=node_ci_df,
        similarity_df=similarity_df,
        path=figures_dir / "uncertainty_time_series.png",
    )
    plot_rank_correlation_timeseries(rank_uncertainty_summary_df, figures_dir / "rank_correlation_time_series.png")
    plot_grouping_similarity(grouping_summary_df, figures_dir / "grouping_similarity_time_series.png")
    plot_community_similarity_timeseries(community_summary_df, figures_dir / "community_similarity_time_series.png")
    plot_concentration_timeseries(concentration_df, figures_dir / "concentration_time_series.png")
    plot_tail_timeseries(tail_df, figures_dir / "tail_metric_time_series.png")
    plot_edge_overlap_timeseries(edge_overlap_df, figures_dir / "topk_overlap_time_series.png")
    plot_stable_edge_counts(stable_edge_summary_df, figures_dir / "stable_edge_counts.png")
    plot_stable_edge_overlap(stable_edge_summary_df, figures_dir / "stable_edge_overlap.png")
    plot_backbone_overlap(backbone_summary_df, figures_dir / "backbone_overlap_time_series.png")
    plot_coverage_rates(coverage_summary_df, figures_dir / "coverage_rates_time_series.png")
    plot_discrepancy_vs_variability(discrepancy_summary_df, figures_dir / "single_vs_bootstrap_variability.png")
    plot_long_distance_panel(distance_class_df, figures_dir / "long_distance_fraction_time_series.png")
    plot_mean_bias_summary(mean_bias_df, figures_dir / "mean_bias_summary.png")
    plot_spatial_clustering(spatial_summary_df, figures_dir / "spatial_clustering_summary.png")
    plot_variance_decomposition(variance_summary_df, figures_dir / "variance_decomposition.png")
    plot_family_similarity_heatmap(family_similarity_df, figures_dir / "family_similarity_heatmap.png")
    plot_family_overlap(family_overlap_df, figures_dir / "family_overlap_summary.png")
    plot_family_ordination(family_ordination_df, figures_dir / "family_ordination.png")
    plot_family_uncertainty(family_uncertainty_df, figures_dir / "family_uncertainty_panel.png")
    plot_ecological_threshold_robustness(ecological_threshold_df, figures_dir / "ecological_threshold_robustness.png")

    selected_times = select_times_for_figures(
        comparison_df=comparison_df,
        max_times=int(config["n_selected_times_for_figures"]),
    )
    for time_index in selected_times:
        archive = np.load(npz_path(Path(config["intermediate_dir"]), time_index, "time_summary"))
        time_label = str(archive["time_label"][0])
        matrices = {
            "single": archive["single"],
            "bootstrap_mean": archive["bootstrap_mean"],
            "difference": archive["difference"],
            "bootstrap_sd": archive["bootstrap_sd"],
            "bootstrap_cv": archive["bootstrap_cv"],
        }
        plot_heatmap_panel(
            matrices=matrices,
            time_label=time_label,
            path=figures_dir / f"heatmap_panel_time_{time_index:02d}.png",
        )
        plot_scatter(
            single_values=archive["scatter_single"],
            mean_values=archive["scatter_mean"],
            time_label=time_label,
            path=figures_dir / f"scatter_time_{time_index:02d}.png",
        )
        plot_distribution_panel(
            single_sampled_edges=archive["scatter_single"],
            mean_sampled_edges=archive["scatter_mean"],
            row_sums_single=archive["row_sums_single"],
            row_sums_mean=archive["row_sums_mean"],
            col_sums_single=archive["col_sums_single"],
            col_sums_mean=archive["col_sums_mean"],
            time_label=time_label,
            path=figures_dir / f"distributions_time_{time_index:02d}.png",
        )
        plot_rank_panel(
            row_sums_single=archive["row_sums_single"],
            row_sums_mean=archive["row_sums_mean"],
            col_sums_single=archive["col_sums_single"],
            col_sums_mean=archive["col_sums_mean"],
            time_label=time_label,
            path=figures_dir / f"rank_panel_time_{time_index:02d}.png",
        )
        plot_top_edge_strength_distribution(
            single_weights=archive["single_top_edge_weights"],
            mean_weights=archive["bootstrap_mean_top_edge_weights"],
            time_label=time_label,
            path=figures_dir / f"top_edge_strength_distribution_time_{time_index:02d}.png",
        )
        plot_rank_change_bars(
            source_changes_df=source_changes_df[source_changes_df["time_index"] == time_index],
            sink_changes_df=sink_changes_df[sink_changes_df["time_index"] == time_index],
            time_label=time_label,
            path=figures_dir / f"rank_change_bars_time_{time_index:02d}.png",
        )
        plot_profile_divergence_hist(
            profile_df=concat_nonempty_frames(
                source_profile_df[source_profile_df["time_index"] == time_index],
                sink_profile_df[sink_profile_df["time_index"] == time_index],
            ),
            time_label=time_label,
            path=figures_dir / f"profile_divergence_hist_time_{time_index:02d}.png",
        )
        plot_profile_top_changes(
            source_df=source_profile_df[source_profile_df["time_index"] == time_index],
            sink_df=sink_profile_df[sink_profile_df["time_index"] == time_index],
            time_label=time_label,
            path=figures_dir / f"profile_top_changes_time_{time_index:02d}.png",
        )
        plot_group_sizes(
            grouping_size_df=grouping_size_df,
            time_label=time_label,
            path=figures_dir / f"group_size_summary_time_{time_index:02d}.png",
        )
        plot_lorenz_curves(
            lorenz_df=lorenz_df,
            time_label=time_label,
            path=figures_dir / f"lorenz_curves_time_{time_index:02d}.png",
        )
        plot_node_strength_uncertainty(
            node_df=node_strength_uncertainty_df,
            time_label=time_label,
            path=figures_dir / f"node_strength_uncertainty_time_{time_index:02d}.png",
        )
        plot_node_cv_hist(
            node_df=node_strength_uncertainty_df,
            time_label=time_label,
            path=figures_dir / f"node_strength_cv_hist_time_{time_index:02d}.png",
        )
        plot_edge_churn(
            edge_churn_df=edge_churn_df,
            time_label=time_label,
            path=figures_dir / f"edge_churn_time_{time_index:02d}.png",
        )
        plot_entropy_strength(
            local_node_metrics_df=local_node_metrics_df,
            time_label=time_label,
            path=figures_dir / f"entropy_strength_time_{time_index:02d}.png",
        )
        plot_dominance_ratio(
            local_node_metrics_df=local_node_metrics_df,
            time_label=time_label,
            path=figures_dir / f"dominance_ratio_time_{time_index:02d}.png",
        )
        plot_effective_links_distribution(
            local_node_metrics_df=local_node_metrics_df,
            time_label=time_label,
            path=figures_dir / f"effective_links_distribution_time_{time_index:02d}.png",
        )
        plot_coverage_examples(
            node_df=node_strength_uncertainty_df,
            time_label=time_label,
            path=figures_dir / f"coverage_examples_time_{time_index:02d}.png",
        )
        rank_detail = concat_nonempty_frames(
            source_rank_uncertainty_df[source_rank_uncertainty_df["time_index"] == time_index],
            sink_rank_uncertainty_df[sink_rank_uncertainty_df["time_index"] == time_index],
        )
        plot_topk_probability(
            rank_df=rank_detail,
            time_label=time_label,
            top_k=int(config["rank_top_k_values"][0]),
            path=figures_dir / f"topk_probability_time_{time_index:02d}.png",
        )
        plot_rank_uncertainty_intervals(
            rank_df=rank_detail,
            time_label=time_label,
            path=figures_dir / f"rank_uncertainty_intervals_time_{time_index:02d}.png",
        )
        if "community_coassignment" in archive:
            plot_coassignment_heatmap(
                coassignment_block=archive["community_coassignment"],
                time_label=time_label,
                path=figures_dir / f"community_coassignment_time_{time_index:02d}.png",
            )
        plot_community_consensus_map(
            assignments_df=community_assignments_df,
            time_label=time_label,
            path=figures_dir / f"community_consensus_time_{time_index:02d}.png",
        )
        plot_community_stability_map(
            assignments_df=community_assignments_df,
            time_label=time_label,
            path=figures_dir / f"community_stability_time_{time_index:02d}.png",
        )
        plot_bridge_map(
            bridge_df=bridge_summary_df,
            time_label=time_label,
            path=figures_dir / f"bridge_map_time_{time_index:02d}.png",
        )
        plot_bridge_uncertainty(
            bridge_df=bridge_summary_df,
            time_label=time_label,
            path=figures_dir / f"bridge_uncertainty_time_{time_index:02d}.png",
        )
        plot_bridge_strength_scatter(
            bridge_df=bridge_summary_df,
            time_label=time_label,
            path=figures_dir / f"bridge_strength_scatter_time_{time_index:02d}.png",
        )
        plot_bridge_removal_impact(
            bridge_df=bridge_summary_df,
            time_label=time_label,
            path=figures_dir / f"bridge_removal_impact_time_{time_index:02d}.png",
        )
        plot_distance_decay(
            distance_curve_df=distance_curve_df,
            time_label=time_label,
            path=figures_dir / f"distance_decay_time_{time_index:02d}.png",
        )
        plot_distance_uncertainty(
            distance_curve_df=distance_curve_df,
            time_label=time_label,
            path=figures_dir / f"distance_uncertainty_time_{time_index:02d}.png",
        )
        plot_regional_exchange_heatmap(
            exchange_df=regional_exchange_df,
            time_label=time_label,
            path=figures_dir / f"regional_exchange_time_{time_index:02d}.png",
        )
        plot_regional_roles(
            role_df=regional_role_df,
            time_label=time_label,
            path=figures_dir / f"regional_roles_time_{time_index:02d}.png",
        )
        plot_regional_pair_intervals(
            exchange_df=regional_exchange_df,
            time_label=time_label,
            path=figures_dir / f"regional_pair_intervals_time_{time_index:02d}.png",
        )
        plot_spatial_hotspots(
            spatial_hotspot_df=spatial_hotspots_df,
            time_label=time_label,
            path=figures_dir / f"spatial_hotspots_time_{time_index:02d}.png",
        )

        coords = ordination_coords_df[ordination_coords_df["time_index"] == time_index]
        meta = ordination_meta_df[ordination_meta_df["time_index"] == time_index]
        plot_ordination(
            coords_df=coords,
            meta_df=meta,
            time_label=time_label,
            path=figures_dir / f"ordination_time_{time_index:02d}.png",
        )
    return selected_times


def main() -> int:
    args = parse_args()
    config = build_runtime_config(args)
    ensure_output_directories(config)

    print("[load] opening datasets")
    bootstrap_ds = open_connectivity_dataset(Path(config["bootstrap_path"]))
    single_ds = open_connectivity_dataset(Path(config["single_path"]))
    bootstrap_da = prepare_connectivity_array(
        bootstrap_ds,
        var_name=str(config["connectivity_var"]),
        treatment_index=int(config["treatment_index"]),
    )
    single_da = prepare_connectivity_array(
        single_ds,
        var_name=str(config["connectivity_var"]),
        treatment_index=int(config["treatment_index"]),
    )

    info = validate_datasets(bootstrap_da, single_da)
    config["n_bootstrap_samples"] = int(info["n_bootstrap_samples"])
    bootstrap_summary = dataset_summary("bootstrap", Path(config["bootstrap_path"]), bootstrap_da)
    single_summary = dataset_summary("single", Path(config["single_path"]), single_da)
    print_summary(bootstrap_summary)
    print_summary(single_summary)
    print(f"[check] n_time={info['n_time']}, n_source={info['n_source']}, n_sink={info['n_sink']}, n_bootstrap_samples={info['n_bootstrap_samples']}")
    for warning in info.get("time_warnings", []):
        print(f"[warning] {warning}")

    metadata_context = load_reef_metadata(
        path=Path(config["reef_metadata_path"]) if config.get("reef_metadata_path") else None,
        n_nodes=int(info["n_source"]),
        community_anchor_column=str(config["community_anchor_column"]),
        region_column=str(config["region_column"]),
    )
    for note in metadata_context.get("notes", []):
        print(f"[warning] {note}")

    distance_matrix, distance_diagnostics = infer_distance_matrix(
        bootstrap_ds,
        variable_candidates=list(config["distance_var_candidates"]),
    )
    if distance_diagnostics.get("selected_variable"):
        print(f"[info] using '{distance_diagnostics['selected_variable']}' as the continuous distance matrix")
    if distance_diagnostics.get("note"):
        print(f"[warning] {distance_diagnostics['note']}")

    all_time_values = resolve_time_values(bootstrap_da, single_da)
    primary_family_name = _family_name_from_single_path(Path(config["single_path"]))
    if config["time_indices"] is None:
        time_indices = list(range(info["n_time"]))
    else:
        time_indices = list(config["time_indices"])
        invalid = [idx for idx in time_indices if idx < 0 or idx >= info["n_time"]]
        if invalid:
            raise ValueError(f"Invalid time indices requested: {invalid}")

    tables: dict[str, list[pd.DataFrame]] = {
        "matrix_metrics": [],
        "network_metrics": [],
        "edge_overlap": [],
        "edge_churn": [],
        "rank_summary": [],
        "source_rank_comparison": [],
        "sink_rank_comparison": [],
        "source_rank_uncertainty_summary": [],
        "sink_rank_uncertainty_summary": [],
        "source_rank_uncertainty": [],
        "sink_rank_uncertainty": [],
        "changing_sources": [],
        "changing_sinks": [],
        "changing_links": [],
        "profile_divergence_summary": [],
        "source_profile_changes": [],
        "sink_profile_changes": [],
        "grouping_summary": [],
        "grouping_sizes": [],
        "community_summary": [],
        "community_assignments": [],
        "community_coassignment": [],
        "community_sample_similarity": [],
        "bridge_summary": [],
        "bridge_time_summary": [],
        "distance_curve_summary": [],
        "distance_class_summary": [],
        "regional_exchange": [],
        "regional_roles": [],
        "regional_summary": [],
        "mean_bias": [],
        "spatial_summary": [],
        "spatial_hotspots": [],
        "variance_input": [],
        "concentration_summary": [],
        "lorenz_curve_points": [],
        "tail_summary": [],
        "bootstrap_sample_tail": [],
        "node_strength_uncertainty": [],
        "top_uncertain_nodes": [],
        "stable_edge_summary": [],
        "stable_edge_detail": [],
        "local_node_metrics": [],
        "backbone_summary": [],
        "coverage_summary": [],
        "discrepancy_summary": [],
        "top_nodes": [],
        "edge_ci": [],
        "node_ci_summary": [],
        "node_ci_detail": [],
        "threshold_summary": [],
        "uncertainty_summary": [],
        "similarity": [],
        "ordination_coords": [],
        "ordination_meta": [],
    }

    for time_index in time_indices:
        time_label = format_time_label(all_time_values[time_index], time_index)
        print(f"[time {time_index:02d}] analysing {time_label}")
        outputs = analyse_time_step(
            bootstrap_da=bootstrap_da,
            single_da=single_da,
            time_index=time_index,
            time_label=time_label,
            config=config,
            metadata_context=metadata_context,
            distance_matrix=distance_matrix,
        )
        for key in tables:
            frame = outputs[key]
            if isinstance(frame, pd.DataFrame) and not frame.empty:
                tables[key].append(frame)

    combined = {key: pd.concat(value, ignore_index=True) if value else pd.DataFrame() for key, value in tables.items()}

    family_specs = discover_family_datasets(
        dataset_dir=Path(config["family_dataset_dir"]),
        primary_bootstrap_path=Path(config["bootstrap_path"]),
        primary_single_path=Path(config["single_path"]),
    )
    family_summary_df, family_similarity_df, family_overlap_df, family_ordination_df, family_uncertainty_df = build_family_comparison(
        family_specs=family_specs,
        metadata_context=metadata_context,
        config=config,
        primary_family_name=primary_family_name,
        primary_uncertainty_df=combined["uncertainty_summary"],
        primary_bridge_time_df=combined["bridge_time_summary"],
        primary_community_df=combined["community_summary"],
        primary_distance_class_df=combined["distance_class_summary"],
    )

    variance_inputs_df = _time_components(combined["variance_input"])
    if not variance_inputs_df.empty:
        variance_inputs_df["family"] = primary_family_name
    family_variance_rows = []
    if not family_summary_df.empty:
        for row in family_summary_df.itertuples(index=False):
            for metric_name, value in (
                ("source_top_mean_strength", row.source_top_mean_strength),
                ("sink_top_mean_strength", row.sink_top_mean_strength),
                ("bridge_top_mean_score", row.bridge_top_mean_score),
                ("long_distance_fraction", row.long_distance_fraction_mean),
            ):
                family_variance_rows.append(
                    {
                        "family": row.family,
                        "time_index": np.nan,
                        "time_label": "family_mean",
                        "year": np.nan,
                        "spawning_period": "family_mean",
                        "metric_name": metric_name,
                        "mean_value": value,
                        "within_bootstrap_var": np.nan,
                    }
                )
    variance_all_df = pd.concat([variance_inputs_df, pd.DataFrame(family_variance_rows)], ignore_index=True) if family_variance_rows else variance_inputs_df
    variance_summary_df, variance_detail_df = compute_variance_decomposition(variance_all_df)

    threshold_sensitivity_summary = combined["threshold_summary"].merge(
        combined["matrix_metrics"][combined["matrix_metrics"]["subset"] == "threshold"][
            ["time_index", "time_label", "threshold", "pearson_r", "rmse", "mae", "changed_fraction"]
        ],
        on=["time_index", "time_label", "threshold"],
        how="left",
    )
    if not combined["edge_ci"].empty:
        threshold_sensitivity_summary = threshold_sensitivity_summary.merge(
            combined["edge_ci"][combined["edge_ci"]["subset"] == "threshold"][
                ["time_index", "time_label", "threshold", "fraction_outside_ci"]
            ],
            on=["time_index", "time_label", "threshold"],
            how="left",
        )

    tables_dir = Path(config["tables_dir"])
    write_dataframe(combined["matrix_metrics"], tables_dir / "per_time_matrix_comparison.csv")
    write_dataframe(combined["network_metrics"], tables_dir / "per_time_network_metrics.csv")
    write_dataframe(combined["edge_overlap"], tables_dir / "top_link_overlap.csv")
    write_dataframe(combined["edge_churn"], tables_dir / "top_link_churn_detail.csv")
    write_dataframe(combined["rank_summary"], tables_dir / "rank_stability_summary.csv")
    write_dataframe(combined["source_rank_comparison"], tables_dir / "source_rank_comparison.csv")
    write_dataframe(combined["sink_rank_comparison"], tables_dir / "sink_rank_comparison.csv")
    write_dataframe(combined["source_rank_uncertainty_summary"], tables_dir / "source_rank_uncertainty_summary.csv")
    write_dataframe(combined["sink_rank_uncertainty_summary"], tables_dir / "sink_rank_uncertainty_summary.csv")
    write_dataframe(combined["source_rank_uncertainty"], tables_dir / "source_rank_uncertainty.csv")
    write_dataframe(combined["sink_rank_uncertainty"], tables_dir / "sink_rank_uncertainty.csv")
    write_dataframe(combined["changing_sources"], tables_dir / "top_changing_sources.csv")
    write_dataframe(combined["changing_sinks"], tables_dir / "top_changing_sinks.csv")
    write_dataframe(combined["changing_links"], tables_dir / "top_changing_links.csv")
    write_dataframe(combined["profile_divergence_summary"], tables_dir / "profile_divergence_summary.csv")
    write_dataframe(combined["source_profile_changes"], tables_dir / "source_profile_changes.csv")
    write_dataframe(combined["sink_profile_changes"], tables_dir / "sink_profile_changes.csv")
    write_dataframe(combined["grouping_summary"], tables_dir / "grouping_similarity_summary.csv")
    write_dataframe(combined["grouping_sizes"], tables_dir / "grouping_size_summary.csv")
    write_dataframe(combined["community_summary"], tables_dir / "community_structure_summary.csv")
    write_dataframe(combined["community_assignments"], tables_dir / "community_assignments.csv")
    write_dataframe(combined["community_coassignment"], tables_dir / "community_coassignment_summary.csv")
    write_dataframe(combined["community_sample_similarity"], tables_dir / "community_bootstrap_similarity.csv")
    write_dataframe(combined["bridge_summary"], tables_dir / "bridge_reef_summary.csv")
    write_dataframe(combined["bridge_time_summary"], tables_dir / "bridge_reef_time_summary.csv")
    write_dataframe(combined["distance_curve_summary"], tables_dir / "distance_uncertainty_curve.csv")
    write_dataframe(combined["distance_class_summary"], tables_dir / "distance_uncertainty_classes.csv")
    write_dataframe(combined["regional_exchange"], tables_dir / "regional_exchange_summary.csv")
    write_dataframe(combined["regional_roles"], tables_dir / "regional_role_summary.csv")
    write_dataframe(combined["regional_summary"], tables_dir / "regional_exchange_time_summary.csv")
    write_dataframe(combined["mean_bias"], tables_dir / "mean_bias_summary.csv")
    write_dataframe(combined["spatial_summary"], tables_dir / "spatial_clustering_summary.csv")
    write_dataframe(combined["spatial_hotspots"], tables_dir / "spatial_hotspots.csv")
    write_dataframe(combined["variance_input"], tables_dir / "variance_inputs.csv")
    write_dataframe(variance_summary_df, tables_dir / "variance_decomposition_summary.csv")
    write_dataframe(variance_detail_df, tables_dir / "variance_decomposition_detail.csv")
    write_dataframe(family_summary_df, tables_dir / "family_summary.csv")
    write_dataframe(family_similarity_df, tables_dir / "family_similarity_summary.csv")
    write_dataframe(family_overlap_df, tables_dir / "family_topk_overlap.csv")
    write_dataframe(family_ordination_df, tables_dir / "family_ordination.csv")
    write_dataframe(family_uncertainty_df, tables_dir / "family_uncertainty_summary.csv")
    write_dataframe(combined["concentration_summary"], tables_dir / "concentration_metrics_summary.csv")
    write_dataframe(combined["lorenz_curve_points"], tables_dir / "lorenz_curve_points.csv")
    write_dataframe(combined["tail_summary"], tables_dir / "tail_metrics_summary.csv")
    write_dataframe(combined["bootstrap_sample_tail"], tables_dir / "bootstrap_sample_tail_summary.csv")
    write_dataframe(combined["node_strength_uncertainty"], tables_dir / "node_strength_uncertainty.csv")
    write_dataframe(combined["top_uncertain_nodes"], tables_dir / "top_uncertain_nodes.csv")
    write_dataframe(combined["stable_edge_summary"], tables_dir / "stable_edge_summary.csv")
    write_dataframe(combined["stable_edge_detail"], tables_dir / "stable_edge_detail.csv")
    write_dataframe(combined["local_node_metrics"], tables_dir / "local_node_metrics.csv")
    write_dataframe(combined["backbone_summary"], tables_dir / "backbone_summary.csv")
    write_dataframe(combined["coverage_summary"], tables_dir / "coverage_summary.csv")
    write_dataframe(combined["discrepancy_summary"], tables_dir / "discrepancy_vs_variability.csv")
    write_dataframe(combined["top_nodes"], tables_dir / "top_sources_and_sinks.csv")
    write_dataframe(combined["edge_ci"], tables_dir / "edge_ci_summary.csv")
    write_dataframe(combined["node_ci_summary"], tables_dir / "node_ci_summary.csv")
    write_dataframe(combined["node_ci_detail"], tables_dir / "node_ci_detail.csv")
    write_dataframe(combined["threshold_summary"], tables_dir / "bootstrap_connectance_summary.csv")
    write_dataframe(threshold_sensitivity_summary, tables_dir / "threshold_sensitivity_summary.csv")
    write_dataframe(combined["uncertainty_summary"], tables_dir / "bootstrap_uncertainty_summary.csv")
    write_dataframe(combined["similarity"], tables_dir / "sample_similarity_summary.csv")
    write_dataframe(combined["ordination_coords"], tables_dir / "ordination_coordinates.csv")
    write_dataframe(combined["ordination_meta"], tables_dir / "ordination_summary.csv")

    selected_time_indices = select_times_for_figures(
        comparison_df=combined["matrix_metrics"],
        max_times=int(config["n_selected_times_for_figures"]),
    )
    ecological_threshold_df = compute_ecological_threshold_suite(
        bootstrap_da=bootstrap_da,
        single_da=single_da,
        time_indices=selected_time_indices[: int(config["n_selected_times_for_threshold_robustness"])],
        all_time_values=all_time_values,
        metadata_context=metadata_context,
        config=config,
    )
    write_dataframe(ecological_threshold_df, tables_dir / "ecological_threshold_summary.csv")

    report_text = build_report(
        comparison_df=combined["matrix_metrics"],
        overlap_df=combined["edge_overlap"],
        uncertainty_df=combined["edge_ci"],
        node_ci_df=combined["node_ci_summary"],
        similarity_df=combined["similarity"],
        concentration_df=combined["concentration_summary"],
        profile_df=combined["profile_divergence_summary"],
        grouping_df=combined["grouping_summary"],
        coverage_df=combined["coverage_summary"],
        discrepancy_df=combined["discrepancy_summary"],
        rank_uncertainty_df=concat_nonempty_frames(
            combined["source_rank_uncertainty_summary"],
            combined["sink_rank_uncertainty_summary"],
        ),
        community_df=combined["community_summary"],
        bridge_df=combined["bridge_time_summary"],
        distance_df=combined["distance_class_summary"],
        regional_df=combined["regional_summary"],
        mean_bias_df=combined["mean_bias"],
        family_df=family_summary_df,
        variance_df=variance_summary_df,
        threshold_df=ecological_threshold_df,
        config=config,
    )
    write_markdown(report_text, Path(config["output_dir"]) / "report_summary.md")
    if not config["skip_plots"]:
        print("[plot] writing figures")
        selected_time_indices = generate_figures(
            comparison_df=combined["matrix_metrics"],
            edge_ci_df=combined["edge_ci"],
            node_ci_df=combined["node_ci_summary"],
            similarity_df=combined["similarity"],
            network_df=combined["network_metrics"],
            rank_uncertainty_summary_df=concat_nonempty_frames(
                combined["source_rank_uncertainty_summary"],
                combined["sink_rank_uncertainty_summary"],
            ),
            source_changes_df=combined["changing_sources"],
            sink_changes_df=combined["changing_sinks"],
            source_profile_df=combined["source_profile_changes"],
            sink_profile_df=combined["sink_profile_changes"],
            grouping_summary_df=combined["grouping_summary"],
            grouping_size_df=combined["grouping_sizes"],
            community_summary_df=combined["community_summary"],
            community_assignments_df=combined["community_assignments"],
            bridge_summary_df=combined["bridge_summary"],
            bridge_time_summary_df=combined["bridge_time_summary"],
            distance_curve_df=combined["distance_curve_summary"],
            distance_class_df=combined["distance_class_summary"],
            regional_exchange_df=combined["regional_exchange"],
            regional_role_df=combined["regional_roles"],
            mean_bias_df=combined["mean_bias"],
            spatial_summary_df=combined["spatial_summary"],
            spatial_hotspots_df=combined["spatial_hotspots"],
            variance_summary_df=variance_summary_df,
            family_similarity_df=family_similarity_df,
            family_overlap_df=family_overlap_df,
            family_ordination_df=family_ordination_df,
            family_uncertainty_df=family_uncertainty_df,
            ecological_threshold_df=ecological_threshold_df,
            concentration_df=combined["concentration_summary"],
            lorenz_df=combined["lorenz_curve_points"],
            tail_df=combined["tail_summary"],
            node_strength_uncertainty_df=combined["node_strength_uncertainty"],
            edge_overlap_df=combined["edge_overlap"],
            edge_churn_df=combined["edge_churn"],
            stable_edge_summary_df=combined["stable_edge_summary"],
            local_node_metrics_df=combined["local_node_metrics"],
            backbone_summary_df=combined["backbone_summary"],
            coverage_summary_df=combined["coverage_summary"],
            discrepancy_summary_df=combined["discrepancy_summary"],
            source_rank_uncertainty_df=combined["source_rank_uncertainty"],
            sink_rank_uncertainty_df=combined["sink_rank_uncertainty"],
            ordination_coords_df=combined["ordination_coords"],
            ordination_meta_df=combined["ordination_meta"],
            config=config,
        )
    analysis_catalog = build_analysis_catalog(selected_time_indices)
    write_markdown(analysis_catalog, Path(config["output_dir"]) / "analysis_catalog.md")

    bootstrap_ds.close()
    single_ds.close()
    print("[done] analysis complete")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover
        print(f"[error] {exc}", file=sys.stderr)
        raise
