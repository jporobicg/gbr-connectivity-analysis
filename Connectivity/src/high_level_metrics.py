"""Higher-level ecological interpretation analyses for connectivity comparison."""

from __future__ import annotations

from typing import Any, Iterable
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.vq import kmeans2
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree

from .extended_metrics import adjusted_rand_index, topk_membership_from_samples
from .metadata import append_node_metadata
from .preprocessing import block_average


def _standardize_features(features: np.ndarray) -> np.ndarray:
    mean = np.nanmean(features, axis=0, keepdims=True)
    std = np.nanstd(features, axis=0, keepdims=True)
    std[std == 0] = 1.0
    return np.nan_to_num((features - mean) / std, nan=0.0)


def _safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2:
        return np.nan
    if float(np.nanstd(x)) == 0.0 or float(np.nanstd(y)) == 0.0:
        return np.nan
    result = stats.spearmanr(x, y, nan_policy="omit")
    return float(result.statistic)


def _safe_pearson(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2:
        return np.nan
    if float(np.nanstd(x)) == 0.0 or float(np.nanstd(y)) == 0.0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def build_directed_profile_features(
    outgoing_anchor: np.ndarray,
    incoming_anchor: np.ndarray,
    outgoing_total: np.ndarray,
    incoming_total: np.ndarray,
) -> np.ndarray:
    """Create a directed flow-profile feature matrix for node clustering."""
    outgoing_profile = np.zeros_like(outgoing_anchor, dtype=np.float64)
    incoming_profile = np.zeros_like(incoming_anchor, dtype=np.float64)
    valid_out = outgoing_total > 0
    valid_in = incoming_total > 0
    outgoing_profile[valid_out] = outgoing_anchor[valid_out] / outgoing_total[valid_out, None]
    incoming_profile[valid_in] = incoming_anchor[valid_in] / incoming_total[valid_in, None]
    return np.column_stack(
        [
            outgoing_profile,
            incoming_profile,
            np.log1p(np.clip(outgoing_total, 0.0, None)),
            np.log1p(np.clip(incoming_total, 0.0, None)),
        ]
    )


def _cluster_labels(features: np.ndarray, n_clusters: int, seed: int) -> np.ndarray:
    if n_clusters <= 1 or features.shape[0] == 0:
        return np.zeros(features.shape[0], dtype=np.int32)
    standardized = _standardize_features(features)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        _, labels = kmeans2(standardized, n_clusters, minit="points", seed=seed)
    return labels.astype(np.int32)


def _calinski_harabasz_score(features: np.ndarray, labels: np.ndarray) -> float:
    unique = np.unique(labels)
    n_points = features.shape[0]
    n_clusters = unique.size
    if n_clusters <= 1 or n_clusters >= n_points:
        return np.nan

    centered = _standardize_features(features)
    overall_mean = np.mean(centered, axis=0)
    between = 0.0
    within = 0.0
    for label in unique:
        cluster = centered[labels == label]
        if cluster.size == 0:
            continue
        cluster_mean = np.mean(cluster, axis=0)
        between += cluster.shape[0] * float(np.sum((cluster_mean - overall_mean) ** 2))
        within += float(np.sum((cluster - cluster_mean) ** 2))
    if within <= 0:
        return np.nan
    return float((between / (n_clusters - 1)) / (within / (n_points - n_clusters)))


def choose_cluster_count(features: np.ndarray, min_clusters: int, max_clusters: int, seed: int) -> int:
    """Choose a practical community count using a simple CH-style criterion."""
    if features.shape[0] == 0:
        return 1
    if np.allclose(np.nanstd(features, axis=0), 0.0):
        return 1

    best_k = max(2, min_clusters)
    best_score = -np.inf
    for n_clusters in range(max(2, min_clusters), max(max_clusters, min_clusters) + 1):
        labels = _cluster_labels(features, n_clusters=n_clusters, seed=seed + n_clusters)
        score = _calinski_harabasz_score(features, labels)
        if np.isfinite(score) and score > best_score:
            best_score = score
            best_k = n_clusters
    return int(best_k)


def align_labels(reference_labels: np.ndarray, candidate_labels: np.ndarray) -> np.ndarray:
    """Align candidate cluster labels to a reference labeling."""
    reference = np.asarray(reference_labels, dtype=np.int32)
    candidate = np.asarray(candidate_labels, dtype=np.int32)
    ref_unique = np.unique(reference)
    cand_unique = np.unique(candidate)
    contingency = np.zeros((cand_unique.size, ref_unique.size), dtype=np.int64)
    for left_idx, left in enumerate(cand_unique):
        for right_idx, right in enumerate(ref_unique):
            contingency[left_idx, right_idx] = int(np.sum((candidate == left) & (reference == right)))

    if contingency.size == 0:
        return candidate

    row_ind, col_ind = linear_sum_assignment(-contingency)
    mapping = {cand_unique[row]: ref_unique[col] for row, col in zip(row_ind, col_ind, strict=True)}
    aligned = np.array([mapping.get(label, label) for label in candidate], dtype=np.int32)
    return aligned


def _value_percentile(sample_values: np.ndarray, value: float) -> float:
    finite = sample_values[np.isfinite(sample_values)]
    if finite.size == 0 or not np.isfinite(value):
        return np.nan
    lower = float(np.mean(finite < value))
    equal = float(np.mean(finite == value))
    return 100.0 * (lower + 0.5 * equal)


def _distribution_rows(
    time_index: int,
    time_label: str,
    analysis_domain: str,
    metric_name: str,
    sample_values: np.ndarray,
    single_value: float,
    mean_value: float,
) -> list[dict[str, Any]]:
    finite = np.asarray(sample_values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    q025 = float(np.nanquantile(finite, 0.025)) if finite.size else np.nan
    q25 = float(np.nanquantile(finite, 0.25)) if finite.size else np.nan
    q50 = float(np.nanmedian(finite)) if finite.size else np.nan
    q75 = float(np.nanquantile(finite, 0.75)) if finite.size else np.nan
    q975 = float(np.nanquantile(finite, 0.975)) if finite.size else np.nan

    rows: list[dict[str, Any]] = []
    for reference_name, value in (("single", single_value), ("bootstrap_mean", mean_value)):
        rows.append(
            {
                "time_index": time_index,
                "time_label": time_label,
                "analysis_domain": analysis_domain,
                "metric_name": metric_name,
                "reference": reference_name,
                "reference_value": float(value),
                "bootstrap_sample_median": q50,
                "bootstrap_q025": q025,
                "bootstrap_q25": q25,
                "bootstrap_q975": q975,
                "reference_percentile": _value_percentile(finite, float(value)),
                "outside_middle_50": bool((value < q25) or (value > q75 if not np.isnan(q75) else False)),
                "outside_95_interval": bool((value < q025) or (value > q975)),
                "reference_minus_sample_median": float(value - q50) if np.isfinite(q50) else np.nan,
            }
        )
    return rows


def _safe_nanmean(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    return float(np.mean(finite)) if finite.size else np.nan


def _safe_nanquantile(values: np.ndarray, quantile: float) -> float:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    return float(np.quantile(finite, quantile)) if finite.size else np.nan


def _top_k_overlap(reference_values: np.ndarray, consensus_set: set[int], top_k: int) -> float:
    if top_k <= 0:
        return np.nan
    top_idx = set(np.argsort(reference_values)[-top_k:])
    return len(top_idx & consensus_set) / float(top_k)


def _top_k_overlap_from_indices(indices: np.ndarray, consensus_set: set[int], top_k: int) -> float:
    if top_k <= 0:
        return np.nan
    return len(set(map(int, indices[:top_k])) & consensus_set) / float(top_k)


def _bridge_components(outgoing_anchor: np.ndarray, incoming_anchor: np.ndarray) -> dict[str, np.ndarray]:
    outgoing_total = np.sum(outgoing_anchor, axis=-1)
    incoming_total = np.sum(incoming_anchor, axis=-1)

    outgoing_profile = np.zeros_like(outgoing_anchor, dtype=np.float64)
    incoming_profile = np.zeros_like(incoming_anchor, dtype=np.float64)
    valid_out = outgoing_total > 0
    valid_in = incoming_total > 0
    outgoing_profile[valid_out] = outgoing_anchor[valid_out] / outgoing_total[valid_out][..., None]
    incoming_profile[valid_in] = incoming_anchor[valid_in] / incoming_total[valid_in][..., None]

    participation_out = 1.0 - np.sum(outgoing_profile**2, axis=-1)
    participation_in = 1.0 - np.sum(incoming_profile**2, axis=-1)
    anchor_switch = 1.0 - np.sum(outgoing_profile * incoming_profile, axis=-1)
    through_flow = np.sqrt(np.clip(outgoing_total, 0.0, None) * np.clip(incoming_total, 0.0, None))
    bridge_score = through_flow * np.sqrt(np.clip(participation_out, 0.0, None) * np.clip(participation_in, 0.0, None)) * np.clip(anchor_switch, 0.0, None)

    total_bridge = np.sum(bridge_score, axis=-1, keepdims=True)
    removal_impact = np.divide(
        bridge_score,
        total_bridge,
        out=np.zeros_like(bridge_score, dtype=np.float64),
        where=total_bridge > 0,
    )

    return {
        "outgoing_total": outgoing_total,
        "incoming_total": incoming_total,
        "participation_out": participation_out,
        "participation_in": participation_in,
        "anchor_switch": anchor_switch,
        "bridge_score": bridge_score,
        "removal_impact": removal_impact,
    }


def compute_community_structure(
    single_out_anchor: np.ndarray,
    single_in_anchor: np.ndarray,
    mean_out_anchor: np.ndarray,
    mean_in_anchor: np.ndarray,
    sample_out_anchor: np.ndarray,
    sample_in_anchor: np.ndarray,
    row_sums_single: np.ndarray,
    row_sums_mean: np.ndarray,
    col_sums_single: np.ndarray,
    col_sums_mean: np.ndarray,
    time_index: int,
    time_label: str,
    reef_metadata: pd.DataFrame | None,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Estimate broad connectivity regions and their uncertainty."""
    single_features = build_directed_profile_features(single_out_anchor, single_in_anchor, row_sums_single, col_sums_single)
    mean_features = build_directed_profile_features(mean_out_anchor, mean_in_anchor, row_sums_mean, col_sums_mean)

    sample_features = np.stack(
        [
            build_directed_profile_features(
                sample_out_anchor[sample_idx],
                sample_in_anchor[sample_idx],
                np.sum(sample_out_anchor[sample_idx], axis=1),
                np.sum(sample_in_anchor[sample_idx], axis=1),
            )
            for sample_idx in range(sample_out_anchor.shape[0])
        ]
    )

    selected_k = choose_cluster_count(
        mean_features,
        min_clusters=int(config["community_min_clusters"]),
        max_clusters=int(config["community_max_clusters"]),
        seed=int(config["random_seed"]) + time_index,
    )

    mean_labels = _cluster_labels(mean_features, selected_k, seed=int(config["random_seed"]) + time_index + 1)
    single_labels = _cluster_labels(single_features, selected_k, seed=int(config["random_seed"]) + time_index + 2)
    sample_labels = np.stack(
        [
            align_labels(
                mean_labels,
                _cluster_labels(sample_features[sample_idx], selected_k, seed=int(config["random_seed"]) + time_index + 10 + sample_idx),
            )
            for sample_idx in range(sample_features.shape[0])
        ]
    )

    consensus_counts = np.zeros((mean_labels.size, selected_k), dtype=np.int32)
    for cluster_id in range(selected_k):
        consensus_counts[:, cluster_id] = np.sum(sample_labels == cluster_id, axis=0)
    consensus_labels = np.argmax(consensus_counts, axis=1).astype(np.int32)
    consensus_probability = np.max(consensus_counts, axis=1) / float(sample_labels.shape[0])
    community_instability = 1.0 - consensus_probability

    sample_vs_mean = np.array([adjusted_rand_index(mean_labels, labels) for labels in sample_labels], dtype=np.float64)
    sample_vs_consensus = np.array([adjusted_rand_index(consensus_labels, labels) for labels in sample_labels], dtype=np.float64)
    single_vs_samples = np.array([adjusted_rand_index(single_labels, labels) for labels in sample_labels], dtype=np.float64)
    occupied_counts = np.array([np.unique(labels).size for labels in sample_labels], dtype=np.float64)

    summary_df = pd.DataFrame(
        [
            {
                "time_index": time_index,
                "time_label": time_label,
                "selected_n_communities": int(selected_k),
                "occupied_communities_single": int(np.unique(single_labels).size),
                "occupied_communities_bootstrap_mean": int(np.unique(mean_labels).size),
                "occupied_communities_consensus": int(np.unique(consensus_labels).size),
                "sample_occupied_communities_median": float(np.nanmedian(occupied_counts)),
                "single_vs_bootstrap_mean_ari": adjusted_rand_index(single_labels, mean_labels),
                "single_vs_consensus_ari": adjusted_rand_index(single_labels, consensus_labels),
                "bootstrap_mean_vs_consensus_ari": adjusted_rand_index(mean_labels, consensus_labels),
                "bootstrap_mean_vs_samples_ari_median": float(np.nanmedian(sample_vs_mean)),
                "bootstrap_mean_vs_samples_ari_q025": float(np.nanquantile(sample_vs_mean, 0.025)),
                "bootstrap_mean_vs_samples_ari_q975": float(np.nanquantile(sample_vs_mean, 0.975)),
                "sample_vs_consensus_ari_median": float(np.nanmedian(sample_vs_consensus)),
                "single_vs_samples_ari_median": float(np.nanmedian(single_vs_samples)),
                "mean_consensus_probability": float(np.nanmean(consensus_probability)),
                "median_consensus_probability": float(np.nanmedian(consensus_probability)),
                "p10_consensus_probability": float(np.nanquantile(consensus_probability, 0.10)),
            }
        ]
    )

    assignments_df = pd.DataFrame(
        {
            "time_index": time_index,
            "time_label": time_label,
            "node_id": np.arange(mean_labels.size, dtype=np.int32),
            "single_community": single_labels,
            "bootstrap_mean_community": mean_labels,
            "consensus_community": consensus_labels,
            "consensus_probability": consensus_probability,
            "community_instability": community_instability,
            "single_matches_consensus": single_labels == consensus_labels,
            "bootstrap_mean_matches_consensus": mean_labels == consensus_labels,
        }
    )
    assignments_df = append_node_metadata(assignments_df, reef_metadata)

    coassignment = np.zeros((mean_labels.size, mean_labels.size), dtype=np.float32)
    for labels in sample_labels:
        coassignment += (labels[:, None] == labels[None, :]).astype(np.float32)
    coassignment /= float(sample_labels.shape[0])

    coassignment_rows: list[dict[str, Any]] = []
    unique_consensus = np.unique(consensus_labels)
    for source_cluster in unique_consensus:
        source_idx = np.flatnonzero(consensus_labels == source_cluster)
        for sink_cluster in unique_consensus:
            sink_idx = np.flatnonzero(consensus_labels == sink_cluster)
            block = coassignment[np.ix_(source_idx, sink_idx)]
            if block.size == 0:
                continue
            if source_cluster == sink_cluster:
                upper = block[np.triu_indices(block.shape[0], k=1)] if block.shape[0] > 1 else block.ravel()
                mean_value = float(np.nanmean(upper)) if upper.size else 1.0
            else:
                mean_value = float(np.nanmean(block))
            coassignment_rows.append(
                {
                    "time_index": time_index,
                    "time_label": time_label,
                    "source_consensus_community": int(source_cluster),
                    "sink_consensus_community": int(sink_cluster),
                    "pair_type": "within" if source_cluster == sink_cluster else "between",
                    "mean_coassignment_probability": mean_value,
                    "source_community_size": int(source_idx.size),
                    "sink_community_size": int(sink_idx.size),
                }
            )
    coassignment_df = pd.DataFrame(coassignment_rows)

    sample_similarity_df = pd.DataFrame(
        {
            "time_index": time_index,
            "time_label": time_label,
            "sample_id": np.arange(sample_labels.shape[0], dtype=np.int32),
            "occupied_communities": occupied_counts,
            "ari_to_bootstrap_mean": sample_vs_mean,
            "ari_to_consensus": sample_vs_consensus,
            "ari_to_single": single_vs_samples,
        }
    )

    order = np.lexsort((community_instability, consensus_labels))
    diagnostics = {
        "selected_k": int(selected_k),
        "single_labels": single_labels,
        "mean_labels": mean_labels,
        "sample_labels": sample_labels,
        "consensus_labels": consensus_labels,
        "consensus_probability": consensus_probability,
        "coassignment_block": block_average(coassignment[np.ix_(order, order)], output_size=int(config["heatmap_size"])),
        "sample_vs_consensus_ari": sample_vs_consensus,
    }
    return summary_df, assignments_df, coassignment_df, sample_similarity_df, diagnostics


def compute_bridge_reefs(
    single_out_anchor: np.ndarray,
    single_in_anchor: np.ndarray,
    mean_out_anchor: np.ndarray,
    mean_in_anchor: np.ndarray,
    sample_out_anchor: np.ndarray,
    sample_in_anchor: np.ndarray,
    row_sums_single: np.ndarray,
    row_sums_mean: np.ndarray,
    col_sums_single: np.ndarray,
    col_sums_mean: np.ndarray,
    community_assignments_df: pd.DataFrame,
    time_index: int,
    time_label: str,
    reef_metadata: pd.DataFrame | None,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Approximate stepping-stone reefs using inter-anchor brokerage and participation."""
    single_components = _bridge_components(single_out_anchor, single_in_anchor)
    mean_components = _bridge_components(mean_out_anchor, mean_in_anchor)
    sample_components = _bridge_components(sample_out_anchor, sample_in_anchor)
    sample_bridge = np.asarray(sample_components["bridge_score"], dtype=np.float64)

    top_k = int(config["bridge_top_k"])
    bridge_topk_probability = topk_membership_from_samples(sample_bridge, [top_k])[top_k]
    q_low = np.nanquantile(sample_bridge, 0.025, axis=0)
    q_high = np.nanquantile(sample_bridge, 0.975, axis=0)
    sample_mean = np.nanmean(sample_bridge, axis=0)
    sample_sd = np.nanstd(sample_bridge, axis=0, ddof=1)
    bridge_cv = np.divide(sample_sd, sample_mean, out=np.full(sample_sd.shape, np.nan), where=np.abs(sample_mean) > 1e-12)

    summary_df = pd.DataFrame(
        {
            "time_index": time_index,
            "time_label": time_label,
            "node_id": np.arange(mean_out_anchor.shape[0], dtype=np.int32),
            "single_bridge_score": single_components["bridge_score"],
            "bootstrap_mean_bridge_score": mean_components["bridge_score"],
            "bootstrap_bridge_mean": sample_mean,
            "bootstrap_bridge_sd": sample_sd,
            "bootstrap_bridge_cv": bridge_cv,
            "bootstrap_bridge_q025": q_low,
            "bootstrap_bridge_q975": q_high,
            f"prob_top_bridge_{top_k}": bridge_topk_probability,
            "single_participation_out": single_components["participation_out"],
            "single_participation_in": single_components["participation_in"],
            "single_anchor_switch": single_components["anchor_switch"],
            "bootstrap_mean_participation_out": mean_components["participation_out"],
            "bootstrap_mean_participation_in": mean_components["participation_in"],
            "bootstrap_mean_anchor_switch": mean_components["anchor_switch"],
            "single_removal_impact": single_components["removal_impact"],
            "bootstrap_mean_removal_impact": mean_components["removal_impact"],
            "single_source_strength": row_sums_single,
            "bootstrap_mean_source_strength": row_sums_mean,
            "single_sink_strength": col_sums_single,
            "bootstrap_mean_sink_strength": col_sums_mean,
        }
    )
    if not community_assignments_df.empty:
        summary_df = summary_df.merge(
            community_assignments_df[
                [
                    "time_index",
                    "time_label",
                    "node_id",
                    "consensus_community",
                    "consensus_probability",
                    "community_instability",
                ]
            ],
            on=["time_index", "time_label", "node_id"],
            how="left",
        )
    summary_df = append_node_metadata(summary_df, reef_metadata)

    time_summary_df = pd.DataFrame(
        [
            {
                "time_index": time_index,
                "time_label": time_label,
                "bridge_top_k": top_k,
                "single_vs_mean_bridge_spearman": _safe_spearman(
                    single_components["bridge_score"],
                    mean_components["bridge_score"],
                ),
                "bridge_vs_source_spearman_single": _safe_spearman(
                    single_components["bridge_score"],
                    row_sums_single,
                ),
                "bridge_vs_sink_spearman_single": _safe_spearman(
                    single_components["bridge_score"],
                    col_sums_single,
                ),
                "bridge_vs_source_spearman_mean": _safe_spearman(
                    mean_components["bridge_score"],
                    row_sums_mean,
                ),
                "bridge_vs_sink_spearman_mean": _safe_spearman(
                    mean_components["bridge_score"],
                    col_sums_mean,
                ),
                "single_top_bridge_overlap_fraction": _top_k_overlap(
                    single_components["bridge_score"],
                    set(np.argsort(sample_mean)[-top_k:]),
                    top_k,
                ),
                "bootstrap_mean_top_bridge_overlap_fraction": _top_k_overlap(
                    mean_components["bridge_score"],
                    set(np.argsort(sample_mean)[-top_k:]),
                    top_k,
                ),
                "bootstrap_sample_top_bridge_overlap_median": float(
                    np.nanmedian(
                        [
                            _top_k_overlap(sample_bridge[sample_idx], set(np.argsort(sample_mean)[-top_k:]), top_k)
                            for sample_idx in range(sample_bridge.shape[0])
                        ]
                    )
                ),
                "median_bridge_cv": _safe_nanquantile(bridge_cv, 0.5),
                "p95_bridge_cv": _safe_nanquantile(bridge_cv, 0.95),
            }
        ]
    )

    diagnostics = {
        "single_bridge_score": np.asarray(single_components["bridge_score"], dtype=np.float64),
        "mean_bridge_score": np.asarray(mean_components["bridge_score"], dtype=np.float64),
        "sample_bridge_score": sample_bridge,
    }
    return summary_df, time_summary_df, diagnostics


def compute_distance_uncertainty(
    curve_edges: np.ndarray,
    curve_counts: np.ndarray,
    single_curve_sums: np.ndarray,
    mean_curve_sums: np.ndarray,
    sample_curve_sums: np.ndarray,
    single_curve_positive: np.ndarray,
    mean_curve_positive: np.ndarray,
    sample_curve_positive: np.ndarray,
    class_labels: list[str],
    single_class_sums: np.ndarray,
    mean_class_sums: np.ndarray,
    sample_class_sums: np.ndarray,
    time_index: int,
    time_label: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Summarize how connectivity and uncertainty vary by dispersal distance."""
    total_single = float(np.sum(single_curve_sums))
    total_mean = float(np.sum(mean_curve_sums))
    total_samples = np.sum(sample_curve_sums, axis=1)

    detail_rows: list[dict[str, Any]] = []
    for bin_index in range(curve_counts.size):
        sample_density = sample_curve_sums[:, bin_index] / curve_counts[bin_index] if curve_counts[bin_index] > 0 else np.full(sample_curve_sums.shape[0], np.nan)
        positive_sample_mean = np.divide(
            sample_curve_sums[:, bin_index],
            sample_curve_positive[:, bin_index],
            out=np.full(sample_curve_positive.shape[0], np.nan),
            where=sample_curve_positive[:, bin_index] > 0,
        )
        density_mean = _safe_nanmean(sample_density)
        detail_rows.append(
            {
                "time_index": time_index,
                "time_label": time_label,
                "distance_bin": bin_index,
                "bin_start_km": float(curve_edges[bin_index]),
                "bin_end_km": float(curve_edges[bin_index + 1]),
                "n_edges": int(curve_counts[bin_index]),
                "single_connectivity_density": float(single_curve_sums[bin_index] / curve_counts[bin_index]) if curve_counts[bin_index] > 0 else np.nan,
                "bootstrap_mean_connectivity_density": float(mean_curve_sums[bin_index] / curve_counts[bin_index]) if curve_counts[bin_index] > 0 else np.nan,
                "bootstrap_sample_density_mean": density_mean,
                "bootstrap_sample_density_q025": _safe_nanquantile(sample_density, 0.025),
                "bootstrap_sample_density_q975": _safe_nanquantile(sample_density, 0.975),
                "bootstrap_sample_density_cv": float(np.nanstd(sample_density, ddof=1) / density_mean)
                if np.isfinite(density_mean) and density_mean > 0
                else np.nan,
                "single_fraction_of_total": float(single_curve_sums[bin_index] / total_single) if total_single > 0 else np.nan,
                "bootstrap_mean_fraction_of_total": float(mean_curve_sums[bin_index] / total_mean) if total_mean > 0 else np.nan,
                "bootstrap_sample_fraction_mean": _safe_nanmean(sample_curve_sums[:, bin_index] / total_samples) if np.all(total_samples > 0) else np.nan,
                "single_positive_mean_weight": float(single_curve_sums[bin_index] / single_curve_positive[bin_index]) if single_curve_positive[bin_index] > 0 else np.nan,
                "bootstrap_mean_positive_mean_weight": float(mean_curve_sums[bin_index] / mean_curve_positive[bin_index]) if mean_curve_positive[bin_index] > 0 else np.nan,
                "bootstrap_sample_positive_mean_weight": _safe_nanmean(positive_sample_mean),
            }
        )
    detail_df = pd.DataFrame(detail_rows)

    class_rows: list[dict[str, Any]] = []
    for class_index, class_label in enumerate(class_labels):
        sample_fraction = np.divide(
            sample_class_sums[:, class_index],
            np.sum(sample_class_sums, axis=1),
            out=np.full(sample_class_sums.shape[0], np.nan),
            where=np.sum(sample_class_sums, axis=1) > 0,
        )
        class_rows.append(
            {
                "time_index": time_index,
                "time_label": time_label,
                "distance_class": class_label,
                "single_fraction_of_total": float(single_class_sums[class_index] / np.sum(single_class_sums)) if np.sum(single_class_sums) > 0 else np.nan,
                "bootstrap_mean_fraction_of_total": float(mean_class_sums[class_index] / np.sum(mean_class_sums)) if np.sum(mean_class_sums) > 0 else np.nan,
                "bootstrap_sample_fraction_mean": float(np.nanmean(sample_fraction)),
                "bootstrap_sample_fraction_q025": float(np.nanquantile(sample_fraction, 0.025)),
                "bootstrap_sample_fraction_q975": float(np.nanquantile(sample_fraction, 0.975)),
            }
        )
    class_df = pd.DataFrame(class_rows)

    diagnostics = {
        "sample_long_fraction": np.divide(
            sample_class_sums[:, -1],
            np.sum(sample_class_sums, axis=1),
            out=np.full(sample_class_sums.shape[0], np.nan),
            where=np.sum(sample_class_sums, axis=1) > 0,
        ),
        "single_long_fraction": float(single_class_sums[-1] / np.sum(single_class_sums)) if np.sum(single_class_sums) > 0 else np.nan,
        "mean_long_fraction": float(mean_class_sums[-1] / np.sum(mean_class_sums)) if np.sum(mean_class_sums) > 0 else np.nan,
    }
    return detail_df, class_df, diagnostics


def compute_regional_exchange(
    single_region_matrix: np.ndarray,
    mean_region_matrix: np.ndarray,
    sample_region_matrices: np.ndarray,
    region_labels: list[str],
    time_index: int,
    time_label: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Summarize region-to-region exchange and import/export roles."""
    pair_rows: list[dict[str, Any]] = []
    for source_region, source_label in enumerate(region_labels):
        for sink_region, sink_label in enumerate(region_labels):
            sample_values = sample_region_matrices[:, source_region, sink_region]
            q_low = float(np.nanquantile(sample_values, 0.025))
            q_high = float(np.nanquantile(sample_values, 0.975))
            pair_rows.append(
                {
                    "time_index": time_index,
                    "time_label": time_label,
                    "source_region": source_label,
                    "sink_region": sink_label,
                    "single_value": float(single_region_matrix[source_region, sink_region]),
                    "bootstrap_mean_value": float(mean_region_matrix[source_region, sink_region]),
                    "bootstrap_sample_mean": float(np.nanmean(sample_values)),
                    "bootstrap_q_low": q_low,
                    "bootstrap_q_high": q_high,
                    "single_outside_ci": bool((single_region_matrix[source_region, sink_region] < q_low) or (single_region_matrix[source_region, sink_region] > q_high)),
                    "bootstrap_mean_outside_ci": bool((mean_region_matrix[source_region, sink_region] < q_low) or (mean_region_matrix[source_region, sink_region] > q_high)),
                }
            )
    exchange_df = pd.DataFrame(pair_rows)

    role_rows: list[dict[str, Any]] = []
    for region_index, region_label in enumerate(region_labels):
        sample_exports = np.sum(sample_region_matrices[:, region_index, :], axis=1)
        sample_imports = np.sum(sample_region_matrices[:, :, region_index], axis=1)
        sample_balance = sample_exports - sample_imports
        role_rows.append(
            {
                "time_index": time_index,
                "time_label": time_label,
                "region": region_label,
                "single_export": float(np.sum(single_region_matrix[region_index, :])),
                "single_import": float(np.sum(single_region_matrix[:, region_index])),
                "single_balance": float(np.sum(single_region_matrix[region_index, :]) - np.sum(single_region_matrix[:, region_index])),
                "bootstrap_mean_export": float(np.sum(mean_region_matrix[region_index, :])),
                "bootstrap_mean_import": float(np.sum(mean_region_matrix[:, region_index])),
                "bootstrap_mean_balance": float(np.sum(mean_region_matrix[region_index, :]) - np.sum(mean_region_matrix[:, region_index])),
                "bootstrap_export_mean": float(np.nanmean(sample_exports)),
                "bootstrap_export_q025": float(np.nanquantile(sample_exports, 0.025)),
                "bootstrap_export_q975": float(np.nanquantile(sample_exports, 0.975)),
                "bootstrap_import_mean": float(np.nanmean(sample_imports)),
                "bootstrap_import_q025": float(np.nanquantile(sample_imports, 0.025)),
                "bootstrap_import_q975": float(np.nanquantile(sample_imports, 0.975)),
                "bootstrap_balance_mean": float(np.nanmean(sample_balance)),
                "bootstrap_balance_q025": float(np.nanquantile(sample_balance, 0.025)),
                "bootstrap_balance_q975": float(np.nanquantile(sample_balance, 0.975)),
            }
        )
    role_df = pd.DataFrame(role_rows)

    single_total = float(np.sum(single_region_matrix))
    mean_total = float(np.sum(mean_region_matrix))
    sample_totals = np.sum(sample_region_matrices, axis=(1, 2))
    single_inter = float(np.sum(single_region_matrix) - np.trace(single_region_matrix))
    mean_inter = float(np.sum(mean_region_matrix) - np.trace(mean_region_matrix))
    sample_inter = np.sum(sample_region_matrices, axis=(1, 2)) - np.trace(sample_region_matrices, axis1=1, axis2=2)

    summary_df = pd.DataFrame(
        [
            {
                "time_index": time_index,
                "time_label": time_label,
                "single_interregional_fraction": single_inter / single_total if single_total > 0 else np.nan,
                "bootstrap_mean_interregional_fraction": mean_inter / mean_total if mean_total > 0 else np.nan,
                "bootstrap_sample_interregional_fraction_median": float(np.nanmedian(sample_inter / sample_totals)) if np.all(sample_totals > 0) else np.nan,
                "single_balance_sd": float(np.nanstd(role_df["single_balance"], ddof=1)),
                "bootstrap_mean_balance_sd": float(np.nanstd(role_df["bootstrap_mean_balance"], ddof=1)),
                "bootstrap_sample_balance_sd_median": float(
                    np.nanmedian(np.nanstd(np.sum(sample_region_matrices, axis=2) - np.sum(sample_region_matrices, axis=1), axis=1, ddof=1))
                ),
            }
        ]
    )
    return exchange_df, role_df, summary_df


def compute_mean_bias(
    row_sums_single: np.ndarray,
    row_sums_mean: np.ndarray,
    row_sums_boot: np.ndarray,
    col_sums_single: np.ndarray,
    col_sums_mean: np.ndarray,
    col_sums_boot: np.ndarray,
    bridge_diagnostics: dict[str, Any],
    community_diagnostics: dict[str, Any],
    sample_top_indices: np.ndarray,
    single_matrix: np.ndarray,
    bootstrap_mean: np.ndarray,
    single_region_matrix: np.ndarray,
    mean_region_matrix: np.ndarray,
    sample_region_matrices: np.ndarray,
    distance_diagnostics: dict[str, Any],
    time_index: int,
    time_label: str,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Package interpretation-focused diagnostics for mean-matrix bias."""
    rows: list[dict[str, Any]] = []
    node_top_k = min(int(config["bias_top_k_nodes"]), row_sums_single.size)
    edge_top_k = min(int(config["bias_top_k_edges"]), single_matrix.size)
    region_top_k = min(int(config["bias_top_k_region_pairs"]), single_region_matrix.size)
    stable_rule = float(config["bias_stable_frequency"])

    source_consensus = set(np.argsort(np.nanmean(row_sums_boot, axis=0))[-node_top_k:])
    source_sample_overlap = np.array(
        [_top_k_overlap(row_sums_boot[sample_idx], source_consensus, node_top_k) for sample_idx in range(row_sums_boot.shape[0])],
        dtype=np.float64,
    )
    rows.extend(
        _distribution_rows(
            time_index,
            time_label,
            "source_ranking",
            f"top_{node_top_k}_overlap_to_consensus",
            source_sample_overlap,
            _top_k_overlap(row_sums_single, source_consensus, node_top_k),
            _top_k_overlap(row_sums_mean, source_consensus, node_top_k),
        )
    )

    sink_consensus = set(np.argsort(np.nanmean(col_sums_boot, axis=0))[-node_top_k:])
    sink_sample_overlap = np.array(
        [_top_k_overlap(col_sums_boot[sample_idx], sink_consensus, node_top_k) for sample_idx in range(col_sums_boot.shape[0])],
        dtype=np.float64,
    )
    rows.extend(
        _distribution_rows(
            time_index,
            time_label,
            "sink_ranking",
            f"top_{node_top_k}_overlap_to_consensus",
            sink_sample_overlap,
            _top_k_overlap(col_sums_single, sink_consensus, node_top_k),
            _top_k_overlap(col_sums_mean, sink_consensus, node_top_k),
        )
    )

    bridge_mean = np.asarray(bridge_diagnostics["mean_bridge_score"], dtype=np.float64)
    bridge_single = np.asarray(bridge_diagnostics["single_bridge_score"], dtype=np.float64)
    bridge_samples = np.asarray(bridge_diagnostics["sample_bridge_score"], dtype=np.float64)
    bridge_consensus = set(np.argsort(np.nanmean(bridge_samples, axis=0))[-node_top_k:])
    bridge_sample_overlap = np.array(
        [_top_k_overlap(bridge_samples[sample_idx], bridge_consensus, node_top_k) for sample_idx in range(bridge_samples.shape[0])],
        dtype=np.float64,
    )
    rows.extend(
        _distribution_rows(
            time_index,
            time_label,
            "bridge_reefs",
            f"top_{node_top_k}_overlap_to_consensus",
            bridge_sample_overlap,
            _top_k_overlap(bridge_single, bridge_consensus, node_top_k),
            _top_k_overlap(bridge_mean, bridge_consensus, node_top_k),
        )
    )

    consensus_labels = np.asarray(community_diagnostics["consensus_labels"], dtype=np.int32)
    mean_labels = np.asarray(community_diagnostics["mean_labels"], dtype=np.int32)
    single_labels = np.asarray(community_diagnostics["single_labels"], dtype=np.int32)
    sample_labels = np.asarray(community_diagnostics["sample_labels"], dtype=np.int32)
    community_sample_ari = np.array([adjusted_rand_index(consensus_labels, labels) for labels in sample_labels], dtype=np.float64)
    rows.extend(
        _distribution_rows(
            time_index,
            time_label,
            "community_structure",
            "ari_to_bootstrap_consensus",
            community_sample_ari,
            adjusted_rand_index(single_labels, consensus_labels),
            adjusted_rand_index(mean_labels, consensus_labels),
        )
    )

    counts = np.zeros(single_matrix.size, dtype=np.uint16)
    np.add.at(counts, sample_top_indices[:, :edge_top_k].ravel(), 1)
    stable_set = set(np.flatnonzero(counts >= int(np.ceil(stable_rule * sample_top_indices.shape[0]))))
    stable_sample_overlap = np.array(
        [_top_k_overlap_from_indices(sample_top_indices[sample_idx], stable_set, edge_top_k) for sample_idx in range(sample_top_indices.shape[0])],
        dtype=np.float64,
    )
    rows.extend(
        _distribution_rows(
            time_index,
            time_label,
            "top_links_backbone",
            f"top_{edge_top_k}_overlap_to_stable_edges",
            stable_sample_overlap,
            _top_k_overlap(single_matrix.ravel(), stable_set, edge_top_k),
            _top_k_overlap(bootstrap_mean.ravel(), stable_set, edge_top_k),
        )
    )

    offdiag_mask = np.ones(single_region_matrix.shape, dtype=bool)
    np.fill_diagonal(offdiag_mask, False)
    flat_mask = offdiag_mask.ravel()
    consensus_pairs = set(np.flatnonzero(flat_mask)[np.argsort(np.nanmean(sample_region_matrices.reshape(sample_region_matrices.shape[0], -1), axis=0)[flat_mask])[-region_top_k:]])
    sample_pair_overlap = np.array(
        [
            len(set(np.flatnonzero(flat_mask)[np.argsort(sample_region_matrices[sample_idx].ravel()[flat_mask])[-region_top_k:]]) & consensus_pairs)
            / float(region_top_k)
            for sample_idx in range(sample_region_matrices.shape[0])
        ],
        dtype=np.float64,
    )
    single_pair_overlap = len(set(np.flatnonzero(flat_mask)[np.argsort(single_region_matrix.ravel()[flat_mask])[-region_top_k:]]) & consensus_pairs) / float(region_top_k)
    mean_pair_overlap = len(set(np.flatnonzero(flat_mask)[np.argsort(mean_region_matrix.ravel()[flat_mask])[-region_top_k:]]) & consensus_pairs) / float(region_top_k)
    rows.extend(
        _distribution_rows(
            time_index,
            time_label,
            "regional_exchange",
            f"top_{region_top_k}_offdiag_pair_overlap_to_consensus",
            sample_pair_overlap,
            single_pair_overlap,
            mean_pair_overlap,
        )
    )

    rows.extend(
        _distribution_rows(
            time_index,
            time_label,
            "long_distance_connectivity",
            "long_distance_fraction_of_total",
            np.asarray(distance_diagnostics["sample_long_fraction"], dtype=np.float64),
            float(distance_diagnostics["single_long_fraction"]),
            float(distance_diagnostics["mean_long_fraction"]),
        )
    )
    return pd.DataFrame(rows)


def _planar_coordinates(latitudes: np.ndarray, longitudes: np.ndarray) -> np.ndarray:
    lat = np.asarray(latitudes, dtype=np.float64)
    lon = np.asarray(longitudes, dtype=np.float64)
    reference_lat = np.nanmean(lat) if np.isfinite(lat).any() else 0.0
    x = lon * np.cos(np.deg2rad(reference_lat))
    y = lat
    return np.column_stack([x, y])


def _neighbor_index(coords: np.ndarray, k_neighbors: int) -> np.ndarray:
    tree = cKDTree(coords)
    _, indices = tree.query(coords, k=min(k_neighbors + 1, coords.shape[0]))
    return indices[:, 1:]


def _morans_i(values: np.ndarray, neighbor_index: np.ndarray) -> float:
    finite = np.isfinite(values)
    if np.sum(finite) < 3:
        return np.nan
    centered = values - np.nanmean(values)
    denom = float(np.nansum(centered**2))
    if denom == 0.0:
        return np.nan
    numerator = 0.0
    n_weights = 0
    for node in np.flatnonzero(finite):
        neighbors = neighbor_index[node]
        neighbors = neighbors[finite[neighbors]]
        if neighbors.size == 0:
            continue
        numerator += float(centered[node] * np.sum(centered[neighbors]))
        n_weights += int(neighbors.size)
    if n_weights == 0:
        return np.nan
    return float((finite.sum() / n_weights) * (numerator / denom))


def compute_spatial_hotspots(
    source_rank_df: pd.DataFrame,
    sink_rank_df: pd.DataFrame,
    bridge_df: pd.DataFrame,
    community_assignments_df: pd.DataFrame,
    reef_metadata: pd.DataFrame | None,
    time_index: int,
    time_label: str,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Test whether disagreement metrics are spatially clustered."""
    if reef_metadata is None or reef_metadata.empty:
        return pd.DataFrame(), pd.DataFrame()
    if reef_metadata["LAT"].isna().all() or reef_metadata["LON"].isna().all():
        return pd.DataFrame(), pd.DataFrame()

    coords = _planar_coordinates(reef_metadata["LAT"].to_numpy(), reef_metadata["LON"].to_numpy())
    neighbors = _neighbor_index(coords, k_neighbors=int(config["spatial_k_neighbors"]))

    metric_frames = []
    if not source_rank_df.empty and f"prob_top_{int(config['rank_top_k_values'][0])}" in source_rank_df:
        metric_frames.append(
            (
                "source_topk_uncertainty",
                1.0 - source_rank_df[f"prob_top_{int(config['rank_top_k_values'][0])}"].to_numpy(dtype=np.float64),
            )
        )
    if not sink_rank_df.empty and f"prob_top_{int(config['rank_top_k_values'][0])}" in sink_rank_df:
        metric_frames.append(
            (
                "sink_topk_uncertainty",
                1.0 - sink_rank_df[f"prob_top_{int(config['rank_top_k_values'][0])}"].to_numpy(dtype=np.float64),
            )
        )
    if not bridge_df.empty and f"prob_top_bridge_{int(config['bridge_top_k'])}" in bridge_df:
        metric_frames.append(
            (
                "bridge_uncertainty",
                1.0 - bridge_df[f"prob_top_bridge_{int(config['bridge_top_k'])}"].to_numpy(dtype=np.float64),
            )
        )
    if not community_assignments_df.empty:
        metric_frames.append(
            (
                "community_instability",
                community_assignments_df["community_instability"].to_numpy(dtype=np.float64),
            )
        )

    summary_rows: list[dict[str, Any]] = []
    node_rows: list[dict[str, Any]] = []
    for metric_name, values in metric_frames:
        standardized = np.zeros_like(values, dtype=np.float64)
        if np.isfinite(values).any():
            mean = float(np.nanmean(values))
            std = float(np.nanstd(values))
            if std > 0:
                standardized = (values - mean) / std
        neighbor_mean = np.array(
            [float(np.nanmean(standardized[neighbors[node]])) for node in range(standardized.size)],
            dtype=np.float64,
        )
        hotspot_score = standardized * neighbor_mean
        summary_rows.append(
            {
                "time_index": time_index,
                "time_label": time_label,
                "metric_name": metric_name,
                "morans_i": _morans_i(values, neighbors),
                "median_hotspot_score": float(np.nanmedian(hotspot_score)),
                "p95_hotspot_score": float(np.nanquantile(hotspot_score, 0.95)),
            }
        )
        node_rows.extend(
            {
                "time_index": time_index,
                "time_label": time_label,
                "node_id": int(node_id),
                "metric_name": metric_name,
                "value": float(values[node_id]),
                "neighbor_mean_z": float(neighbor_mean[node_id]),
                "hotspot_score": float(hotspot_score[node_id]),
            }
            for node_id in range(values.size)
        )

    node_df = append_node_metadata(pd.DataFrame(node_rows), reef_metadata)
    return pd.DataFrame(summary_rows), node_df


def compute_variance_decomposition(metric_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Partition variation into bootstrap, spawning-period, year, and family components."""
    if metric_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    detail_df = metric_df.copy()
    summary_rows: list[dict[str, Any]] = []

    for metric_name, subset in metric_df.groupby("metric_name", sort=False):
        within_bootstrap = float(np.nanmean(subset["within_bootstrap_var"])) if "within_bootstrap_var" in subset else np.nan

        family_means = subset.groupby("family", dropna=False)["mean_value"].mean()
        among_family = float(np.nanvar(family_means.to_numpy(), ddof=0)) if family_means.size > 1 else np.nan

        year_subset = subset[subset["year"].notna()].copy() if "year" in subset else pd.DataFrame()
        year_group = year_subset.groupby(["family", "year"], dropna=False)["mean_value"].mean().reset_index()
        family_lookup = family_means.to_dict()
        year_residuals = np.array(
            [
                row.mean_value - family_lookup.get(row.family, np.nan)
                for row in year_group.itertuples(index=False)
            ],
            dtype=np.float64,
        )
        among_year = float(np.nanvar(year_residuals, ddof=0)) if np.isfinite(year_residuals).sum() > 1 else np.nan

        period_subset = subset[subset["year"].notna()].copy() if "year" in subset else pd.DataFrame()
        period_group = period_subset.groupby(["family", "year", "spawning_period"], dropna=False)["mean_value"].mean().reset_index()
        year_lookup = {
            (row.family, row.year): row.mean_value
            for row in year_group.itertuples(index=False)
        }
        period_residuals = np.array(
            [
                row.mean_value - year_lookup.get((row.family, row.year), np.nan)
                for row in period_group.itertuples(index=False)
            ],
            dtype=np.float64,
        )
        among_period = float(np.nanvar(period_residuals, ddof=0)) if np.isfinite(period_residuals).sum() > 1 else np.nan

        components = {
            "within_time_bootstrap": within_bootstrap,
            "among_spawning_periods": among_period,
            "among_years": among_year,
            "among_families": among_family,
        }
        finite_components = {key: value for key, value in components.items() if np.isfinite(value)}
        total = float(np.nansum(list(finite_components.values()))) if finite_components else np.nan
        dominant = max(finite_components, key=finite_components.get) if finite_components else "insufficient_data"

        summary_rows.append(
            {
                "metric_name": metric_name,
                "within_time_bootstrap_var": within_bootstrap,
                "among_spawning_periods_var": among_period,
                "among_years_var": among_year,
                "among_families_var": among_family,
                "approx_total_var": total,
                "within_vs_temporal_ratio": within_bootstrap / (among_period + among_year)
                if np.isfinite(within_bootstrap) and np.isfinite(among_period + among_year) and (among_period + among_year) > 0
                else np.nan,
                "dominant_source": dominant,
                "n_families": int(subset["family"].nunique(dropna=True)),
                "n_time_points": int(subset["time_label"].nunique()),
            }
        )
    return pd.DataFrame(summary_rows), detail_df


def compute_ecological_thresholds(
    bootstrap_mean: np.ndarray,
    threshold_presence_counts: dict[float, np.ndarray],
    anchor_onehot: np.ndarray,
    time_index: int,
    time_label: str,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Evaluate higher-level structural conclusions under ecological threshold modes."""
    mean_positive = np.where(np.isfinite(bootstrap_mean), np.clip(bootstrap_mean, 0.0, None), 0.0)
    row_sums_mean = np.sum(mean_positive, axis=1)

    baseline_out = mean_positive @ anchor_onehot
    baseline_in = mean_positive.T @ anchor_onehot
    baseline_features = build_directed_profile_features(baseline_out, baseline_in, row_sums_mean, np.sum(mean_positive, axis=0))
    baseline_k = choose_cluster_count(
        baseline_features,
        min_clusters=int(config["community_min_clusters"]),
        max_clusters=int(config["community_max_clusters"]),
        seed=int(config["random_seed"]) + time_index + 500,
    )
    baseline_labels = _cluster_labels(baseline_features, baseline_k, seed=int(config["random_seed"]) + time_index + 501)
    baseline_bridge = _bridge_components(baseline_out, baseline_in)["bridge_score"]
    baseline_bridge_top = set(np.argsort(baseline_bridge)[-int(config["bridge_top_k"]):])
    baseline_mask = mean_positive > 0

    flat_values = mean_positive.ravel()
    positive_values = flat_values[flat_values > 0]
    percentile_cut = float(np.nanquantile(positive_values, 0.995)) if positive_values.size else np.inf
    sorted_positive = np.sort(positive_values)[::-1]
    cumulative_target = 0.9 * float(np.sum(sorted_positive)) if sorted_positive.size else np.nan
    cumulative_cut = np.inf
    if sorted_positive.size and np.isfinite(cumulative_target):
        cumulative_index = int(np.searchsorted(np.cumsum(sorted_positive), cumulative_target, side="left"))
        cumulative_cut = float(sorted_positive[min(cumulative_index, sorted_positive.size - 1)])

    persistence_counts = threshold_presence_counts.get(0.0)
    stable_counts = threshold_presence_counts.get(1e-8, persistence_counts)
    rows: list[dict[str, Any]] = []
    mode_masks = [
        ("absolute", mean_positive > 1e-8),
        ("percentile", mean_positive >= percentile_cut),
        ("cumulative_flow", mean_positive >= cumulative_cut),
        (
            "persistence",
            (persistence_counts / max(1, int(config.get("n_bootstrap_samples", 1)))) >= 0.5 if persistence_counts is not None else np.zeros_like(mean_positive, dtype=bool),
        ),
        (
            "row_share",
            np.divide(mean_positive, row_sums_mean[:, None], out=np.zeros_like(mean_positive), where=row_sums_mean[:, None] > 0) >= float(config["local_share_thresholds"][0]),
        ),
        (
            "stable_edge",
            (stable_counts / max(1, int(config.get("n_bootstrap_samples", 1)))) >= float(config["bias_stable_frequency"]) if stable_counts is not None else np.zeros_like(mean_positive, dtype=bool),
        ),
    ]

    for mode_name, mode_mask in mode_masks:
        masked = np.where(mode_mask, mean_positive, 0.0)
        out_anchor = masked @ anchor_onehot
        in_anchor = masked.T @ anchor_onehot
        row_sums = np.sum(masked, axis=1)
        col_sums = np.sum(masked, axis=0)
        features = build_directed_profile_features(out_anchor, in_anchor, row_sums, col_sums)
        selected_k = choose_cluster_count(
            features,
            min_clusters=int(config["community_min_clusters"]),
            max_clusters=int(config["community_max_clusters"]),
            seed=int(config["random_seed"]) + time_index + 600 + len(rows),
        )
        labels = _cluster_labels(features, selected_k, seed=int(config["random_seed"]) + time_index + 700 + len(rows))
        bridge_score = _bridge_components(out_anchor, in_anchor)["bridge_score"]
        bridge_top = set(np.argsort(bridge_score)[-int(config["bridge_top_k"]):])
        overlap = int(np.sum(mode_mask & baseline_mask))
        union = int(np.sum(mode_mask | baseline_mask))
        rows.append(
            {
                "time_index": time_index,
                "time_label": time_label,
                "threshold_mode": mode_name,
                "active_edge_count": int(np.sum(mode_mask)),
                "active_edge_fraction": float(np.mean(mode_mask)),
                "community_count": int(np.unique(labels).size),
                "community_ari_vs_baseline": adjusted_rand_index(labels, baseline_labels),
                "bridge_top_overlap_fraction": len(bridge_top & baseline_bridge_top) / float(max(1, int(config["bridge_top_k"]))),
                "backbone_jaccard_vs_baseline": overlap / union if union else np.nan,
                "inter_anchor_fraction": float((np.sum(masked) - np.trace(anchor_onehot.T @ masked @ anchor_onehot)) / np.sum(masked))
                if np.sum(masked) > 0
                else np.nan,
            }
        )
    return pd.DataFrame(rows)
