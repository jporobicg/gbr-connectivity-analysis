"""Extended structure, strength, backbone, and uncertainty analyses."""

from __future__ import annotations

from typing import Iterable
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.vq import kmeans2


def _safe_pearson(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2:
        return np.nan
    if float(np.std(x)) == 0.0 or float(np.std(y)) == 0.0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def rank_desc(values: np.ndarray) -> np.ndarray:
    """Return descending ranks using average ties."""
    return stats.rankdata(-np.asarray(values, dtype=np.float64), method="average")


def topk_membership_from_samples(values: np.ndarray, top_k_values: Iterable[int]) -> dict[int, np.ndarray]:
    """Return top-k membership frequencies from bootstrap sample vectors."""
    values = np.asarray(values, dtype=np.float64)
    n_samples, n_nodes = values.shape
    membership: dict[int, np.ndarray] = {}
    order = np.argsort(values, axis=1)[:, ::-1]
    for top_k in sorted(set(int(value) for value in top_k_values)):
        count = np.zeros(n_nodes, dtype=np.float64)
        top_idx = order[:, : min(top_k, n_nodes)]
        for sample_idx in range(n_samples):
            count[top_idx[sample_idx]] += 1
        membership[top_k] = count / n_samples
    return membership


def compute_rank_uncertainty(
    single_values: np.ndarray,
    bootstrap_mean_values: np.ndarray,
    bootstrap_samples: np.ndarray,
    entity_type: str,
    top_k_values: Iterable[int],
    time_index: int,
    time_label: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarize rank stability against the full bootstrap distribution."""
    single_rank = rank_desc(single_values)
    bootstrap_mean_rank = rank_desc(bootstrap_mean_values)

    sample_ranks = np.vstack([rank_desc(sample) for sample in bootstrap_samples])
    prob_topk = topk_membership_from_samples(bootstrap_samples, top_k_values)

    q_low = np.nanquantile(sample_ranks, 0.025, axis=0)
    q_high = np.nanquantile(sample_ranks, 0.975, axis=0)
    outside = (single_rank < q_low) | (single_rank > q_high)

    detail = pd.DataFrame(
        {
            "time_index": time_index,
            "time_label": time_label,
            "entity_type": entity_type,
            "node_id": np.arange(single_values.size, dtype=int),
            "single_strength": single_values,
            "bootstrap_mean_strength": bootstrap_mean_values,
            "single_rank": single_rank,
            "bootstrap_mean_rank": bootstrap_mean_rank,
            "bootstrap_rank_mean": np.nanmean(sample_ranks, axis=0),
            "bootstrap_rank_median": np.nanmedian(sample_ranks, axis=0),
            "bootstrap_rank_sd": np.nanstd(sample_ranks, axis=0, ddof=1),
            "bootstrap_rank_q025": q_low,
            "bootstrap_rank_q975": q_high,
            "single_rank_outside_ci": outside,
            "single_vs_bootstrap_mean_rank_shift": bootstrap_mean_rank - single_rank,
        }
    )
    for top_k, probability in prob_topk.items():
        detail[f"prob_top_{top_k}"] = probability

    summary_rows: list[dict[str, float | int | str]] = []
    for top_k, probability in prob_topk.items():
        sample_spearman = np.array(
            [stats.spearmanr(single_values, sample).statistic for sample in bootstrap_samples],
            dtype=np.float64,
        )
        mean_sample_spearman = np.array(
            [stats.spearmanr(bootstrap_mean_values, sample).statistic for sample in bootstrap_samples],
            dtype=np.float64,
        )
        summary_rows.append(
            {
                "time_index": time_index,
                "time_label": time_label,
                "entity_type": entity_type,
                "top_k": int(top_k),
                "single_vs_mean_spearman": float(stats.spearmanr(single_values, bootstrap_mean_values).statistic),
                "single_vs_samples_spearman_median": float(np.nanmedian(sample_spearman)),
                "mean_vs_samples_spearman_median": float(np.nanmedian(mean_sample_spearman)),
                "single_rank_outside_ci_fraction": float(np.mean(outside)),
                "mean_prob_top_k": float(np.mean(probability)),
                "median_prob_top_k": float(np.median(probability)),
                "max_prob_top_k": float(np.max(probability)),
            }
        )

    return pd.DataFrame(summary_rows), detail.sort_values(
        ["single_rank_outside_ci", "bootstrap_rank_sd"],
        ascending=[False, False],
    ).reset_index(drop=True)


def _js_divergence(p: np.ndarray, q: np.ndarray) -> float:
    m = 0.5 * (p + q)
    with np.errstate(divide="ignore", invalid="ignore"):
        kl_pm = np.where((p > 0) & (m > 0), p * np.log(p / m), 0.0)
        kl_qm = np.where((q > 0) & (m > 0), q * np.log(q / m), 0.0)
    return float(0.5 * np.sum(kl_pm) + 0.5 * np.sum(kl_qm))


def _hellinger_distance(p: np.ndarray, q: np.ndarray) -> float:
    return float(np.sqrt(0.5 * np.sum((np.sqrt(p) - np.sqrt(q)) ** 2)))


def compute_profile_divergence(
    single_matrix: np.ndarray,
    bootstrap_mean: np.ndarray,
    row_sums_single: np.ndarray,
    row_sums_mean: np.ndarray,
    col_sums_single: np.ndarray,
    col_sums_mean: np.ndarray,
    time_index: int,
    time_label: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare normalized destination and origin profiles."""
    rows: list[dict[str, float | int | str]] = []
    for node_id in range(single_matrix.shape[0]):
        if row_sums_single[node_id] > 0 and row_sums_mean[node_id] > 0:
            p = np.clip(single_matrix[node_id] / row_sums_single[node_id], 0.0, None)
            q = np.clip(bootstrap_mean[node_id] / row_sums_mean[node_id], 0.0, None)
            rows.append(
                {
                    "time_index": time_index,
                    "time_label": time_label,
                    "entity_type": "source",
                    "node_id": node_id,
                    "js_divergence": _js_divergence(p, q),
                    "hellinger_distance": _hellinger_distance(p, q),
                    "single_total_strength": float(row_sums_single[node_id]),
                    "bootstrap_mean_total_strength": float(row_sums_mean[node_id]),
                }
            )

    for node_id in range(single_matrix.shape[1]):
        if col_sums_single[node_id] > 0 and col_sums_mean[node_id] > 0:
            p = np.clip(single_matrix[:, node_id] / col_sums_single[node_id], 0.0, None)
            q = np.clip(bootstrap_mean[:, node_id] / col_sums_mean[node_id], 0.0, None)
            rows.append(
                {
                    "time_index": time_index,
                    "time_label": time_label,
                    "entity_type": "sink",
                    "node_id": node_id,
                    "js_divergence": _js_divergence(p, q),
                    "hellinger_distance": _hellinger_distance(p, q),
                    "single_total_strength": float(col_sums_single[node_id]),
                    "bootstrap_mean_total_strength": float(col_sums_mean[node_id]),
                }
            )

    detail = pd.DataFrame(rows)
    summary = (
        detail.groupby(["time_index", "time_label", "entity_type"], as_index=False)
        .agg(
            mean_js_divergence=("js_divergence", "mean"),
            median_js_divergence=("js_divergence", "median"),
            p95_js_divergence=("js_divergence", lambda values: float(np.quantile(values, 0.95))),
            mean_hellinger_distance=("hellinger_distance", "mean"),
            median_hellinger_distance=("hellinger_distance", "median"),
        )
    )
    return summary, detail.sort_values("js_divergence", ascending=False).reset_index(drop=True)


def _edge_distribution_stats(
    values: np.ndarray,
    top_share_fracs: Iterable[float],
    top_k_values: Iterable[int],
) -> dict[str, float]:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {}

    sorted_values = np.sort(values)
    total = float(np.sum(sorted_values))
    positive = sorted_values[sorted_values > 0]
    stats_row = {
        "gini": _gini(sorted_values),
        "max_edge": float(sorted_values[-1]),
        "positive_edge_fraction": float(np.mean(sorted_values > 0)),
        "overall_mean": float(np.mean(sorted_values)),
        "overall_median": float(np.median(sorted_values)),
        "positive_median": float(np.median(positive)) if positive.size else np.nan,
        "positive_p90": float(np.quantile(positive, 0.90)) if positive.size else np.nan,
        "positive_p95": float(np.quantile(positive, 0.95)) if positive.size else np.nan,
        "positive_p99": float(np.quantile(positive, 0.99)) if positive.size else np.nan,
    }

    for frac in top_share_fracs:
        top_n = max(1, int(np.ceil(sorted_values.size * frac)))
        share = float(np.sum(sorted_values[-top_n:]) / total) if total > 0 else np.nan
        stats_row[f"share_top_{int(frac * 100)}pct"] = share

    for top_k in top_k_values:
        k = min(int(top_k), sorted_values.size)
        top_values = sorted_values[-k:]
        stats_row[f"top_{k}_mean"] = float(np.mean(top_values))
        stats_row[f"top_{k}_sum"] = float(np.sum(top_values))

    if positive.size:
        stats_row["strongest_to_positive_median"] = float(sorted_values[-1] / np.median(positive))
        top_1pct_n = max(1, int(np.ceil(sorted_values.size * 0.01)))
        stats_row["top_1pct_mean_to_overall_mean"] = (
            float(np.mean(sorted_values[-top_1pct_n:]) / stats_row["overall_mean"]) if stats_row["overall_mean"] > 0 else np.nan
        )
    else:
        stats_row["strongest_to_positive_median"] = np.nan
        stats_row["top_1pct_mean_to_overall_mean"] = np.nan
    return stats_row


def _gini(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    if np.all(values == 0):
        return 0.0
    values = np.sort(values)
    index = np.arange(1, values.size + 1)
    return float((2 * np.sum(index * values) / (values.size * np.sum(values))) - (values.size + 1) / values.size)


def lorenz_curve_points(
    values: np.ndarray,
    dataset: str,
    time_index: int,
    time_label: str,
    n_points: int = 100,
) -> pd.DataFrame:
    """Compute Lorenz curve coordinates."""
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return pd.DataFrame()

    values = np.sort(np.clip(values, 0.0, None))
    cum_values = np.cumsum(values)
    total = cum_values[-1] if cum_values.size else 0.0
    edge_share = np.linspace(0.0, 1.0, values.size + 1)
    if total > 0:
        weight_share = np.concatenate([[0.0], cum_values / total])
    else:
        weight_share = edge_share.copy()

    sample_points = np.linspace(0, values.size, num=min(n_points, values.size + 1), dtype=int)
    return pd.DataFrame(
        {
            "time_index": time_index,
            "time_label": time_label,
            "dataset": dataset,
            "edge_fraction": edge_share[sample_points],
            "weight_fraction": weight_share[sample_points],
        }
    )


def compute_concentration_and_tail(
    single_matrix: np.ndarray,
    bootstrap_mean: np.ndarray,
    sampled_bootstrap: np.ndarray,
    sample_top_values: np.ndarray,
    time_index: int,
    time_label: str,
    top_share_fracs: Iterable[float],
    top_k_values: Iterable[int],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Summarize concentration and upper-tail behavior."""
    rows: list[dict[str, float | int | str]] = []
    lorenz_rows = [
        lorenz_curve_points(single_matrix.ravel(), "single", time_index, time_label),
        lorenz_curve_points(bootstrap_mean.ravel(), "bootstrap_mean", time_index, time_label),
    ]

    for dataset, values in (("single", single_matrix.ravel()), ("bootstrap_mean", bootstrap_mean.ravel())):
        row = {
            "time_index": time_index,
            "time_label": time_label,
            "dataset": dataset,
        }
        row.update(_edge_distribution_stats(values, top_share_fracs, top_k_values))
        rows.append(row)

    sample_rows: list[dict[str, float | int | str]] = []
    positive_bootstrap = np.where(sampled_bootstrap > 0, sampled_bootstrap, np.nan)
    for sample_idx in range(sampled_bootstrap.shape[0]):
        row = {
            "time_index": time_index,
            "time_label": time_label,
            "dataset": "bootstrap_sample",
            "sample_id": sample_idx,
            "sampled_edge_gini": _gini(np.nan_to_num(sampled_bootstrap[sample_idx], nan=0.0)),
            "sampled_positive_median": float(np.nanmedian(positive_bootstrap[sample_idx])),
            "sampled_positive_p95": float(np.nanquantile(positive_bootstrap[sample_idx], 0.95)),
            "max_edge": float(np.max(sample_top_values[sample_idx])),
        }
        for top_k in top_k_values:
            k = min(int(top_k), sample_top_values.shape[1])
            row[f"top_{k}_mean"] = float(np.mean(np.sort(sample_top_values[sample_idx])[-k:]))
        sample_rows.append(row)

    return pd.DataFrame(rows), pd.concat(lorenz_rows, ignore_index=True), pd.DataFrame(sample_rows)


def compute_node_strength_uncertainty(
    single_values: np.ndarray,
    bootstrap_mean_values: np.ndarray,
    bootstrap_samples: np.ndarray,
    entity_type: str,
    time_index: int,
    time_label: str,
) -> pd.DataFrame:
    """Return full node-level uncertainty summaries."""
    bootstrap_sample_mean = np.nanmean(bootstrap_samples, axis=0)
    bootstrap_sample_std = np.nanstd(bootstrap_samples, axis=0, ddof=1)
    q_low = np.nanquantile(bootstrap_samples, 0.025, axis=0)
    q_high = np.nanquantile(bootstrap_samples, 0.975, axis=0)
    cv = np.full(bootstrap_sample_std.shape, np.nan, dtype=np.float64)
    valid = np.abs(bootstrap_sample_mean) > 1e-12
    cv[valid] = bootstrap_sample_std[valid] / np.abs(bootstrap_sample_mean[valid])
    z_score = np.full(bootstrap_sample_std.shape, np.nan, dtype=np.float64)
    valid_std = bootstrap_sample_std > 0
    z_score[valid_std] = (single_values[valid_std] - bootstrap_sample_mean[valid_std]) / bootstrap_sample_std[valid_std]

    return pd.DataFrame(
        {
            "time_index": time_index,
            "time_label": time_label,
            "entity_type": entity_type,
            "node_id": np.arange(single_values.size, dtype=int),
            "single_value": single_values,
            "bootstrap_mean_value": bootstrap_mean_values,
            "bootstrap_sample_mean": bootstrap_sample_mean,
            "bootstrap_sample_sd": bootstrap_sample_std,
            "bootstrap_sample_cv": cv,
            "bootstrap_q_low": q_low,
            "bootstrap_q_high": q_high,
            "single_outside_ci": (single_values < q_low) | (single_values > q_high),
            "single_minus_bootstrap_mean": single_values - bootstrap_sample_mean,
            "single_z_score": z_score,
        }
    )


def compute_local_node_metrics(
    matrix: np.ndarray,
    row_sums: np.ndarray,
    col_sums: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute entropy, effective degree, and dominance for rows and columns."""
    positive = np.where(matrix > 0, matrix, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        w_log_w = np.where(positive > 0, positive * np.log(positive), 0.0)
    row_w_log_w = np.sum(w_log_w, axis=1)
    col_w_log_w = np.sum(w_log_w, axis=0)
    row_max = np.max(positive, axis=1)
    col_max = np.max(positive, axis=0)

    row_entropy = np.full(row_sums.shape, np.nan, dtype=np.float64)
    col_entropy = np.full(col_sums.shape, np.nan, dtype=np.float64)
    valid_row = row_sums > 0
    valid_col = col_sums > 0
    row_entropy[valid_row] = np.log(row_sums[valid_row]) - (row_w_log_w[valid_row] / row_sums[valid_row])
    col_entropy[valid_col] = np.log(col_sums[valid_col]) - (col_w_log_w[valid_col] / col_sums[valid_col])

    row_effective = np.exp(row_entropy, where=np.isfinite(row_entropy), out=np.full(row_entropy.shape, np.nan))
    col_effective = np.exp(col_entropy, where=np.isfinite(col_entropy), out=np.full(col_entropy.shape, np.nan))
    row_dominance = np.full(row_sums.shape, np.nan, dtype=np.float64)
    col_dominance = np.full(col_sums.shape, np.nan, dtype=np.float64)
    row_dominance[valid_row] = row_max[valid_row] / row_sums[valid_row]
    col_dominance[valid_col] = col_max[valid_col] / col_sums[valid_col]
    return row_entropy, row_effective, row_dominance, col_entropy, col_effective, col_dominance


def summarize_local_node_metrics(
    single_row_metrics: tuple[np.ndarray, np.ndarray, np.ndarray],
    mean_row_metrics: tuple[np.ndarray, np.ndarray, np.ndarray],
    bootstrap_row_metrics: tuple[np.ndarray, np.ndarray, np.ndarray],
    single_col_metrics: tuple[np.ndarray, np.ndarray, np.ndarray],
    mean_col_metrics: tuple[np.ndarray, np.ndarray, np.ndarray],
    bootstrap_col_metrics: tuple[np.ndarray, np.ndarray, np.ndarray],
    time_index: int,
    time_label: str,
) -> pd.DataFrame:
    """Combine local-structure summaries for sources and sinks."""
    tables = []
    for entity_type, single_metrics, mean_metrics, boot_metrics in (
        ("source", single_row_metrics, mean_row_metrics, bootstrap_row_metrics),
        ("sink", single_col_metrics, mean_col_metrics, bootstrap_col_metrics),
    ):
        entropy_single, effective_single, dominance_single = single_metrics
        entropy_mean, effective_mean, dominance_mean = mean_metrics
        entropy_boot, effective_boot, dominance_boot = boot_metrics
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            entropy_boot_mean = np.nanmean(entropy_boot, axis=0)
            entropy_boot_sd = np.nanstd(entropy_boot, axis=0, ddof=1)
            effective_boot_mean = np.nanmean(effective_boot, axis=0)
            effective_boot_sd = np.nanstd(effective_boot, axis=0, ddof=1)
            dominance_boot_mean = np.nanmean(dominance_boot, axis=0)
            dominance_boot_sd = np.nanstd(dominance_boot, axis=0, ddof=1)
        table = pd.DataFrame(
            {
                "time_index": time_index,
                "time_label": time_label,
                "entity_type": entity_type,
                "node_id": np.arange(entropy_single.size, dtype=int),
                "single_entropy": entropy_single,
                "bootstrap_mean_entropy": entropy_mean,
                "bootstrap_sample_entropy_mean": entropy_boot_mean,
                "bootstrap_sample_entropy_sd": entropy_boot_sd,
                "single_effective_n": effective_single,
                "bootstrap_mean_effective_n": effective_mean,
                "bootstrap_sample_effective_n_mean": effective_boot_mean,
                "bootstrap_sample_effective_n_sd": effective_boot_sd,
                "single_dominance": dominance_single,
                "bootstrap_mean_dominance": dominance_mean,
                "bootstrap_sample_dominance_mean": dominance_boot_mean,
                "bootstrap_sample_dominance_sd": dominance_boot_sd,
            }
        )
        tables.append(table)
    return pd.concat(tables, ignore_index=True)


def _standardize_features(features: np.ndarray) -> np.ndarray:
    mean = np.nanmean(features, axis=0, keepdims=True)
    std = np.nanstd(features, axis=0, keepdims=True)
    std[std == 0] = 1.0
    return np.nan_to_num((features - mean) / std, nan=0.0)


def adjusted_rand_index(labels_a: np.ndarray, labels_b: np.ndarray) -> float:
    """Compute adjusted Rand index without an extra dependency."""
    labels_a = np.asarray(labels_a, dtype=int)
    labels_b = np.asarray(labels_b, dtype=int)
    n = labels_a.size
    if n == 0:
        return np.nan

    _, inverse_a = np.unique(labels_a, return_inverse=True)
    _, inverse_b = np.unique(labels_b, return_inverse=True)
    contingency = np.zeros((inverse_a.max() + 1, inverse_b.max() + 1), dtype=np.int64)
    np.add.at(contingency, (inverse_a, inverse_b), 1)

    def comb2(values: np.ndarray) -> np.ndarray:
        return values * (values - 1) / 2

    sum_comb = np.sum(comb2(contingency))
    sum_a = np.sum(comb2(contingency.sum(axis=1)))
    sum_b = np.sum(comb2(contingency.sum(axis=0)))
    total = comb2(np.array([n], dtype=np.int64))[0]
    expected = (sum_a * sum_b) / total if total else 0.0
    max_index = 0.5 * (sum_a + sum_b)
    denom = max_index - expected
    if denom == 0:
        return np.nan
    return float((sum_comb - expected) / denom)


def cluster_node_features(features: np.ndarray, n_groups: int, seed: int) -> np.ndarray:
    """Cluster node features with a simple k-means approach."""
    standardized = _standardize_features(features)
    _, labels = kmeans2(standardized, n_groups, minit="points", seed=seed)
    return labels.astype(int)


def compute_grouping_summary(
    single_features: np.ndarray,
    mean_features: np.ndarray,
    sample_features: np.ndarray,
    n_groups: int,
    time_index: int,
    time_label: str,
    random_seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare broad node groupings between products and across bootstrap samples."""
    single_labels = cluster_node_features(single_features, n_groups=n_groups, seed=random_seed + time_index)
    mean_labels = cluster_node_features(mean_features, n_groups=n_groups, seed=random_seed + time_index + 1)
    sample_labels = [
        cluster_node_features(sample_features[sample_idx], n_groups=n_groups, seed=random_seed + time_index + 10 + sample_idx)
        for sample_idx in range(sample_features.shape[0])
    ]

    mean_vs_samples = np.array([adjusted_rand_index(mean_labels, labels) for labels in sample_labels], dtype=np.float64)
    single_vs_samples = np.array([adjusted_rand_index(single_labels, labels) for labels in sample_labels], dtype=np.float64)

    summary = pd.DataFrame(
        [
            {
                "time_index": time_index,
                "time_label": time_label,
                "n_groups": n_groups,
                "single_vs_bootstrap_mean_ari": adjusted_rand_index(single_labels, mean_labels),
                "bootstrap_mean_vs_samples_ari_median": float(np.nanmedian(mean_vs_samples)),
                "bootstrap_mean_vs_samples_ari_q025": float(np.nanquantile(mean_vs_samples, 0.025)),
                "bootstrap_mean_vs_samples_ari_q975": float(np.nanquantile(mean_vs_samples, 0.975)),
                "single_vs_samples_ari_median": float(np.nanmedian(single_vs_samples)),
                "single_vs_samples_ari_q025": float(np.nanquantile(single_vs_samples, 0.025)),
                "single_vs_samples_ari_q975": float(np.nanquantile(single_vs_samples, 0.975)),
            }
        ]
    )

    rows: list[dict[str, float | int | str]] = []
    for dataset, labels in (("single", single_labels), ("bootstrap_mean", mean_labels)):
        counts = np.bincount(labels, minlength=n_groups)
        for group_id, size in enumerate(np.sort(counts)[::-1], start=1):
            rows.append(
                {
                    "time_index": time_index,
                    "time_label": time_label,
                    "dataset": dataset,
                    "group_order": group_id,
                    "group_size": int(size),
                }
            )

    sample_size_matrix = np.vstack([np.sort(np.bincount(labels, minlength=n_groups))[::-1] for labels in sample_labels])
    for group_id in range(n_groups):
        rows.append(
            {
                "time_index": time_index,
                "time_label": time_label,
                "dataset": "bootstrap_samples_median",
                "group_order": group_id + 1,
                "group_size": float(np.median(sample_size_matrix[:, group_id])),
            }
        )
    return summary, pd.DataFrame(rows)


def compute_edge_churn(
    single_matrix: np.ndarray,
    bootstrap_mean: np.ndarray,
    top_k_values: Iterable[int],
    time_index: int,
    time_label: str,
    top_n_per_status: int = 50,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute top-k edge overlap plus entering/leaving links."""
    valid = np.isfinite(single_matrix) & np.isfinite(bootstrap_mean)
    single_flat = np.where(valid, single_matrix, -np.inf).ravel()
    mean_flat = np.where(valid, bootstrap_mean, -np.inf).ravel()
    rows: list[dict[str, float | int | str]] = []
    detail_rows: list[dict[str, float | int | str]] = []

    for top_k in sorted(set(int(value) for value in top_k_values)):
        single_idx = np.argpartition(single_flat, -top_k)[-top_k:]
        mean_idx = np.argpartition(mean_flat, -top_k)[-top_k:]
        single_set = set(map(int, single_idx))
        mean_set = set(map(int, mean_idx))
        overlap = single_set & mean_set
        entering = list(mean_set - single_set)
        leaving = list(single_set - mean_set)
        rows.append(
            {
                "time_index": time_index,
                "time_label": time_label,
                "top_k": top_k,
                "overlap_count": len(overlap),
                "overlap_fraction": len(overlap) / top_k,
                "jaccard_index": len(overlap) / len(single_set | mean_set) if (single_set | mean_set) else np.nan,
                "entering_count": len(entering),
                "leaving_count": len(leaving),
                "churn_fraction": (len(entering) + len(leaving)) / (2 * top_k),
            }
        )
        for status, index_values, weight_array in (
            ("entering", entering, mean_flat),
            ("leaving", leaving, single_flat),
        ):
            if not index_values:
                continue
            ordered = sorted(index_values, key=lambda idx: weight_array[idx], reverse=True)[:top_n_per_status]
            source_idx, sink_idx = np.unravel_index(np.array(ordered, dtype=np.int64), single_matrix.shape)
            for source, sink, flat_idx in zip(source_idx, sink_idx, ordered, strict=True):
                detail_rows.append(
                    {
                        "time_index": time_index,
                        "time_label": time_label,
                        "top_k": top_k,
                        "status": status,
                        "source": int(source),
                        "sink": int(sink),
                        "single_value": float(single_matrix[source, sink]),
                        "bootstrap_mean_value": float(bootstrap_mean[source, sink]),
                    }
                )

    return pd.DataFrame(rows), pd.DataFrame(detail_rows)


def compute_stable_edge_summary(
    single_matrix: np.ndarray,
    bootstrap_mean: np.ndarray,
    sample_top_indices: np.ndarray,
    threshold_presence_counts: dict[float, np.ndarray],
    top_k_values: Iterable[int],
    stable_frequency_rules: Iterable[float],
    time_index: int,
    time_label: str,
    top_n_details: int = 100,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarize stable edges across bootstrap samples."""
    total_edges = single_matrix.size
    n_samples = sample_top_indices.shape[0]
    rows: list[dict[str, float | int | str]] = []
    detail_rows: list[dict[str, float | int | str]] = []

    single_flat = single_matrix.ravel()
    mean_flat = bootstrap_mean.ravel()
    for top_k in sorted(set(int(value) for value in top_k_values)):
        counts = np.zeros(total_edges, dtype=np.uint8)
        np.add.at(counts, sample_top_indices[:, :top_k].ravel(), 1)
        single_top = set(np.argpartition(single_flat, -top_k)[-top_k:])
        mean_top = set(np.argpartition(mean_flat, -top_k)[-top_k:])
        for rule in stable_frequency_rules:
            min_count = int(np.ceil(rule * n_samples))
            stable_idx = np.flatnonzero(counts >= min_count)
            stable_set = set(map(int, stable_idx))
            rows.append(
                {
                    "time_index": time_index,
                    "time_label": time_label,
                    "rule_type": "top_k_frequency",
                    "reference_value": top_k,
                    "frequency_rule": float(rule),
                    "stable_edge_count": int(stable_idx.size),
                    "single_overlap_count": len(stable_set & single_top),
                    "single_overlap_fraction": len(stable_set & single_top) / top_k if top_k else np.nan,
                    "bootstrap_mean_overlap_count": len(stable_set & mean_top),
                    "bootstrap_mean_overlap_fraction": len(stable_set & mean_top) / top_k if top_k else np.nan,
                }
            )
            if stable_idx.size:
                ordered = stable_idx[np.argsort(mean_flat[stable_idx])[::-1][:top_n_details]]
                source_idx, sink_idx = np.unravel_index(ordered, single_matrix.shape)
                for source, sink, flat_idx in zip(source_idx, sink_idx, ordered, strict=True):
                    detail_rows.append(
                        {
                            "time_index": time_index,
                            "time_label": time_label,
                            "rule_type": "top_k_frequency",
                            "reference_value": top_k,
                            "frequency_rule": float(rule),
                            "source": int(source),
                            "sink": int(sink),
                            "bootstrap_frequency": int(counts[flat_idx]),
                            "single_value": float(single_matrix[source, sink]),
                            "bootstrap_mean_value": float(bootstrap_mean[source, sink]),
                        }
                    )

    for threshold, counts_matrix in threshold_presence_counts.items():
        counts = counts_matrix.ravel()
        single_active = set(np.flatnonzero(single_flat > threshold))
        mean_active = set(np.flatnonzero(mean_flat > threshold))
        for rule in stable_frequency_rules:
            min_count = int(np.ceil(rule * n_samples))
            stable_idx = np.flatnonzero(counts >= min_count)
            stable_set = set(map(int, stable_idx))
            rows.append(
                {
                    "time_index": time_index,
                    "time_label": time_label,
                    "rule_type": "threshold_frequency",
                    "reference_value": float(threshold),
                    "frequency_rule": float(rule),
                    "stable_edge_count": int(stable_idx.size),
                    "single_overlap_count": len(stable_set & single_active),
                    "single_overlap_fraction": len(stable_set & single_active) / len(single_active) if single_active else np.nan,
                    "bootstrap_mean_overlap_count": len(stable_set & mean_active),
                    "bootstrap_mean_overlap_fraction": len(stable_set & mean_active) / len(mean_active) if mean_active else np.nan,
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(detail_rows)


def compute_backbone_summary(
    single_matrix: np.ndarray,
    bootstrap_mean: np.ndarray,
    row_sums_single: np.ndarray,
    row_sums_mean: np.ndarray,
    thresholds: Iterable[float],
    local_share_thresholds: Iterable[float],
    time_index: int,
    time_label: str,
) -> pd.DataFrame:
    """Use a simple row-share backbone proxy."""
    rows: list[dict[str, float | int | str]] = []
    for threshold in thresholds:
        single_mask = single_matrix > threshold
        mean_mask = bootstrap_mean > threshold
        for share in local_share_thresholds:
            single_share_mask = np.zeros(single_matrix.shape, dtype=bool)
            mean_share_mask = np.zeros(bootstrap_mean.shape, dtype=bool)
            valid_single = row_sums_single > 0
            valid_mean = row_sums_mean > 0
            single_share_mask[valid_single] = (
                single_matrix[valid_single] / row_sums_single[valid_single, None]
            ) >= share
            mean_share_mask[valid_mean] = (
                bootstrap_mean[valid_mean] / row_sums_mean[valid_mean, None]
            ) >= share
            single_backbone = single_mask & single_share_mask
            mean_backbone = mean_mask & mean_share_mask
            overlap = int(np.sum(single_backbone & mean_backbone))
            union = int(np.sum(single_backbone | mean_backbone))
            rows.append(
                {
                    "time_index": time_index,
                    "time_label": time_label,
                    "threshold": float(threshold),
                    "local_share": float(share),
                    "single_backbone_edges": int(np.sum(single_backbone)),
                    "bootstrap_mean_backbone_edges": int(np.sum(mean_backbone)),
                    "overlap_count": overlap,
                    "overlap_fraction_single": overlap / int(np.sum(single_backbone)) if np.sum(single_backbone) else np.nan,
                    "overlap_fraction_mean": overlap / int(np.sum(mean_backbone)) if np.sum(mean_backbone) else np.nan,
                    "jaccard_index": overlap / union if union else np.nan,
                }
            )
    return pd.DataFrame(rows)


def compute_system_coverage_summary(
    single_metrics: dict[str, float],
    bootstrap_mean_metrics: dict[str, float],
    bootstrap_sample_metrics: dict[str, np.ndarray],
    time_index: int,
    time_label: str,
) -> pd.DataFrame:
    """Compare system-level single and bootstrap-mean values against sample intervals."""
    rows: list[dict[str, float | int | str | bool]] = []
    for metric_name, sample_values in bootstrap_sample_metrics.items():
        sample_values = np.asarray(sample_values, dtype=np.float64)
        q_low = float(np.nanquantile(sample_values, 0.025))
        q_high = float(np.nanquantile(sample_values, 0.975))
        for reference_name, metric_dict in (("single", single_metrics), ("bootstrap_mean", bootstrap_mean_metrics)):
            value = float(metric_dict.get(metric_name, np.nan))
            rows.append(
                {
                    "time_index": time_index,
                    "time_label": time_label,
                    "metric_name": metric_name,
                    "reference": reference_name,
                    "value": value,
                    "bootstrap_q_low": q_low,
                    "bootstrap_q_high": q_high,
                    "outside_ci": bool((value < q_low) or (value > q_high)),
                }
            )
    return pd.DataFrame(rows)


def compute_discrepancy_vs_variability(
    sampled_single: np.ndarray,
    sampled_mean: np.ndarray,
    sampled_bootstrap: np.ndarray,
    row_sums_single: np.ndarray,
    row_sums_mean: np.ndarray,
    row_sums_boot: np.ndarray,
    col_sums_single: np.ndarray,
    col_sums_mean: np.ndarray,
    col_sums_boot: np.ndarray,
    sample_top_indices: np.ndarray,
    single_matrix: np.ndarray,
    bootstrap_mean: np.ndarray,
    top_k_values: Iterable[int],
    time_index: int,
    time_label: str,
) -> pd.DataFrame:
    """Compare single-vs-mean discrepancy to within-bootstrap variability."""
    valid_edges = np.isfinite(sampled_single) & np.isfinite(sampled_mean) & np.all(np.isfinite(sampled_bootstrap), axis=0)
    sampled_single = sampled_single[valid_edges]
    sampled_mean = sampled_mean[valid_edges]
    sampled_bootstrap = sampled_bootstrap[:, valid_edges]

    pairwise_corr = []
    pairwise_rmse = []
    for left in range(sampled_bootstrap.shape[0]):
        for right in range(left + 1, sampled_bootstrap.shape[0]):
            pairwise_corr.append(_safe_pearson(sampled_bootstrap[left], sampled_bootstrap[right]))
            pairwise_rmse.append(float(np.sqrt(np.mean((sampled_bootstrap[left] - sampled_bootstrap[right]) ** 2))))
    pairwise_corr = np.asarray(pairwise_corr, dtype=np.float64)
    pairwise_rmse = np.asarray(pairwise_rmse, dtype=np.float64)

    row_sample_vs_mean = np.array(
        [float(stats.spearmanr(row_sums_mean, sample).statistic) for sample in row_sums_boot],
        dtype=np.float64,
    )
    col_sample_vs_mean = np.array(
        [float(stats.spearmanr(col_sums_mean, sample).statistic) for sample in col_sums_boot],
        dtype=np.float64,
    )

    rows = [
        {
            "time_index": time_index,
            "time_label": time_label,
            "metric_name": "sampled_edge_pearson",
            "single_vs_bootstrap_mean": _safe_pearson(sampled_single, sampled_mean),
            "within_bootstrap_median": float(np.nanmedian(pairwise_corr)),
            "within_bootstrap_q025": float(np.nanquantile(pairwise_corr, 0.025)),
            "within_bootstrap_q975": float(np.nanquantile(pairwise_corr, 0.975)),
        },
        {
            "time_index": time_index,
            "time_label": time_label,
            "metric_name": "sampled_edge_rmse",
            "single_vs_bootstrap_mean": float(np.sqrt(np.mean((sampled_single - sampled_mean) ** 2))),
            "within_bootstrap_median": float(np.nanmedian(pairwise_rmse)),
            "within_bootstrap_q025": float(np.nanquantile(pairwise_rmse, 0.025)),
            "within_bootstrap_q975": float(np.nanquantile(pairwise_rmse, 0.975)),
        },
        {
            "time_index": time_index,
            "time_label": time_label,
            "metric_name": "row_rank_spearman",
            "single_vs_bootstrap_mean": float(stats.spearmanr(row_sums_single, row_sums_mean).statistic),
            "within_bootstrap_median": float(np.nanmedian(row_sample_vs_mean)),
            "within_bootstrap_q025": float(np.nanquantile(row_sample_vs_mean, 0.025)),
            "within_bootstrap_q975": float(np.nanquantile(row_sample_vs_mean, 0.975)),
        },
        {
            "time_index": time_index,
            "time_label": time_label,
            "metric_name": "col_rank_spearman",
            "single_vs_bootstrap_mean": float(stats.spearmanr(col_sums_single, col_sums_mean).statistic),
            "within_bootstrap_median": float(np.nanmedian(col_sample_vs_mean)),
            "within_bootstrap_q025": float(np.nanquantile(col_sample_vs_mean, 0.025)),
            "within_bootstrap_q975": float(np.nanquantile(col_sample_vs_mean, 0.975)),
        },
    ]

    single_flat = single_matrix.ravel()
    mean_flat = bootstrap_mean.ravel()
    for top_k in sorted(set(int(value) for value in top_k_values)):
        single_top = set(np.argpartition(single_flat, -top_k)[-top_k:])
        mean_top = set(np.argpartition(mean_flat, -top_k)[-top_k:])
        sample_overlap = []
        for sample_idx in range(sample_top_indices.shape[0]):
            sample_top = set(map(int, sample_top_indices[sample_idx, :top_k]))
            sample_overlap.append(len(mean_top & sample_top) / top_k)
        rows.append(
            {
                "time_index": time_index,
                "time_label": time_label,
                "metric_name": f"top_{top_k}_overlap",
                "single_vs_bootstrap_mean": len(single_top & mean_top) / top_k,
                "within_bootstrap_median": float(np.nanmedian(sample_overlap)),
                "within_bootstrap_q025": float(np.nanquantile(sample_overlap, 0.025)),
                "within_bootstrap_q975": float(np.nanquantile(sample_overlap, 0.975)),
            }
        )
    return pd.DataFrame(rows)


def sample_topk_tracker_update(
    current_values: np.ndarray,
    current_indices: np.ndarray,
    block_values: np.ndarray,
    global_block_indices: np.ndarray,
    keep_k: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Update per-sample exact top-k edges using a new block."""
    flat_block = block_values.reshape(block_values.shape[0], -1)
    local_k = min(keep_k, flat_block.shape[1])
    part_idx = np.argpartition(flat_block, -local_k, axis=1)[:, -local_k:]
    part_values = np.take_along_axis(flat_block, part_idx, axis=1)
    part_indices = global_block_indices[part_idx]

    merged_values = np.concatenate([current_values, part_values], axis=1)
    merged_indices = np.concatenate([current_indices, part_indices], axis=1)
    keep_idx = np.argpartition(merged_values, -keep_k, axis=1)[:, -keep_k:]
    next_values = np.take_along_axis(merged_values, keep_idx, axis=1)
    next_indices = np.take_along_axis(merged_indices, keep_idx, axis=1)
    sort_idx = np.argsort(next_values, axis=1)[:, ::-1]
    return (
        np.take_along_axis(next_values, sort_idx, axis=1),
        np.take_along_axis(next_indices, sort_idx, axis=1),
    )
