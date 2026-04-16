"""Matrix comparison metrics and bootstrap uncertainty summaries."""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats


def _safe_pearson(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2:
        return np.nan
    x_std = float(np.std(x))
    y_std = float(np.std(y))
    if x_std == 0.0 or y_std == 0.0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def _safe_spearman(
    x: np.ndarray,
    y: np.ndarray,
    max_points: int = 200000,
    rng: np.random.Generator | None = None,
) -> tuple[float, int]:
    if x.size < 2:
        return np.nan, int(x.size)
    if x.size > max_points:
        if rng is None:
            rng = np.random.default_rng(42)
        sample_idx = rng.choice(x.size, size=max_points, replace=False)
        x = x[sample_idx]
        y = y[sample_idx]
    result = stats.spearmanr(x, y, nan_policy="omit")
    return float(result.statistic), int(x.size)


def _subset_definitions(thresholds: Iterable[float], default_change_threshold: float) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = [
        {"subset": "all", "threshold": 0.0, "change_threshold": default_change_threshold},
        {"subset": "nonzero", "threshold": 0.0, "change_threshold": default_change_threshold},
    ]
    for threshold in thresholds:
        rows.append({"subset": "threshold", "threshold": float(threshold), "change_threshold": float(threshold)})
    return rows


def compute_matrix_comparison_metrics(
    single: np.ndarray,
    bootstrap_mean: np.ndarray,
    thresholds: Iterable[float],
    time_index: int,
    time_label: str,
    row_sums_single: np.ndarray,
    row_sums_mean: np.ndarray,
    col_sums_single: np.ndarray,
    col_sums_mean: np.ndarray,
    default_change_threshold: float = 1e-10,
    relative_diff_epsilon: float = 1e-12,
    spearman_max_points: int = 200000,
    random_seed: int = 42,
) -> pd.DataFrame:
    """Compare two dense matrices under multiple masking rules."""
    rng = np.random.default_rng(random_seed + time_index)
    valid = np.isfinite(single) & np.isfinite(bootstrap_mean)
    rows: list[dict[str, float | int | str]] = []

    for subset in _subset_definitions(thresholds, default_change_threshold):
        subset_name = str(subset["subset"])
        threshold = float(subset["threshold"])
        change_threshold = float(subset["change_threshold"])

        if subset_name == "all":
            mask = valid
        elif subset_name == "nonzero":
            mask = valid & ((single != 0) | (bootstrap_mean != 0))
        else:
            mask = valid & ((single > threshold) | (bootstrap_mean > threshold))

        x = single[mask]
        y = bootstrap_mean[mask]
        if x.size == 0:
            row = {
                "time_index": time_index,
                "time_label": time_label,
                "subset": subset_name,
                "threshold": threshold,
                "n_edges": 0,
            }
            rows.append(row)
            continue

        diff = x - y
        abs_diff = np.abs(diff)
        pearson_r = _safe_pearson(x, y)
        spearman_r, spearman_n = _safe_spearman(x, y, max_points=spearman_max_points, rng=rng)

        denom = float(np.sum(np.abs(y)))
        relative_abs_diff = np.nan if denom <= relative_diff_epsilon else float(np.sum(abs_diff) / denom)

        row = {
            "time_index": time_index,
            "time_label": time_label,
            "subset": subset_name,
            "threshold": threshold,
            "n_edges": int(x.size),
            "pearson_r": pearson_r,
            "spearman_r": spearman_r,
            "spearman_sample_size": spearman_n,
            "rmse": float(np.sqrt(np.mean(diff**2))),
            "mae": float(np.mean(abs_diff)),
            "relative_abs_diff": relative_abs_diff,
            "sum_abs_diff": float(np.sum(abs_diff)),
            "frobenius_norm": float(np.sqrt(np.sum(diff**2))),
            "changed_fraction": float(np.mean(abs_diff > change_threshold)),
            "single_total_weight": float(np.sum(x)),
            "bootstrap_mean_total_weight": float(np.sum(y)),
            "row_sum_pearson_r": _safe_pearson(row_sums_single, row_sums_mean),
            "row_sum_spearman_r": float(stats.spearmanr(row_sums_single, row_sums_mean).statistic),
            "row_sum_rmse": float(np.sqrt(np.mean((row_sums_single - row_sums_mean) ** 2))),
            "col_sum_pearson_r": _safe_pearson(col_sums_single, col_sums_mean),
            "col_sum_spearman_r": float(stats.spearmanr(col_sums_single, col_sums_mean).statistic),
            "col_sum_rmse": float(np.sqrt(np.mean((col_sums_single - col_sums_mean) ** 2))),
        }
        rows.append(row)

    return pd.DataFrame(rows)


def compute_topk_edge_overlap(
    single: np.ndarray,
    bootstrap_mean: np.ndarray,
    top_k_values: Iterable[int],
    time_index: int,
    time_label: str,
) -> pd.DataFrame:
    """Summarize overlap among the strongest links."""
    valid = np.isfinite(single) & np.isfinite(bootstrap_mean)
    single_flat = np.where(valid, single, -np.inf).ravel()
    mean_flat = np.where(valid, bootstrap_mean, -np.inf).ravel()
    n_total = int(np.isfinite(single_flat).sum())
    rows: list[dict[str, float | int | str]] = []

    for top_k in top_k_values:
        k = min(int(top_k), n_total)
        if k <= 0:
            continue
        single_idx = np.argpartition(single_flat, -k)[-k:]
        mean_idx = np.argpartition(mean_flat, -k)[-k:]
        single_set = set(map(int, single_idx))
        mean_set = set(map(int, mean_idx))
        overlap = single_set & mean_set
        union = single_set | mean_set
        rows.append(
            {
                "time_index": time_index,
                "time_label": time_label,
                "top_k": k,
                "overlap_count": len(overlap),
                "overlap_fraction": len(overlap) / k,
                "jaccard_index": len(overlap) / len(union) if union else np.nan,
            }
        )
    return pd.DataFrame(rows)


def compute_rank_stability(
    single_values: np.ndarray,
    bootstrap_mean_values: np.ndarray,
    top_k_values: Iterable[int],
    entity_label: str,
    time_index: int,
    time_label: str,
    top_n_changes: int = 50,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare strength rankings and identify the largest changes."""
    single_rank = pd.Series(single_values).rank(method="dense", ascending=False).to_numpy()
    mean_rank = pd.Series(bootstrap_mean_values).rank(method="dense", ascending=False).to_numpy()
    abs_rank_change = np.abs(single_rank - mean_rank)

    summary_rows: list[dict[str, float | int | str]] = []
    for top_k in top_k_values:
        single_top = set(np.argsort(single_values)[-top_k:])
        mean_top = set(np.argsort(bootstrap_mean_values)[-top_k:])
        overlap = len(single_top & mean_top)
        summary_rows.append(
            {
                "time_index": time_index,
                "time_label": time_label,
                "entity_type": entity_label,
                "top_k": int(top_k),
                "spearman_r": float(stats.spearmanr(single_values, bootstrap_mean_values).statistic),
                "pearson_r": _safe_pearson(single_values, bootstrap_mean_values),
                "overlap_count": overlap,
                "overlap_fraction": overlap / float(top_k),
                "max_rank_shift": float(abs_rank_change.max()),
                "median_rank_shift": float(np.median(abs_rank_change)),
            }
        )

    detail = pd.DataFrame(
        {
            "time_index": time_index,
            "time_label": time_label,
            "entity_type": entity_label,
            "node_id": np.arange(single_values.size, dtype=int),
            "single_strength": single_values,
            "bootstrap_mean_strength": bootstrap_mean_values,
            "single_rank": single_rank,
            "bootstrap_mean_rank": mean_rank,
            "rank_shift": mean_rank - single_rank,
            "abs_rank_shift": abs_rank_change,
            "abs_strength_difference": np.abs(single_values - bootstrap_mean_values),
        }
    ).sort_values(["abs_rank_shift", "abs_strength_difference"], ascending=[False, False])

    return pd.DataFrame(summary_rows), detail.head(top_n_changes).reset_index(drop=True)


def top_changing_links(
    single: np.ndarray,
    bootstrap_mean: np.ndarray,
    bootstrap_std: np.ndarray | None,
    bootstrap_cv: np.ndarray | None,
    q_low: np.ndarray | None,
    q_high: np.ndarray | None,
    time_index: int,
    time_label: str,
    top_n: int = 100,
    relative_diff_epsilon: float = 1e-12,
) -> pd.DataFrame:
    """Return the links with the largest absolute differences."""
    valid = np.isfinite(single) & np.isfinite(bootstrap_mean)
    diff = np.where(valid, np.abs(single - bootstrap_mean), -np.inf)
    n_keep = min(top_n, int(np.isfinite(diff).sum()))
    if n_keep == 0:
        return pd.DataFrame()

    flat_idx = np.argpartition(diff.ravel(), -n_keep)[-n_keep:]
    source_idx, sink_idx = np.unravel_index(flat_idx, diff.shape)

    rows: list[dict[str, float | int | str | bool]] = []
    for source, sink in zip(source_idx, sink_idx, strict=True):
        single_value = float(single[source, sink])
        mean_value = float(bootstrap_mean[source, sink])
        abs_diff = abs(single_value - mean_value)
        denom = max(abs(mean_value), relative_diff_epsilon)
        row = {
            "time_index": time_index,
            "time_label": time_label,
            "source": int(source),
            "sink": int(sink),
            "single_value": single_value,
            "bootstrap_mean_value": mean_value,
            "abs_difference": abs_diff,
            "relative_difference": abs_diff / denom,
        }
        if bootstrap_std is not None:
            row["bootstrap_std"] = float(bootstrap_std[source, sink])
        if bootstrap_cv is not None:
            row["bootstrap_cv"] = float(bootstrap_cv[source, sink])
        if q_low is not None and q_high is not None:
            q_low_value = float(q_low[source, sink])
            q_high_value = float(q_high[source, sink])
            row["bootstrap_q_low"] = q_low_value
            row["bootstrap_q_high"] = q_high_value
            row["single_outside_ci"] = bool((single_value < q_low_value) or (single_value > q_high_value))
        rows.append(row)

    return pd.DataFrame(rows).sort_values("abs_difference", ascending=False).reset_index(drop=True)


def summarize_single_vs_bootstrap_ci(
    single: np.ndarray,
    q_low: np.ndarray,
    q_high: np.ndarray,
    thresholds: Iterable[float],
    time_index: int,
    time_label: str,
) -> pd.DataFrame:
    """Quantify where the single matrix falls outside bootstrap intervals."""
    valid = np.isfinite(single) & np.isfinite(q_low) & np.isfinite(q_high)
    rows: list[dict[str, float | int | str]] = []

    for subset in ["all", "nonzero"]:
        if subset == "all":
            mask = valid
        else:
            mask = valid & (single != 0)
        values = single[mask]
        lower = q_low[mask]
        upper = q_high[mask]
        if values.size == 0:
            continue
        rows.append(
            {
                "time_index": time_index,
                "time_label": time_label,
                "subset": subset,
                "threshold": 0.0,
                "n_edges": int(values.size),
                "fraction_outside_ci": float(np.mean((values < lower) | (values > upper))),
                "fraction_below_ci": float(np.mean(values < lower)),
                "fraction_above_ci": float(np.mean(values > upper)),
            }
        )

    for threshold in thresholds:
        mask = valid & (single > threshold)
        values = single[mask]
        lower = q_low[mask]
        upper = q_high[mask]
        if values.size == 0:
            continue
        rows.append(
            {
                "time_index": time_index,
                "time_label": time_label,
                "subset": "threshold",
                "threshold": float(threshold),
                "n_edges": int(values.size),
                "fraction_outside_ci": float(np.mean((values < lower) | (values > upper))),
                "fraction_below_ci": float(np.mean(values < lower)),
                "fraction_above_ci": float(np.mean(values > upper)),
            }
        )
    return pd.DataFrame(rows)


def bootstrap_interval_summary(
    single_values: np.ndarray,
    bootstrap_samples: np.ndarray,
    entity_type: str,
    time_index: int,
    time_label: str,
    top_n: int = 50,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarize whether node-level single values lie inside bootstrap intervals."""
    bootstrap_mean = np.nanmean(bootstrap_samples, axis=0)
    bootstrap_std = np.nanstd(bootstrap_samples, axis=0, ddof=1)
    q_low = np.nanquantile(bootstrap_samples, 0.025, axis=0)
    q_high = np.nanquantile(bootstrap_samples, 0.975, axis=0)
    outside = (single_values < q_low) | (single_values > q_high)

    z_score = np.full(single_values.shape, np.nan, dtype=np.float64)
    valid = bootstrap_std > 0
    z_score[valid] = (single_values[valid] - bootstrap_mean[valid]) / bootstrap_std[valid]

    summary = pd.DataFrame(
        [
            {
                "time_index": time_index,
                "time_label": time_label,
                "entity_type": entity_type,
                "fraction_outside_ci": float(np.mean(outside)),
                "mean_abs_zscore": float(np.nanmean(np.abs(z_score))),
                "median_abs_zscore": float(np.nanmedian(np.abs(z_score))),
            }
        ]
    )

    detail = pd.DataFrame(
        {
            "time_index": time_index,
            "time_label": time_label,
            "entity_type": entity_type,
            "node_id": np.arange(single_values.size, dtype=int),
            "single_value": single_values,
            "bootstrap_mean": bootstrap_mean,
            "bootstrap_std": bootstrap_std,
            "bootstrap_q_low": q_low,
            "bootstrap_q_high": q_high,
            "outside_ci": outside,
            "z_score": z_score,
            "abs_difference": np.abs(single_values - bootstrap_mean),
        }
    ).sort_values(["outside_ci", "abs_difference"], ascending=[False, False])

    return summary, detail.head(top_n).reset_index(drop=True)


def summarize_sample_similarity(
    sampled_bootstrap: np.ndarray,
    sampled_single: np.ndarray,
    time_index: int,
    time_label: str,
) -> pd.DataFrame:
    """Compare within-bootstrap similarity to single-vs-bootstrap similarity."""
    valid_columns = np.all(np.isfinite(sampled_bootstrap), axis=0) & np.isfinite(sampled_single)
    sampled_bootstrap = sampled_bootstrap[:, valid_columns]
    sampled_single = sampled_single[valid_columns]
    if sampled_bootstrap.shape[1] < 2:
        return pd.DataFrame(
            [
                {
                    "time_index": time_index,
                    "time_label": time_label,
                    "n_sampled_edges": int(sampled_bootstrap.shape[1]),
                }
            ]
        )

    pairwise = np.corrcoef(sampled_bootstrap)
    upper = pairwise[np.triu_indices(pairwise.shape[0], k=1)]
    single_vs_samples = np.array(
        [_safe_pearson(sampled_single, sampled_bootstrap[idx]) for idx in range(sampled_bootstrap.shape[0])],
        dtype=np.float64,
    )
    bootstrap_mean_sampled = np.nanmean(sampled_bootstrap, axis=0)
    single_vs_bootstrap_mean = _safe_pearson(sampled_single, bootstrap_mean_sampled)

    return pd.DataFrame(
        [
            {
                "time_index": time_index,
                "time_label": time_label,
                "n_sampled_edges": int(sampled_bootstrap.shape[1]),
                "mean_pairwise_bootstrap_r": float(np.nanmean(upper)),
                "median_pairwise_bootstrap_r": float(np.nanmedian(upper)),
                "q025_pairwise_bootstrap_r": float(np.nanquantile(upper, 0.025)),
                "q975_pairwise_bootstrap_r": float(np.nanquantile(upper, 0.975)),
                "mean_single_vs_bootstrap_r": float(np.nanmean(single_vs_samples)),
                "median_single_vs_bootstrap_r": float(np.nanmedian(single_vs_samples)),
                "single_vs_bootstrap_mean_r": float(single_vs_bootstrap_mean),
            }
        ]
    )


def pca_ordination(
    bootstrap_features: np.ndarray,
    single_features: np.ndarray,
    time_index: int,
    time_label: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Project bootstrap samples and the single estimate into a 2D PCA space."""
    feature_mean = np.nanmean(bootstrap_features, axis=0, keepdims=True)
    centered_bootstrap = bootstrap_features - feature_mean
    centered_single = single_features - feature_mean[0]
    bootstrap_mean_features = np.nanmean(bootstrap_features, axis=0)
    centered_bootstrap_mean = bootstrap_mean_features - feature_mean[0]

    u, singular_values, vt = np.linalg.svd(centered_bootstrap, full_matrices=False)
    coords = u[:, :2] * singular_values[:2]
    single_coords = centered_single @ vt[:2].T
    mean_coords = centered_bootstrap_mean @ vt[:2].T

    variance = singular_values**2
    explained_ratio = variance / variance.sum() if variance.sum() > 0 else np.array([np.nan, np.nan])

    coord_df = pd.DataFrame(
        {
            "time_index": time_index,
            "time_label": time_label,
            "sample_label": [f"bootstrap_{idx:03d}" for idx in range(bootstrap_features.shape[0])] + ["bootstrap_mean", "single"],
            "sample_type": ["bootstrap"] * bootstrap_features.shape[0] + ["bootstrap_mean", "single"],
            "pc1": np.concatenate([coords[:, 0], [mean_coords[0], single_coords[0]]]),
            "pc2": np.concatenate([coords[:, 1], [mean_coords[1], single_coords[1]]]),
        }
    )
    meta_df = pd.DataFrame(
        [
            {
                "time_index": time_index,
                "time_label": time_label,
                "pc1_explained_ratio": float(explained_ratio[0]) if explained_ratio.size > 0 else np.nan,
                "pc2_explained_ratio": float(explained_ratio[1]) if explained_ratio.size > 1 else np.nan,
            }
        ]
    )
    return coord_df, meta_df
