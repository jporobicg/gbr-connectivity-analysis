"""Weighted directed-network summaries for connectivity matrices."""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components


def gini(values: np.ndarray) -> float:
    """Compute the Gini coefficient for non-negative values."""
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    if np.all(values == 0):
        return 0.0
    values = np.sort(values)
    index = np.arange(1, values.size + 1)
    return float((2 * np.sum(index * values) / (values.size * np.sum(values))) - (values.size + 1) / values.size)


def component_metrics(
    matrix: np.ndarray,
    threshold: float,
    max_edges: int = 2_000_000,
) -> tuple[float, float]:
    """Return strongly connected component counts for thresholded links."""
    active = np.isfinite(matrix) & (matrix > threshold)
    rows, cols = np.where(active)
    n_edges = rows.size
    if n_edges == 0:
        return np.nan, np.nan
    if n_edges > max_edges:
        return np.nan, np.nan
    graph = sparse.csr_matrix((np.ones(n_edges, dtype=np.int8), (rows, cols)), shape=matrix.shape)
    n_components, labels = connected_components(graph, directed=True, connection="strong")
    component_sizes = np.bincount(labels)
    largest_component = component_sizes.max() if component_sizes.size else np.nan
    return float(n_components), float(largest_component)


def network_metric_rows(
    matrix: np.ndarray,
    row_sums: np.ndarray,
    col_sums: np.ndarray,
    thresholds: Iterable[float],
    component_thresholds: Iterable[float],
    time_index: int,
    time_label: str,
    dataset_label: str,
    max_component_edges: int = 2_000_000,
) -> pd.DataFrame:
    """Summarize a matrix as a weighted directed network."""
    finite = np.isfinite(matrix)
    valid_edge_count = int(finite.sum())
    diagonal = np.diag(matrix)
    diagonal_sum = float(np.nansum(diagonal))
    total_weight = float(np.nansum(matrix))
    outgoing_gini = gini(row_sums)
    incoming_gini = gini(col_sums)
    nonzero_values = matrix[finite & (matrix > 0)]

    rows: list[dict[str, float | int | str]] = []
    for threshold in sorted(set(float(value) for value in thresholds)):
        active = finite & (matrix > threshold)
        active_count = int(active.sum())
        active_weights = matrix[active]
        out_degree = active.sum(axis=1).astype(np.int32)
        in_degree = active.sum(axis=0).astype(np.int32)

        row = {
            "time_index": time_index,
            "time_label": time_label,
            "dataset": dataset_label,
            "threshold": threshold,
            "n_valid_edges": valid_edge_count,
            "n_active_edges": active_count,
            "connectance": active_count / valid_edge_count if valid_edge_count else np.nan,
            "total_weight": total_weight,
            "self_recruitment_sum": diagonal_sum,
            "self_recruitment_fraction": diagonal_sum / total_weight if total_weight else np.nan,
            "outgoing_gini": outgoing_gini,
            "incoming_gini": incoming_gini,
            "mean_outgoing_strength": float(np.nanmean(row_sums)),
            "median_outgoing_strength": float(np.nanmedian(row_sums)),
            "p95_outgoing_strength": float(np.nanquantile(row_sums, 0.95)),
            "mean_incoming_strength": float(np.nanmean(col_sums)),
            "median_incoming_strength": float(np.nanmedian(col_sums)),
            "p95_incoming_strength": float(np.nanquantile(col_sums, 0.95)),
            "mean_out_degree": float(np.mean(out_degree)),
            "median_out_degree": float(np.median(out_degree)),
            "p95_out_degree": float(np.quantile(out_degree, 0.95)),
            "mean_in_degree": float(np.mean(in_degree)),
            "median_in_degree": float(np.median(in_degree)),
            "p95_in_degree": float(np.quantile(in_degree, 0.95)),
            "mean_active_weight": float(np.nanmean(active_weights)) if active_weights.size else np.nan,
            "median_active_weight": float(np.nanmedian(active_weights)) if active_weights.size else np.nan,
            "p95_active_weight": float(np.nanquantile(active_weights, 0.95)) if active_weights.size else np.nan,
            "max_active_weight": float(np.nanmax(active_weights)) if active_weights.size else np.nan,
            "global_nonzero_mean_weight": float(np.nanmean(nonzero_values)) if nonzero_values.size else np.nan,
            "global_nonzero_median_weight": float(np.nanmedian(nonzero_values)) if nonzero_values.size else np.nan,
        }

        if threshold in set(float(value) for value in component_thresholds):
            n_components, largest_component = component_metrics(
                matrix,
                threshold=threshold,
                max_edges=max_component_edges,
            )
            row["n_strong_components"] = n_components
            row["largest_strong_component"] = largest_component
        else:
            row["n_strong_components"] = np.nan
            row["largest_strong_component"] = np.nan
        rows.append(row)

    return pd.DataFrame(rows)


def top_nodes_table(
    single_values: np.ndarray,
    bootstrap_mean_values: np.ndarray,
    entity_type: str,
    time_index: int,
    time_label: str,
    top_n: int = 20,
) -> pd.DataFrame:
    """Return the strongest exporters or importers for both products."""
    single_rank = np.argsort(single_values)[::-1]
    bootstrap_rank = np.argsort(bootstrap_mean_values)[::-1]
    keep = np.unique(np.concatenate([single_rank[:top_n], bootstrap_rank[:top_n]]))

    rows: list[dict[str, float | int | str]] = []
    for node_id in keep:
        row = {
            "time_index": time_index,
            "time_label": time_label,
            "entity_type": entity_type,
            "node_id": int(node_id),
            "single_strength": float(single_values[node_id]),
            "bootstrap_mean_strength": float(bootstrap_mean_values[node_id]),
            "single_rank": int(np.where(single_rank == node_id)[0][0] + 1),
            "bootstrap_mean_rank": int(np.where(bootstrap_rank == node_id)[0][0] + 1),
        }
        rows.append(row)

    sort_col = "single_strength" if entity_type == "source" else "bootstrap_mean_strength"
    return pd.DataFrame(rows).sort_values(sort_col, ascending=False).reset_index(drop=True)
