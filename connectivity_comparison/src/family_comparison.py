"""Family-level comparison utilities for multiple connectivity products."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .high_level_metrics import (
    _bridge_components,
    _cluster_labels,
    _safe_pearson,
    build_directed_profile_features,
    choose_cluster_count,
)
from .extended_metrics import adjusted_rand_index
from .io_utils import open_connectivity_dataset, prepare_connectivity_array
from .metadata import infer_distance_matrix


def _family_name_from_single_path(path: Path) -> str:
    stem = path.stem
    if stem.startswith("connectivity_"):
        stem = stem[len("connectivity_") :]
    if stem.endswith("_single"):
        stem = stem[: -len("_single")]
    return stem


def discover_family_datasets(
    dataset_dir: Path,
    primary_bootstrap_path: Path,
    primary_single_path: Path,
) -> list[dict[str, Any]]:
    """Discover single and bootstrap NetCDF files by family name."""
    family_specs: dict[str, dict[str, Any]] = {}

    for single_path in sorted(dataset_dir.glob("connectivity_*_single.nc")):
        family_name = _family_name_from_single_path(single_path)
        family_specs.setdefault(family_name, {})["family"] = family_name
        family_specs[family_name]["single_path"] = single_path

    for bootstrap_path in sorted(dataset_dir.glob("connectivity_*.nc")):
        if bootstrap_path.name.endswith("_single.nc"):
            continue
        family_name = bootstrap_path.stem[len("connectivity_") :]
        family_specs.setdefault(family_name, {})["family"] = family_name
        family_specs[family_name]["bootstrap_path"] = bootstrap_path

    primary_family = _family_name_from_single_path(primary_single_path)
    family_specs.setdefault(primary_family, {})["family"] = primary_family
    family_specs[primary_family]["single_path"] = primary_single_path
    family_specs[primary_family]["bootstrap_path"] = primary_bootstrap_path

    for spec in family_specs.values():
        spec["has_single"] = "single_path" in spec
        spec["has_bootstrap"] = "bootstrap_path" in spec
    return sorted(family_specs.values(), key=lambda item: item["family"])


def summarize_single_family(
    family_name: str,
    single_path: Path,
    metadata_context: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    """Build a time-mean structural summary for one family single estimate."""
    family_seed = sum((index + 1) * ord(char) for index, char in enumerate(family_name))
    dataset = open_connectivity_dataset(single_path)
    data_array = prepare_connectivity_array(
        dataset,
        var_name=str(config["connectivity_var"]),
        treatment_index=int(config["treatment_index"]),
    )
    n_time = int(data_array.sizes["time"])
    n_nodes = int(data_array.sizes["source"])
    anchor_onehot = np.asarray(metadata_context["community_anchor_onehot"], dtype=np.float32)
    n_anchor = anchor_onehot.shape[1]

    distance_matrix, _ = infer_distance_matrix(
        dataset,
        variable_candidates=list(config["distance_var_candidates"]),
    )
    long_distance_threshold = float(config["distance_class_edges_km"][-2])
    long_distance_mask = (
        np.asarray(distance_matrix > long_distance_threshold, dtype=bool)
        if distance_matrix is not None
        else None
    )

    row_sums_total = np.zeros(n_nodes, dtype=np.float64)
    col_sums_total = np.zeros(n_nodes, dtype=np.float64)
    out_anchor_total = np.zeros((n_nodes, n_anchor), dtype=np.float64)
    in_anchor_total = np.zeros((n_nodes, n_anchor), dtype=np.float64)
    anchor_matrix_total = np.zeros((n_anchor, n_anchor), dtype=np.float64)
    long_distance_fraction_values: list[float] = []

    for time_index in range(n_time):
        matrix = np.asarray(data_array.isel(time=time_index, sample=0).values, dtype=np.float32)
        matrix = np.nan_to_num(matrix, nan=0.0)
        row_sums = np.sum(matrix, axis=1, dtype=np.float64)
        col_sums = np.sum(matrix, axis=0, dtype=np.float64)
        row_sums_total += row_sums
        col_sums_total += col_sums
        out_anchor_total += matrix @ anchor_onehot
        in_anchor_total += matrix.T @ anchor_onehot
        anchor_matrix_total += anchor_onehot.T @ matrix @ anchor_onehot
        if long_distance_mask is not None:
            total = float(np.sum(matrix))
            long_fraction = float(np.sum(matrix[long_distance_mask]) / total) if total > 0 else np.nan
            long_distance_fraction_values.append(long_fraction)

    row_sums_mean = row_sums_total / max(n_time, 1)
    col_sums_mean = col_sums_total / max(n_time, 1)
    out_anchor_mean = out_anchor_total / max(n_time, 1)
    in_anchor_mean = in_anchor_total / max(n_time, 1)
    anchor_matrix_mean = anchor_matrix_total / max(n_time, 1)

    features = build_directed_profile_features(
        out_anchor_mean,
        in_anchor_mean,
        row_sums_mean,
        col_sums_mean,
    )
    selected_k = choose_cluster_count(
        features,
        min_clusters=int(config["community_min_clusters"]),
        max_clusters=int(config["community_max_clusters"]),
        seed=int(config["random_seed"]) + family_seed % 10_000,
    )
    community_labels = _cluster_labels(
        features,
        n_clusters=selected_k,
        seed=int(config["random_seed"]) + family_seed % 10_000 + 1,
    )
    bridge_score = _bridge_components(out_anchor_mean, in_anchor_mean)["bridge_score"]

    anchor_vector = anchor_matrix_mean.ravel()
    anchor_total = float(np.sum(anchor_vector))
    if anchor_total > 0:
        anchor_vector = anchor_vector / anchor_total

    summary_vector = np.concatenate(
        [
            anchor_vector,
            np.quantile(row_sums_mean, [0.5, 0.9, 0.99]),
            np.quantile(col_sums_mean, [0.5, 0.9, 0.99]),
            np.quantile(bridge_score, [0.5, 0.9, 0.99]),
            np.array(
                [
                    np.nanmean(long_distance_fraction_values) if long_distance_fraction_values else np.nan,
                ],
                dtype=np.float64,
            ),
        ]
    )

    dataset.close()
    return {
        "family": family_name,
        "row_sums_mean": row_sums_mean,
        "col_sums_mean": col_sums_mean,
        "bridge_score": bridge_score,
        "community_labels": community_labels,
        "selected_n_communities": int(selected_k),
        "summary_vector": summary_vector,
        "source_top_mean_strength": float(np.mean(np.sort(row_sums_mean)[-int(config["top_k_nodes"]):])),
        "sink_top_mean_strength": float(np.mean(np.sort(col_sums_mean)[-int(config["top_k_nodes"]):])),
        "bridge_top_mean_score": float(np.mean(np.sort(bridge_score)[-int(config["bridge_top_k"]):])),
        "long_distance_fraction_mean": float(np.nanmean(long_distance_fraction_values)) if long_distance_fraction_values else np.nan,
    }


def build_family_comparison(
    family_specs: list[dict[str, Any]],
    metadata_context: dict[str, Any],
    config: dict[str, Any],
    primary_family_name: str,
    primary_uncertainty_df: pd.DataFrame,
    primary_bridge_time_df: pd.DataFrame,
    primary_community_df: pd.DataFrame,
    primary_distance_class_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compare time-mean family structure using all available single datasets."""
    summaries = [
        summarize_single_family(
            family_name=str(spec["family"]),
            single_path=Path(spec["single_path"]),
            metadata_context=metadata_context,
            config=config,
        )
        for spec in family_specs
        if spec.get("has_single")
    ]
    if not summaries:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    top_k = int(config["family_top_k"])
    summary_lookup = {item["family"]: item for item in summaries}
    family_rows = []
    for spec in family_specs:
        family_name = str(spec["family"])
        if family_name not in summary_lookup:
            continue
        item = summary_lookup[family_name]
        family_rows.append(
            {
                "family": family_name,
                "has_bootstrap": bool(spec.get("has_bootstrap", False)),
                "selected_n_communities": item["selected_n_communities"],
                "source_top_mean_strength": item["source_top_mean_strength"],
                "sink_top_mean_strength": item["sink_top_mean_strength"],
                "bridge_top_mean_score": item["bridge_top_mean_score"],
                "long_distance_fraction_mean": item["long_distance_fraction_mean"],
            }
        )
    family_summary_df = pd.DataFrame(family_rows)

    pair_rows = []
    overlap_rows = []
    family_names = [item["family"] for item in summaries]
    for left_idx, left_name in enumerate(family_names):
        left = summary_lookup[left_name]
        for right_name in family_names[left_idx + 1 :]:
            right = summary_lookup[right_name]
            pair_rows.append(
                {
                    "family_a": left_name,
                    "family_b": right_name,
                    "source_strength_spearman": float(stats.spearmanr(left["row_sums_mean"], right["row_sums_mean"]).statistic),
                    "source_strength_pearson": _safe_pearson(left["row_sums_mean"], right["row_sums_mean"]),
                    "sink_strength_spearman": float(stats.spearmanr(left["col_sums_mean"], right["col_sums_mean"]).statistic),
                    "sink_strength_pearson": _safe_pearson(left["col_sums_mean"], right["col_sums_mean"]),
                    "bridge_score_spearman": float(stats.spearmanr(left["bridge_score"], right["bridge_score"]).statistic),
                    "bridge_score_pearson": _safe_pearson(left["bridge_score"], right["bridge_score"]),
                    "community_ari": float(adjusted_rand_index(left["community_labels"], right["community_labels"]))
                    if left["community_labels"].shape == right["community_labels"].shape
                    else np.nan,
                    "long_distance_fraction_gap": float(left["long_distance_fraction_mean"] - right["long_distance_fraction_mean"]),
                }
            )
            for metric_name, left_values, right_values in (
                ("source_strength", left["row_sums_mean"], right["row_sums_mean"]),
                ("sink_strength", left["col_sums_mean"], right["col_sums_mean"]),
                ("bridge_score", left["bridge_score"], right["bridge_score"]),
            ):
                left_top = set(np.argsort(left_values)[-top_k:])
                right_top = set(np.argsort(right_values)[-top_k:])
                overlap_rows.append(
                    {
                        "family_a": left_name,
                        "family_b": right_name,
                        "metric_name": metric_name,
                        "top_k": top_k,
                        "overlap_fraction": len(left_top & right_top) / float(top_k),
                    }
                )
    pairwise_similarity_df = pd.DataFrame(pair_rows)
    topk_overlap_df = pd.DataFrame(overlap_rows)

    vector_matrix = np.vstack([summary_lookup[name]["summary_vector"] for name in family_names])
    centered = vector_matrix - np.nanmean(vector_matrix, axis=0, keepdims=True)
    centered = np.nan_to_num(centered, nan=0.0)
    if centered.shape[0] >= 2:
        _, singular_values, vt = np.linalg.svd(centered, full_matrices=False)
        coords = centered @ vt[:2].T
        explained = singular_values**2
        explained = explained / explained.sum() if explained.sum() > 0 else np.array([np.nan, np.nan])
    else:
        coords = np.zeros((centered.shape[0], 2), dtype=np.float64)
        explained = np.array([np.nan, np.nan], dtype=np.float64)
    ordination_df = pd.DataFrame(
        {
            "family": family_names,
            "pc1": coords[:, 0],
            "pc2": coords[:, 1],
            "pc1_explained_ratio": explained[0] if explained.size > 0 else np.nan,
            "pc2_explained_ratio": explained[1] if explained.size > 1 else np.nan,
        }
    )

    uncertainty_rows = []
    for spec in family_specs:
        family_name = str(spec["family"])
        row = {
            "family": family_name,
            "has_bootstrap": bool(spec.get("has_bootstrap", False)),
            "median_edge_cv_nonzero": np.nan,
            "fraction_edge_values_outside_ci": np.nan,
            "median_bridge_cv": np.nan,
            "mean_consensus_probability": np.nan,
            "long_distance_fraction_mean": summary_lookup.get(family_name, {}).get("long_distance_fraction_mean", np.nan),
        }
        if family_name == primary_family_name and not primary_uncertainty_df.empty:
            row["median_edge_cv_nonzero"] = float(primary_uncertainty_df["median_edge_cv_nonzero"].median())
            row["fraction_edge_values_outside_ci"] = float(primary_uncertainty_df["fraction_edge_values_outside_ci"].median())
        if family_name == primary_family_name and not primary_bridge_time_df.empty:
            row["median_bridge_cv"] = float(primary_bridge_time_df["median_bridge_cv"].median())
        if family_name == primary_family_name and not primary_community_df.empty:
            row["mean_consensus_probability"] = float(primary_community_df["mean_consensus_probability"].median())
        if family_name == primary_family_name and not primary_distance_class_df.empty:
            long_rows = primary_distance_class_df[primary_distance_class_df["distance_class"] == "long"]
            if not long_rows.empty:
                row["long_distance_fraction_mean"] = float(long_rows["bootstrap_mean_fraction_of_total"].median())
        uncertainty_rows.append(row)
    uncertainty_df = pd.DataFrame(uncertainty_rows)
    return family_summary_df, pairwise_similarity_df, topk_overlap_df, ordination_df, uncertainty_df
