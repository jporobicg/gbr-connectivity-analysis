"""Metadata and ancillary-data helpers for connectivity comparison analyses."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import xarray as xr


def _sanitize_label(value: object, prefix: str, index: int) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return f"{prefix}_{index + 1:02d}"
    text = str(value).strip()
    if not text:
        return f"{prefix}_{index + 1:02d}"
    safe = text.replace("/", "_").replace(" ", "_")
    return safe


def _encode_groups(values: pd.Series, prefix: str) -> tuple[np.ndarray, list[str], np.ndarray]:
    filled = values.copy()
    missing = filled.isna() | (filled.astype(str).str.strip() == "")
    if missing.any():
        filled.loc[missing] = [f"{prefix}_{index + 1:02d}" for index in np.flatnonzero(missing)]

    categories = pd.Index(pd.unique(filled))
    labels = [_sanitize_label(value, prefix, index) for index, value in enumerate(categories)]
    mapping = dict(zip(categories, range(len(categories)), strict=True))
    codes = filled.map(mapping).to_numpy(dtype=np.int32, na_value=-1)
    onehot = np.zeros((len(codes), len(labels)), dtype=np.float32)
    valid = codes >= 0
    onehot[np.flatnonzero(valid), codes[valid]] = 1.0
    return codes, labels, onehot


def load_reef_metadata(
    path: Path | None,
    n_nodes: int,
    community_anchor_column: str,
    region_column: str,
) -> dict[str, Any]:
    """Load reef metadata and derive reusable grouping arrays."""
    notes: list[str] = []
    fallback = pd.DataFrame({"node_id": np.arange(n_nodes, dtype=np.int32)})

    if path is None or not Path(path).exists():
        notes.append("Reef metadata file was not found; using node-index fallbacks for regions and anchors.")
        fallback["LAT"] = np.nan
        fallback["LON"] = np.nan
        fallback["community_anchor"] = [f"anchor_{index // max(1, n_nodes // 11):02d}" for index in range(n_nodes)]
        fallback["region_name"] = [f"region_{index // max(1, n_nodes // 5):02d}" for index in range(n_nodes)]
        metadata = fallback
    else:
        metadata = pd.read_csv(Path(path))
        metadata = metadata.reset_index(drop=True)
        if len(metadata) != n_nodes:
            notes.append(
                "Reef metadata length does not match the connectivity matrix; truncating or padding by node index order."
            )
            metadata = metadata.iloc[: min(len(metadata), n_nodes)].copy()
            if len(metadata) < n_nodes:
                pad = pd.DataFrame({"node_id": np.arange(len(metadata), n_nodes, dtype=np.int32)})
                metadata = pd.concat([metadata, pad], ignore_index=True)

        metadata["node_id"] = np.arange(n_nodes, dtype=np.int32)

        if community_anchor_column not in metadata:
            notes.append(
                f"Community anchor column '{community_anchor_column}' was not found; falling back to '{region_column}' or node bins."
            )
            if region_column in metadata:
                metadata["community_anchor"] = metadata[region_column]
            else:
                metadata["community_anchor"] = [f"anchor_{index // max(1, n_nodes // 11):02d}" for index in range(n_nodes)]
        else:
            metadata["community_anchor"] = metadata[community_anchor_column]

        if region_column not in metadata:
            notes.append(f"Region column '{region_column}' was not found; falling back to community anchors.")
            metadata["region_name"] = metadata["community_anchor"]
        else:
            metadata["region_name"] = metadata[region_column]

        if "LAT" not in metadata:
            notes.append("Reef latitude column 'LAT' was not found in the metadata.")
            metadata["LAT"] = np.nan
        if "LON" not in metadata:
            notes.append("Reef longitude column 'LON' was not found in the metadata.")
            metadata["LON"] = np.nan

    community_codes, community_labels, community_onehot = _encode_groups(
        metadata["community_anchor"].astype("object"),
        prefix="anchor",
    )
    region_codes, region_labels, region_onehot = _encode_groups(
        metadata["region_name"].astype("object"),
        prefix="region",
    )

    metadata["community_anchor_label"] = [community_labels[max(code, 0)] if code >= 0 else None for code in community_codes]
    metadata["region_label"] = [region_labels[max(code, 0)] if code >= 0 else None for code in region_codes]

    return {
        "reef_metadata": metadata[
            [
                "node_id",
                "LAT",
                "LON",
                "community_anchor_label",
                "region_label",
            ]
        ].copy(),
        "community_anchor_codes": community_codes,
        "community_anchor_labels": community_labels,
        "community_anchor_onehot": community_onehot,
        "region_codes": region_codes,
        "region_labels": region_labels,
        "region_onehot": region_onehot,
        "notes": notes,
        "community_anchor_column": community_anchor_column,
        "region_column": region_column,
    }


def _distance_candidate_score(values: np.ndarray) -> tuple[float, dict[str, float]]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return -np.inf, {"max": np.nan, "non_integer_fraction": np.nan, "n_unique": np.nan}

    sample = finite[: min(200_000, finite.size)]
    n_unique = float(len(np.unique(sample)))
    non_integer_fraction = float(np.mean(~np.isclose(sample, np.round(sample))))
    max_value = float(np.max(sample))
    score = n_unique + (50.0 if max_value > 100 else 0.0) + (50.0 * non_integer_fraction)
    return score, {
        "max": max_value,
        "non_integer_fraction": non_integer_fraction,
        "n_unique": n_unique,
    }


def infer_distance_matrix(
    dataset: xr.Dataset,
    variable_candidates: list[str],
) -> tuple[np.ndarray | None, dict[str, Any]]:
    """Choose the most plausible continuous reef-to-reef distance matrix."""
    chosen_name = None
    chosen_values: np.ndarray | None = None
    chosen_score = -np.inf
    diagnostics: dict[str, Any] = {"candidates": {}}

    for name in variable_candidates:
        if name not in dataset:
            continue
        data = dataset[name]
        if tuple(data.dims) != ("source", "sink"):
            continue
        values = np.asarray(data.values, dtype=np.float32)
        score, stats = _distance_candidate_score(values)
        stats.update({"attrs": dict(data.attrs)})
        diagnostics["candidates"][name] = stats
        if score > chosen_score:
            chosen_name = name
            chosen_values = values
            chosen_score = score

    diagnostics["selected_variable"] = chosen_name
    if chosen_name == "direction" and "distance" in diagnostics["candidates"]:
        diagnostics["note"] = (
            "The ancillary variable named 'direction' was selected as the continuous distance matrix because "
            "'distance' appears to contain discrete distance bins in the current NetCDF files."
        )
    return chosen_values, diagnostics


def append_node_metadata(df: pd.DataFrame, metadata: pd.DataFrame | None) -> pd.DataFrame:
    """Attach node coordinates and labels when the table contains `node_id`."""
    if metadata is None or metadata.empty or df.empty or "node_id" not in df:
        return df
    extra_cols = [column for column in metadata.columns if column != "node_id"]
    return df.merge(metadata[["node_id", *extra_cols]], on="node_id", how="left")
