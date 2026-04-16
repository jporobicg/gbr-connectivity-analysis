"""Preprocessing and convenience helpers for large connectivity matrices."""

from __future__ import annotations

from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd


def iter_slices(length: int, block_size: int) -> Iterator[slice]:
    """Yield consecutive slices of fixed size."""
    for start in range(0, length, block_size):
        stop = min(start + block_size, length)
        yield slice(start, stop)


def iter_matrix_blocks(n_source: int, n_sink: int, block_size: int) -> Iterator[tuple[slice, slice]]:
    """Yield source/sink block slices for dense matrices."""
    for source_slice in iter_slices(n_source, block_size):
        for sink_slice in iter_slices(n_sink, block_size):
            yield source_slice, sink_slice


def clean_array(values: np.ndarray) -> np.ndarray:
    """Return a float32 array with masked values converted to NaN."""
    if np.ma.isMaskedArray(values):
        values = values.filled(np.nan)
    return np.asarray(values, dtype=np.float32)


def coefficient_of_variation(std: np.ndarray, mean: np.ndarray, epsilon: float = 1e-12) -> np.ndarray:
    """Compute coefficient of variation while avoiding divide-by-zero."""
    cv = np.full(mean.shape, np.nan, dtype=np.float32)
    valid = np.isfinite(std) & np.isfinite(mean) & (np.abs(mean) > epsilon)
    cv[valid] = std[valid] / np.abs(mean[valid])
    return cv


def format_time_label(value: object, index: int) -> str:
    """Convert a time coordinate to a stable string label."""
    if pd.isna(value):
        return f"time_{index:02d}_NaT"
    if isinstance(value, np.datetime64):
        return pd.Timestamp(value).strftime("%Y-%m-%d")
    return str(value)


def sampled_edge_indices(
    n_source: int,
    n_sink: int,
    n_edges: int,
    rng: np.random.Generator,
    include_diagonal: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample unique source/sink index pairs."""
    total = n_source * n_sink
    n_edges = min(n_edges, total)
    flat_indices = rng.choice(total, size=n_edges, replace=False)
    source_idx = flat_indices // n_sink
    sink_idx = flat_indices % n_sink

    if include_diagonal:
        return source_idx.astype(np.int32), sink_idx.astype(np.int32)

    valid = source_idx != sink_idx
    source_idx = source_idx[valid]
    sink_idx = sink_idx[valid]
    while source_idx.size < n_edges:
        extra_needed = n_edges - source_idx.size
        extra_flat = rng.choice(total, size=extra_needed * 2, replace=False)
        extra_source = extra_flat // n_sink
        extra_sink = extra_flat % n_sink
        extra_valid = extra_source != extra_sink
        source_idx = np.concatenate([source_idx, extra_source[extra_valid]])
        sink_idx = np.concatenate([sink_idx, extra_sink[extra_valid]])
        pairs = np.unique(np.column_stack([source_idx, sink_idx]), axis=0)
        source_idx = pairs[:, 0]
        sink_idx = pairs[:, 1]
    return source_idx[:n_edges].astype(np.int32), sink_idx[:n_edges].astype(np.int32)


def assign_sampled_edges_from_block(
    source_idx: np.ndarray,
    sink_idx: np.ndarray,
    source_slice: slice,
    sink_slice: slice,
    single_block: np.ndarray,
    bootstrap_block: np.ndarray,
    sampled_single: np.ndarray,
    sampled_bootstrap: np.ndarray,
) -> None:
    """Extract sampled edge values that fall inside the current matrix block."""
    source_start = int(source_slice.start)
    source_stop = int(source_slice.stop)
    sink_start = int(sink_slice.start)
    sink_stop = int(sink_slice.stop)

    in_block = (
        (source_idx >= source_start)
        & (source_idx < source_stop)
        & (sink_idx >= sink_start)
        & (sink_idx < sink_stop)
    )
    if not np.any(in_block):
        return

    local_source = source_idx[in_block] - source_start
    local_sink = sink_idx[in_block] - sink_start
    sampled_single[in_block] = single_block[local_source, local_sink]
    sampled_bootstrap[:, in_block] = bootstrap_block[:, local_source, local_sink]


def block_average(matrix: np.ndarray, output_size: int = 120) -> np.ndarray:
    """Downsample a large matrix by averaging over coarse blocks."""
    source_groups = np.array_split(np.arange(matrix.shape[0]), output_size)
    sink_groups = np.array_split(np.arange(matrix.shape[1]), output_size)
    reduced = np.full((len(source_groups), len(sink_groups)), np.nan, dtype=np.float32)

    for i, source_group in enumerate(source_groups):
        for j, sink_group in enumerate(sink_groups):
            block = matrix[np.ix_(source_group, sink_group)]
            if np.isfinite(block).any():
                reduced[i, j] = np.nanmean(block)

    return reduced


def npz_path(directory: Path, time_index: int, stem: str) -> Path:
    """Return a stable path for intermediate NumPy archives."""
    return directory / f"{stem}_time_{time_index:02d}.npz"
