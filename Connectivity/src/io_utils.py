"""Input/output helpers for the connectivity comparison workflow."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import xarray as xr

from .preprocessing import format_time_label


def open_connectivity_dataset(path: Path) -> xr.Dataset:
    """Open a NetCDF dataset with masked values decoded."""
    if not path.exists():
        raise FileNotFoundError(f"NetCDF file not found: {path}")
    return xr.open_dataset(path, mask_and_scale=True)


def prepare_connectivity_array(
    dataset: xr.Dataset,
    var_name: str = "connectivity",
    treatment_index: int = 0,
) -> xr.DataArray:
    """Return the connectivity data in a consistent dimension order."""
    if var_name not in dataset:
        raise KeyError(f"Variable '{var_name}' not found. Available variables: {list(dataset.data_vars)}")

    data = dataset[var_name]

    if "treatment" in data.dims:
        treatment_size = int(data.sizes["treatment"])
        if treatment_size <= treatment_index:
            raise IndexError(
                f"Requested treatment index {treatment_index} but dataset only has {treatment_size} treatment(s)."
            )
        data = data.isel(treatment=treatment_index, drop=True)

    if "sample" not in data.dims:
        data = data.expand_dims(sample=[0])

    expected_dims = ["time", "source", "sink", "sample"]
    missing = [dim for dim in expected_dims if dim not in data.dims]
    if missing:
        raise ValueError(f"Connectivity variable is missing expected dimensions: {missing}")

    return data.transpose("time", "source", "sink", "sample")


def validate_datasets(bootstrap: xr.DataArray, single: xr.DataArray) -> dict[str, Any]:
    """Check that the two datasets can be compared directly."""
    for dim in ("time", "source", "sink"):
        if int(bootstrap.sizes[dim]) != int(single.sizes[dim]):
            raise ValueError(
                f"Dimension mismatch for '{dim}': bootstrap={bootstrap.sizes[dim]}, single={single.sizes[dim]}"
            )
    if int(single.sizes["sample"]) != 1:
        raise ValueError(f"Single dataset should have one sample, found {single.sizes['sample']}")
    if int(bootstrap.sizes["sample"]) < 2:
        raise ValueError("Bootstrap dataset should have more than one sample.")

    bootstrap_times = np.asarray(bootstrap["time"].values)
    single_times = np.asarray(single["time"].values)
    if bootstrap_times.shape != single_times.shape:
        raise ValueError("Time coordinates differ between the bootstrap and single datasets.")

    warnings: list[str] = []
    for index, (boot, sing) in enumerate(zip(bootstrap_times, single_times, strict=True)):
        if (pd.isna(boot) and pd.isna(sing)) or (boot == sing):
            continue
        if pd.isna(boot) and not pd.isna(sing):
            warnings.append(
                f"time index {index} is NaT in bootstrap but {format_time_label(sing, index)} in single; using the valid label."
            )
            continue
        if pd.isna(sing) and not pd.isna(boot):
            warnings.append(
                f"time index {index} is NaT in single but {format_time_label(boot, index)} in bootstrap; using the valid label."
            )
            continue
        raise ValueError(
            "Time coordinates differ between the bootstrap and single datasets at "
            f"index {index}: bootstrap={boot}, single={sing}"
        )

    return {
        "n_time": int(bootstrap.sizes["time"]),
        "n_source": int(bootstrap.sizes["source"]),
        "n_sink": int(bootstrap.sizes["sink"]),
        "n_bootstrap_samples": int(bootstrap.sizes["sample"]),
        "time_warnings": warnings,
    }


def resolve_time_values(bootstrap: xr.DataArray, single: xr.DataArray) -> np.ndarray:
    """Return one time label per index, preferring non-missing coordinates."""
    bootstrap_times = np.asarray(bootstrap["time"].values)
    single_times = np.asarray(single["time"].values)
    resolved: list[object] = []
    for boot, sing in zip(bootstrap_times, single_times, strict=True):
        if not pd.isna(boot):
            resolved.append(boot)
        elif not pd.isna(sing):
            resolved.append(sing)
        else:
            resolved.append(boot)
    return np.asarray(resolved, dtype=object)


def dataset_summary(label: str, path: Path, data: xr.DataArray) -> dict[str, Any]:
    """Return a lightweight summary without scanning the full dataset."""
    time_values = np.asarray(data["time"].values)
    first_valid_index = 0
    if time_values.size:
        for index, value in enumerate(time_values):
            if not pd.isna(value):
                first_valid_index = index
                break

    preview = data.isel(time=first_valid_index, sample=0).values.astype(np.float32)
    finite = np.isfinite(preview)
    nonzero = finite & (preview != 0)

    summary = {
        "label": label,
        "path": str(path),
        "shape": tuple(int(data.sizes[dim]) for dim in data.dims),
        "dims": {dim: int(data.sizes[dim]) for dim in data.dims},
        "time_preview_index": first_valid_index,
        "time_preview_label": format_time_label(time_values[first_valid_index], first_valid_index),
        "finite_fraction_preview": float(finite.mean()),
        "zero_fraction_preview": float((preview[finite] == 0).mean()) if finite.any() else np.nan,
        "min_preview": float(np.nanmin(preview)) if finite.any() else np.nan,
        "max_preview": float(np.nanmax(preview)) if finite.any() else np.nan,
        "median_nonzero_preview": float(np.nanmedian(preview[nonzero])) if nonzero.any() else np.nan,
    }
    return summary


def ensure_output_directories(config: dict[str, Any]) -> None:
    """Create output directories if needed."""
    for key in ("output_dir", "tables_dir", "figures_dir", "intermediate_dir"):
        Path(config[key]).mkdir(parents=True, exist_ok=True)


def write_dataframe(df: pd.DataFrame, path: Path) -> None:
    """Write a dataframe to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def write_markdown(text: str, path: Path) -> None:
    """Write Markdown text to disk."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
