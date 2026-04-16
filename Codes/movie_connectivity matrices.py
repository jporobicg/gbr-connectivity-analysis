#!/usr/bin/env python3
"""Create a movie of connectivity matrices by date."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

DEFAULT_INPUT = Path(
    "/home/por07g/Documents/Projects/GBR_modeling/Connectivity_analysis/datasets/connectivity_matrices/connectivity_merulinidae_single.nc"
)
DEFAULT_OUTPUT = Path(
    "/home/por07g/Documents/Projects/GBR_modeling/Connectivity_analysis/figures/connectivity_matrices_by_date.mov"
)
FILL_THRESHOLD = 1e30


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a .mov movie of connectivity matrices across all dates."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Input NetCDF path.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Output .mov file path.",
    )
    parser.add_argument(
        "--fps",
        type=float,
        default=2.0,
        help="Frames per second for the movie.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="DPI for each frame image.",
    )
    parser.add_argument(
        "--cmap",
        type=str,
        default="inferno",
        help="Matplotlib colormap name.",
    )
    return parser.parse_args()


def select_connectivity(ds: xr.Dataset) -> xr.DataArray:
    if "connectivity" not in ds.data_vars:
        raise KeyError("Dataset does not contain a 'connectivity' variable.")

    conn = ds["connectivity"]

    for dim in ("treatment", "sample"):
        if dim in conn.dims:
            conn = conn.isel({dim: 0})

    expected_dims = {"time", "source", "sink"}
    if not expected_dims.issubset(set(conn.dims)):
        raise ValueError(
            f"'connectivity' must include dims {expected_dims}, found {conn.dims}."
        )

    return conn.transpose("time", "source", "sink")


def mask_invalid(matrix: np.ndarray, fill_value: float | None) -> np.ma.MaskedArray:
    arr = np.asarray(matrix, dtype=np.float32)
    invalid = ~np.isfinite(arr) | (arr >= FILL_THRESHOLD) | (arr <= 0.0)
    if fill_value is not None:
        invalid |= arr == fill_value
    return np.ma.array(arr, mask=invalid)


def compute_color_limits(conn: xr.DataArray, fill_value: float | None) -> tuple[float, float]:
    vmin = np.inf
    vmax = -np.inf
    n_time = conn.sizes["time"]

    for i in range(n_time):
        matrix = mask_invalid(conn.isel(time=i).values, fill_value)
        if matrix.count() == 0:
            continue

        local_min = float(matrix.min())
        local_max = float(matrix.max())
        vmin = min(vmin, local_min)
        vmax = max(vmax, local_max)
        print(
            f"[{i + 1}/{n_time}] range scan: local_min={local_min:.6g}, local_max={local_max:.6g}",
            flush=True,
        )

    if not np.isfinite(vmin) or not np.isfinite(vmax):
        raise ValueError("No valid connectivity values found to define color limits.")
    if vmax <= vmin:
        vmax = vmin + 1e-12

    return vmin, vmax


def date_label(value: np.datetime64) -> str:
    return np.datetime_as_string(np.datetime64(value), unit="D")


def render_frames(
    conn: xr.DataArray,
    fill_value: float | None,
    vmin: float,
    vmax: float,
    temp_dir: Path,
    dpi: int,
    cmap: str,
) -> None:
    n_time = conn.sizes["time"]
    time_values = conn["time"].values

    for i in range(n_time):
        matrix = mask_invalid(conn.isel(time=i).values, fill_value)
        title = f"Connectivity matrix - {date_label(time_values[i])}"

        fig, ax = plt.subplots(figsize=(8, 8), constrained_layout=True)
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")
        cmap_obj = plt.get_cmap(cmap).copy()
        cmap_obj.set_bad(color="white")
        image = ax.imshow(
            matrix,
            cmap=cmap_obj,
            origin="lower",
            interpolation="nearest",
            aspect="equal",
            vmin=vmin,
            vmax=vmax,
        )

        ax.set_title(title)
        ax.set_xlabel("Sink index")
        ax.set_ylabel("Source index")
        ax.set_xticks([])
        ax.set_yticks([])

        cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Connectivity")

        frame_path = temp_dir / f"frame_{i:04d}.png"
        fig.savefig(frame_path, dpi=dpi)
        plt.close(fig)
        print(f"[{i + 1}/{n_time}] frame saved: {frame_path.name}", flush=True)


def encode_movie_with_ffmpeg(frames_dir: Path, output_path: Path, fps: float) -> None:
    ffmpeg_bin = shutil.which("ffmpeg")
    if ffmpeg_bin is None:
        raise RuntimeError(
            "ffmpeg is required but was not found in PATH. Install ffmpeg and rerun."
        )

    encoders_info = subprocess.run(
        [ffmpeg_bin, "-hide_banner", "-encoders"],
        check=True,
        capture_output=True,
        text=True,
    )
    available_encoders = encoders_info.stdout + encoders_info.stderr

    codec_args: list[str]
    if "libx264" in available_encoders:
        codec_args = ["-c:v", "libx264", "-pix_fmt", "yuv420p"]
    elif "libopenh264" in available_encoders:
        codec_args = ["-c:v", "libopenh264", "-pix_fmt", "yuv420p", "-b:v", "8M"]
    elif "mpeg4" in available_encoders:
        codec_args = ["-c:v", "mpeg4", "-q:v", "3"]
    else:
        raise RuntimeError(
            "No suitable ffmpeg encoder found (tried libx264, libopenh264, mpeg4)."
        )

    cmd = [
        ffmpeg_bin,
        "-y",
        "-framerate",
        str(fps),
        "-i",
        str(frames_dir / "frame_%04d.png"),
        *codec_args,
        str(output_path),
    ]
    subprocess.run(cmd, check=True)


def main() -> None:
    args = parse_args()

    if not args.input.exists():
        raise FileNotFoundError(f"Input file not found: {args.input}")
    if shutil.which("ffmpeg") is None:
        raise RuntimeError(
            "ffmpeg is required to create .mov files, but it is not available in this environment."
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with xr.open_dataset(args.input) as ds:
        conn = select_connectivity(ds)
        fill_value = conn.attrs.get("_FillValue", conn.encoding.get("_FillValue"))

        print("Scanning global connectivity range for fixed colorbar...", flush=True)
        vmin, vmax = compute_color_limits(conn, fill_value)
        print(f"Using fixed colorbar limits: vmin={vmin:.6g}, vmax={vmax:.6g}", flush=True)

        with tempfile.TemporaryDirectory(
            prefix="tmp_connectivity_frames_",
            dir=args.output.parent,
        ) as tmp_dir:
            temp_dir = Path(tmp_dir)
            print(f"Creating temporary frames in: {temp_dir}", flush=True)
            render_frames(
                conn=conn,
                fill_value=fill_value,
                vmin=vmin,
                vmax=vmax,
                temp_dir=temp_dir,
                dpi=args.dpi,
                cmap=args.cmap,
            )

            print(f"Encoding movie to: {args.output}", flush=True)
            encode_movie_with_ffmpeg(temp_dir, args.output, args.fps)

    print(f"Movie created successfully: {args.output}", flush=True)


if __name__ == "__main__":
    main()
