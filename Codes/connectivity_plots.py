"""
Plot connectivity matrices and connection maps from NetCDF files.
"""
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


def _get_connectivity_matrix(
    ds: xr.Dataset,
    time_idx=None,
    treatment_idx=None,
    sample_idx=None,
    aggregate="sum",
):
    """
    Extract a 2D connectivity matrix from an xarray Dataset.

    Handles different NetCDF layouts:
    - (source, sink) or (source, target)
    - with optional time, treatment, sample dimensions (reduced by aggregate or indices).
    """
    if "connectivity" not in ds.data_vars:
        raise ValueError("Dataset has no 'connectivity' variable")

    conn = ds["connectivity"]

    # Reduce extra dimensions
    dims = set(conn.dims)
    # Standard names for the two spatial dimensions
    spatial = []
    for name in ["source", "sink", "target"]:
        if name in dims:
            spatial.append(name)
    if len(spatial) < 2:
        raise ValueError(
            f"Connectivity must have two spatial dimensions (e.g. source/sink or source/target). Found: {list(dims)}"
        )

    # Aggregate or select time
    if "time" in dims:
        if time_idx is not None:
            conn = conn.isel(time=time_idx)
        elif aggregate == "sum":
            conn = conn.sum(dim="time")
        elif aggregate == "mean":
            conn = conn.mean(dim="time")
        else:
            conn = conn.isel(time=0)
    if "treatment" in conn.dims:
        conn = conn.isel(treatment=treatment_idx if treatment_idx is not None else 0)
    if "sample" in conn.dims:
        conn = conn.isel(sample=sample_idx if sample_idx is not None else 0)

    return conn.values


def plot_connectivity_matrix(
    nc_path,
    time_idx=None,
    aggregate="sum",
    ax=None,
    title=None,
    cmap="viridis",
    log_scale=False,
    colorbar=True,
    **kwargs,
):
    """
    Plot the connectivity matrix from a NetCDF file.

    Parameters
    ----------
    nc_path : str or Path
        Path to the NetCDF file (e.g. connectivity matrix).
    time_idx : int, optional
        If given, use this time index instead of aggregating over time.
    aggregate : str, optional
        How to reduce time dimension if time_idx is None: "sum", "mean", or "first".
    ax : matplotlib axes, optional
        Axes to plot on. If None, current axes or new figure is used.
    title : str, optional
        Plot title.
    cmap : str, optional
        Colormap name (default "viridis").
    log_scale : bool, optional
        If True, use log10(1 + value) for display.
    colorbar : bool, optional
        Whether to add a colorbar.
    **kwargs
        Passed to matplotlib imshow (e.g. aspect, interpolation).

    Returns
    -------
    fig, ax, im
        Figure, axes, and image artist (ax may be None if ax was passed).
    """
    nc_path = Path(nc_path)
    if not nc_path.exists():
        raise FileNotFoundError(f"NetCDF file not found: {nc_path}")

    ds = xr.open_dataset(nc_path)
    try:
        matrix = _get_connectivity_matrix(ds, time_idx=time_idx, aggregate=aggregate)
        # Mask NetCDF fill/missing and zeros so color scale reflects positive data only
        fill = getattr(ds["connectivity"].attrs, "_FillValue", None)
        if fill is not None:
            mask = (matrix >= fill) | ~np.isfinite(matrix) | (matrix == 0)
        else:
            mask = (matrix >= 1e30) | ~np.isfinite(matrix) | (matrix == 0)
        matrix = np.ma.masked_where(mask, matrix)
    finally:
        ds.close()

    mask = np.ma.getmaskarray(matrix)
    if log_scale:
        matrix = np.ma.filled(matrix, 0)
        matrix = np.log10(1.0 + np.maximum(matrix, 0))
        matrix = np.ma.array(matrix, mask=mask)

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 8)))
        return_fig = True
    else:
        fig = ax.get_figure()
        return_fig = False

    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    n_src, n_tgt = matrix.shape
    extent = (0, n_tgt, 0, n_src)  # x: 0 to n_tgt, y: 0 to n_src

    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad(color="lightgray", alpha=0.8)

    vmin = kwargs.pop("vmin", None)
    vmax = kwargs.pop("vmax", None)
    if vmin is None or vmax is None:
        valid = np.ma.filled(matrix, np.nan)
        valid = valid[np.isfinite(valid)]
        if len(valid) > 0:
            if vmin is None:
                vmin = float(np.min(valid))
            if vmax is None:
                vmax = float(np.max(valid))
            if vmax <= vmin:
                vmax = vmin + 1.0

    imshow_kw = dict(
        cmap=cmap_obj,
        aspect=kwargs.pop("aspect", "equal"),
        interpolation=kwargs.pop("interpolation", "nearest"),
        origin="lower",
        extent=extent,
        **kwargs,
    )
    if vmin is not None:
        imshow_kw["vmin"] = vmin
    if vmax is not None:
        imshow_kw["vmax"] = vmax

    im = ax.imshow(matrix, **imshow_kw)
    if colorbar:
        plt.colorbar(im, ax=ax, label="Connectivity")
    if title is not None:
        ax.set_title(title)
    else:
        ax.set_title(nc_path.name)

    ax.set_xlabel("Target")
    ax.set_ylabel("Source")

    if return_fig:
        return fig, ax, im
    return fig, ax, im


def plot_connections_per_source_histogram(
    nc_path,
    time_idx=0,
    treatment_idx=None,
    sample_idx=None,
    valid_max=1e16,
    ax=None,
    title=None,
    **kwargs,
):
    """
    Histogram: x = source index, y = count of valid positive connections per source.
    Valid connection = value > 0 and value < valid_max (and not fill, finite).

    Parameters
    ----------
    nc_path : str or Path
        Path to the NetCDF connectivity file.
    time_idx : int, optional
        Time slice to use (default 0).
    treatment_idx : int, optional
        Treatment slice (default 0). Used only if dataset has treatment dimension.
    sample_idx : int, optional
        Sample slice (default 0). Used only if dataset has sample dimension.
    valid_max : float, optional
        Values >= valid_max are treated as invalid (default 1e16).
    ax : matplotlib axes, optional
    title : str, optional
    **kwargs
        Passed to ax.bar() or ax.step().

    Returns
    -------
    fig, ax
    """
    nc_path = Path(nc_path)
    if not nc_path.exists():
        raise FileNotFoundError(f"NetCDF file not found: {nc_path}")

    ds = xr.open_dataset(nc_path)
    try:
        matrix = _get_connectivity_matrix(
            ds,
            time_idx=time_idx,
            treatment_idx=treatment_idx,
            sample_idx=sample_idx,
            aggregate=None,
        )
        fill = getattr(ds["connectivity"].attrs, "_FillValue", 1e30)
    finally:
        ds.close()

    valid = np.isfinite(matrix) & (matrix < fill) & (matrix > 0) & (matrix < valid_max)
    count_per_source = np.sum(valid, axis=1)  # shape (n_sources,)

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (12, 4)))
        return_fig = True
    else:
        fig = ax.get_figure()
        return_fig = False

    sources = np.arange(len(count_per_source))
    ax.bar(sources, count_per_source, width=1.0, **kwargs)
    ax.set_xlabel("Source")
    ax.set_ylabel("Count of valid positive connections")
    if title is not None:
        ax.set_title(title)
    else:
        ax.set_title(f"Connections per source ({nc_path.name})")

    if return_fig:
        return fig, ax
    return fig, ax


def count_sources_with_valid_connections_per_slice(
    nc_path,
    valid_max=1e16,
    time_indices=None,
    treatment_indices=None,
    sample_indices=None,
):
    """
    Count how many sources have at least one valid connection (0 < value < valid_max,
    not fill, finite) for each (time, treatment, sample) slice.

    Parameters
    ----------
    nc_path : str or Path
        Path to the NetCDF connectivity file.
    valid_max : float, optional
        Values >= valid_max are invalid (default 1e16).
    time_indices, treatment_indices, sample_indices : list or None, optional
        If given, only these indices are computed (e.g. time_indices=[0, 1]).
        If None, all indices for that dimension are used.

    Returns
    -------
    counts : ndarray
        Shape (n_time, n_treatment, n_sample). counts[t,tr,s] = number of sources
        with at least one valid connection in that slice. Uncomputed cells are -1.
    """
    nc_path = Path(nc_path)
    if not nc_path.exists():
        raise FileNotFoundError(f"NetCDF file not found: {nc_path}")

    ds = xr.open_dataset(nc_path)
    conn = ds["connectivity"]
    fill = getattr(conn.attrs, "_FillValue", 1e30)
    dims = list(conn.dims)
    n_time = conn.sizes.get("time", 1)
    n_treatment = conn.sizes.get("treatment", 1)
    n_sample = conn.sizes.get("sample", 1)

    times = list(range(n_time)) if time_indices is None else time_indices
    treatments = list(range(n_treatment)) if treatment_indices is None else treatment_indices
    samples = list(range(n_sample)) if sample_indices is None else sample_indices

    counts = np.full((n_time, n_treatment, n_sample), -1, dtype=int)
    for t in times:
        for tr in treatments:
            for s in samples:
                sel = {}
                if "time" in dims:
                    sel["time"] = t
                if "treatment" in dims:
                    sel["treatment"] = tr
                if "sample" in dims:
                    sel["sample"] = s
                mat = conn.isel(sel).values
                if mat.ndim > 2:
                    mat = mat.reshape(mat.shape[0], mat.shape[1], -1)[:, :, 0]
                valid = (
                    np.isfinite(mat)
                    & (mat < fill)
                    & (mat > 0)
                    & (mat < valid_max)
                )
                counts[t, tr, s] = int(np.sum(np.any(valid, axis=1)))
    ds.close()
    return counts


def plot_connectivity_map(
    nc_path,
    shapefile_path,
    time_idx=None,
    aggregate="sum",
    min_value=0,
    max_connections=None,
    min_percentile=None,
    ax=None,
    title=None,
    line_alpha_scale=1.0,
    line_color="black",
    base_facecolor="lightgray",
    base_edgecolor="gray",
    **kwargs,
):
    """
    Plot a map of connectivity: shapefile polygons and lines between centroids
    for source–target pairs with connectivity above a threshold.

    Matrix indices are mapped to shapefile rows by position: index i
    corresponds to the i-th polygon in the shapefile. If the matrix is smaller
    than the number of polygons (e.g. 30x30), only the first N polygons are
    used as nodes.

    For large matrices (e.g. 3806x3806), use max_connections or min_percentile
    to limit the number of lines drawn.

    Parameters
    ----------
    nc_path : str or Path
        Path to the NetCDF connectivity file.
    shapefile_path : str or Path
        Path to the shapefile (e.g. GBR reef polygons).
    time_idx : int, optional
        Time index to use; if None, time is reduced using aggregate.
    aggregate : str, optional
        How to reduce time: "sum", "mean", or "first".
    min_value : float, optional
        Only draw connections with connectivity > min_value.
    max_connections : int, optional
        If set, only draw the top max_connections strongest links (by value).
    min_percentile : float, optional
        If set (e.g. 99), only draw connections above this percentile of non-zero values.
    ax : matplotlib axes, optional
        Axes to plot on (must support GeoDataFrame.plot). If None, new figure.
    title : str, optional
        Plot title.
    line_alpha_scale : float, optional
        Scale for line alpha (alpha = min(1, value * line_alpha_scale)).
    line_color : str, optional
        Color of connection lines.
    base_facecolor, base_edgecolor : str, optional
        Fill and edge color for the shapefile polygons.
    **kwargs
        Passed to GeoDataFrame.plot for the base map.

    Returns
    -------
    fig, ax
        Figure and axes.
    """
    nc_path = Path(nc_path)
    shapefile_path = Path(shapefile_path)
    if not nc_path.exists():
        raise FileNotFoundError(f"NetCDF file not found: {nc_path}")
    if not shapefile_path.exists():
        raise FileNotFoundError(f"Shapefile not found: {shapefile_path}")

    ds = xr.open_dataset(nc_path)
    try:
        matrix = _get_connectivity_matrix(ds, time_idx=time_idx, aggregate=aggregate)
    finally:
        ds.close()

    gdf = gpd.read_file(shapefile_path)
    n_poly = len(gdf)
    n_nodes = matrix.shape[0]

    if n_nodes > n_poly:
        raise ValueError(
            f"Connectivity matrix size ({n_nodes}) is larger than number of "
            f"polygons in shapefile ({n_poly}). Cannot map indices to geometries."
        )

    # Use first n_nodes polygons (by row order) for centroids
    gdf_nodes = gdf.iloc[:n_nodes]
    # Reproject to a suitable CRS for mapping (e.g. Web Mercator for GBR)
    if gdf_nodes.crs and gdf_nodes.crs != "EPSG:3857":
        gdf_nodes = gdf_nodes.to_crs("EPSG:3857")
    centroids = gdf_nodes.geometry.centroid

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 10)))
        return_fig = True
    else:
        fig = ax.get_figure()
        return_fig = False

    # Base map: full shapefile for context (same CRS)
    if gdf.crs != gdf_nodes.crs:
        gdf_plot = gdf.to_crs(gdf_nodes.crs)
    else:
        gdf_plot = gdf
    gdf_plot.plot(ax=ax, facecolor=base_facecolor, edgecolor=base_edgecolor, **kwargs)

    # Collect (i, j, val) for connections above threshold (vectorized)
    i_idx, j_idx = np.where(matrix > min_value)
    vals = matrix[i_idx, j_idx]
    if len(vals) == 0:
        connections = []
    else:
        if min_percentile is not None:
            thresh = np.nanpercentile(vals, min_percentile)
            keep = vals >= thresh
            i_idx, j_idx, vals = i_idx[keep], j_idx[keep], vals[keep]
        if max_connections is not None and len(vals) > max_connections:
            top = np.argsort(-vals)[:max_connections]
            i_idx, j_idx, vals = i_idx[top], j_idx[top], vals[top]
        connections = list(zip(i_idx.tolist(), j_idx.tolist(), vals.tolist()))

    # Draw connections (only the filtered set for performance)
    for i, j, val in connections:
        c_src = centroids.iloc[i]
        c_tgt = centroids.iloc[j]
        alpha = min(1.0, val * line_alpha_scale)
        ax.plot(
            [c_src.x, c_tgt.x],
            [c_src.y, c_tgt.y],
            color=line_color,
            alpha=alpha,
            linewidth=0.5,
        )

    ax.set_axis_off()
    if title is not None:
        ax.set_title(title)
    else:
        ax.set_title(f"Connectivity map: {nc_path.name}")

    if return_fig:
        return fig, ax
    return fig, ax


def main():
    """Generate connectivity matrix plots for the requested NetCDF files."""
    base = Path(__file__).resolve().parent
    data_dir = base.parent / "datasets"
    nc_files = [
        data_dir / "connectivity_matrices" / "connectivity_acroporidae_single.nc",
        data_dir / "connectivity_matrices" / "old_connectivity_matrices.nc",
    ]

    for nc_path in nc_files:
        if not nc_path.exists():
            print(f"Skip (not found): {nc_path}")
            continue
        name = nc_path.stem

        # Plot single slice connectivity(:, :, 1, 1, 1) (0-based: time=0, treatment=0, sample=0)
        fig1, ax1, _ = plot_connectivity_matrix(
            nc_path,
            time_idx=0,
            title=f"Connectivity matrix: {name}",
            log_scale=False,
        )
        fig1.tight_layout()
        out1 = base / f"connectivity_matrix_{name}.png"
        fig1.savefig(out1, dpi=150, bbox_inches="tight")
        plt.close(fig1)
        print(f"Saved: {out1}")

    print("Done.")


if __name__ == "__main__":
    main()
