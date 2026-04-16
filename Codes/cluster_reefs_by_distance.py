#!/usr/bin/env python3
"""
Cluster reef polygons by edge-to-edge distance and export cluster outputs.

Outputs:
1) reef-level shapefile: keeps reef IDs and adds cluster fields
2) cluster polygon shapefile: dissolved cluster geometries
3) cluster centroid shapefile: one centroid point per cluster
4) reef-linked cluster polygon shapefile: keeps reef_id with cluster polygon geometry
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import numpy as np


class UnionFind:
    """Simple disjoint-set structure for connected-component labeling."""

    def __init__(self, size: int) -> None:
        self.parent = np.arange(size, dtype=int)
        self.rank = np.zeros(size, dtype=int)

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            self.parent[ra] = rb
        elif self.rank[ra] > self.rank[rb]:
            self.parent[rb] = ra
        else:
            self.parent[rb] = ra
            self.rank[ra] += 1


def delete_shapefile(path: Path) -> None:
    """Delete existing shapefile sidecar files if they exist."""
    for ext in (".shp", ".shx", ".dbf", ".prj", ".cpg", ".qmd", ".qix"):
        sidecar = path.with_suffix(ext)
        if sidecar.exists():
            sidecar.unlink()


def make_cluster_ids(gdf_metric: gpd.GeoDataFrame, distance_m: float) -> np.ndarray:
    """Return cluster ids using connected components of distance-based adjacency."""
    geoms = gdf_metric.geometry.values
    sindex = gdf_metric.sindex
    n = len(gdf_metric)
    uf = UnionFind(n)

    for i, geom in enumerate(geoms):
        if geom is None or geom.is_empty:
            continue
        minx, miny, maxx, maxy = geom.bounds
        hits = sindex.intersection((minx - distance_m, miny - distance_m, maxx + distance_m, maxy + distance_m))
        for j in hits:
            if j <= i:
                continue
            other = geoms[j]
            if other is None or other.is_empty:
                continue
            if geom.distance(other) <= distance_m:
                uf.union(i, j)

    roots = np.array([uf.find(i) for i in range(n)], dtype=int)
    unique_roots = {root: idx + 1 for idx, root in enumerate(np.unique(roots))}
    return np.array([unique_roots[r] for r in roots], dtype=int)


def split_cluster_indices_by_area(
    geoms_metric: gpd.GeoSeries,
    cent_x: np.ndarray,
    cent_y: np.ndarray,
    indices: np.ndarray,
    max_area_km2: float,
) -> list[np.ndarray]:
    """
    Recursively split one cluster until each part has area <= max_area_km2.

    Splits are done by centroid-based bisection along the largest spatial span.
    """
    finished: list[np.ndarray] = []
    stack: list[np.ndarray] = [indices]

    while stack:
        idx = stack.pop()
        subset = geoms_metric.iloc[idx]
        if hasattr(subset, "union_all"):
            area_m2 = subset.union_all().area
        else:
            area_m2 = subset.unary_union.area
        area_km2 = area_m2 / 1_000_000.0
        if area_km2 <= max_area_km2 or len(idx) <= 1:
            finished.append(idx)
            continue

        xvals = cent_x[idx]
        yvals = cent_y[idx]
        xspan = float(xvals.max() - xvals.min())
        yspan = float(yvals.max() - yvals.min())
        split_axis_vals = xvals if xspan >= yspan else yvals

        order = np.argsort(split_axis_vals, kind="mergesort")
        mid = len(idx) // 2
        left = idx[order[:mid]]
        right = idx[order[mid:]]

        if len(left) == 0 or len(right) == 0:
            finished.append(idx)
            continue

        stack.append(left)
        stack.append(right)

    return finished


def enforce_max_cluster_area(
    gdf_metric: gpd.GeoDataFrame,
    cluster_ids: np.ndarray,
    max_area_km2: float,
) -> np.ndarray:
    """Split oversized clusters and return new cluster ids."""
    geoms = gdf_metric.geometry
    centroids = geoms.centroid
    cent_x = centroids.x.to_numpy()
    cent_y = centroids.y.to_numpy()

    out = np.zeros(len(gdf_metric), dtype=int)
    next_id = 1
    for cid in np.unique(cluster_ids):
        idx = np.flatnonzero(cluster_ids == cid)
        parts = split_cluster_indices_by_area(geoms, cent_x, cent_y, idx, max_area_km2)
        for part in parts:
            out[part] = next_id
            next_id += 1
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Cluster reef polygons by distance.")
    parser.add_argument(
        "--input",
        default="/home/por07g/Documents/Projects/GBR_modeling/Connectivity_analysis/datasets/shapefiles/gbr1_coral_1m_merged_buffer0p001.shp",
        help="Input reef polygon shapefile.",
    )
    parser.add_argument(
        "--output",
        default="/home/por07g/Documents/Projects/GBR_modeling/Connectivity_analysis/datasets/shapefiles/cluster_reefs.shp",
        help="Output reef-level shapefile with cluster fields.",
    )
    parser.add_argument(
        "--distance-m",
        type=float,
        default=0.0,
        help="Maximum edge-to-edge distance (meters) to connect reefs into a cluster.",
    )
    parser.add_argument(
        "--id-field",
        default="FID",
        help="Reef ID field name to preserve (default: FID).",
    )
    parser.add_argument(
        "--metric-crs",
        default="EPSG:3577",
        help="Projected CRS used for distance/area calculations.",
    )
    parser.add_argument(
        "--max-area-km2",
        type=float,
        default=None,
        help="Optional max cluster area (km2). Oversized clusters are recursively split.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    cluster_shape_path = output_path.with_name(f"{output_path.stem}_clusters.shp")
    cluster_centroid_path = output_path.with_name(f"{output_path.stem}_centroids.shp")
    reef_cluster_poly_path = output_path.with_name(f"{output_path.stem}_reef_cluster_poly.shp")

    reefs = gpd.read_file(input_path)
    if reefs.empty:
        raise ValueError(f"No features found in {input_path}")

    if args.id_field in reefs.columns:
        reefs["reef_id"] = reefs[args.id_field]
    else:
        reefs["reef_id"] = np.arange(1, len(reefs) + 1, dtype=int)

    reefs_metric = reefs.to_crs(args.metric_crs)
    initial_clusters = make_cluster_ids(reefs_metric, max(0.0, args.distance_m))
    if args.max_area_km2 is not None and args.max_area_km2 > 0:
        reefs["cluster_id"] = enforce_max_cluster_area(reefs_metric, initial_clusters, args.max_area_km2)
    else:
        reefs["cluster_id"] = initial_clusters

    clusters = reefs.dissolve(by="cluster_id", as_index=False)[["cluster_id", "geometry"]]
    clusters_metric = clusters.to_crs(args.metric_crs)
    # Store area in km2 to avoid DBF width overflow with very large m2 values.
    clusters["cl_ar_km2"] = clusters_metric.geometry.area / 1_000_000.0

    cent_metric = clusters_metric.geometry.centroid
    centroids = gpd.GeoSeries(cent_metric, crs=args.metric_crs).to_crs(reefs.crs)
    clusters["cl_cx"] = centroids.x
    clusters["cl_cy"] = centroids.y

    reefs = reefs.merge(clusters[["cluster_id", "cl_ar_km2", "cl_cx", "cl_cy"]], on="cluster_id", how="left")

    cluster_points = gpd.GeoDataFrame(
        clusters[["cluster_id", "cl_ar_km2", "cl_cx", "cl_cy"]].copy(),
        geometry=centroids,
        crs=reefs.crs,
    )
    cluster_geom_by_id = clusters.set_index("cluster_id").geometry
    reef_cluster_polys = gpd.GeoDataFrame(
        reefs[["reef_id", "cluster_id", "cl_ar_km2", "cl_cx", "cl_cy"]].copy(),
        geometry=reefs["cluster_id"].map(cluster_geom_by_id),
        crs=reefs.crs,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    delete_shapefile(output_path)
    delete_shapefile(cluster_shape_path)
    delete_shapefile(cluster_centroid_path)
    delete_shapefile(reef_cluster_poly_path)

    reefs.to_file(output_path, driver="ESRI Shapefile")
    clusters.to_file(cluster_shape_path, driver="ESRI Shapefile")
    cluster_points.to_file(cluster_centroid_path, driver="ESRI Shapefile")
    reef_cluster_polys.to_file(reef_cluster_poly_path, driver="ESRI Shapefile")

    print(f"Reef clusters written to: {output_path}")
    print(f"Cluster shapes written to: {cluster_shape_path}")
    print(f"Cluster centroids written to: {cluster_centroid_path}")
    print(f"Reef-linked cluster polygons written to: {reef_cluster_poly_path}")
    print(f"Distance threshold (m): {max(0.0, args.distance_m)}")
    if args.max_area_km2 is not None and args.max_area_km2 > 0:
        print(f"Max cluster area limit (km2): {args.max_area_km2}")
    print(f"Total reefs: {len(reefs)}")
    print(f"Total clusters: {clusters['cluster_id'].nunique()}")


if __name__ == "__main__":
    main()
