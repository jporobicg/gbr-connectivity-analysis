"""Project configuration for connectivity dataset comparison."""

from __future__ import annotations

from pathlib import Path
from typing import Any


PROJECT_DIR = Path(__file__).resolve().parent
REPO_DIR = PROJECT_DIR.parent

DEFAULT_BOOTSTRAP_PATH = REPO_DIR / "datasets" / "connectivity_matrices" / "connectivity_acroporidae.nc"
DEFAULT_SINGLE_PATH = REPO_DIR / "datasets" / "connectivity_matrices" / "connectivity_acroporidae_single.nc"
DEFAULT_REEF_METADATA_PATH = REPO_DIR / "datasets" / "reefs" / "Reefs2024.csv"
DEFAULT_FAMILY_DATASET_DIR = REPO_DIR / "datasets" / "connectivity_matrices"

DEFAULT_OUTPUT_DIR = PROJECT_DIR / "outputs"
DEFAULT_TABLES_DIR = DEFAULT_OUTPUT_DIR / "tables"
DEFAULT_FIGURES_DIR = DEFAULT_OUTPUT_DIR / "figures"
DEFAULT_INTERMEDIATE_DIR = DEFAULT_OUTPUT_DIR / "intermediate"


def get_default_config() -> dict[str, Any]:
    """Return the default runtime configuration."""
    return {
        "bootstrap_path": DEFAULT_BOOTSTRAP_PATH,
        "single_path": DEFAULT_SINGLE_PATH,
        "output_dir": DEFAULT_OUTPUT_DIR,
        "tables_dir": DEFAULT_TABLES_DIR,
        "figures_dir": DEFAULT_FIGURES_DIR,
        "intermediate_dir": DEFAULT_INTERMEDIATE_DIR,
        "connectivity_var": "connectivity",
        "treatment_index": 0,
        "reef_metadata_path": DEFAULT_REEF_METADATA_PATH,
        "family_dataset_dir": DEFAULT_FAMILY_DATASET_DIR,
        "community_anchor_column": "AIMS_sector",
        "region_column": "AREA_DESCR",
        "distance_var_candidates": ["direction", "distance"],
        "block_size": 512,
        "thresholds": [0.0, 1e-10, 1e-8, 1e-6, 1e-4],
        "component_thresholds": [1e-8, 1e-6, 1e-4],
        "top_k_edges": [100, 500, 1000],
        "rank_top_k_values": [10, 50, 100],
        "top_k_nodes": 20,
        "bridge_top_k": 50,
        "family_top_k": 50,
        "bias_top_k_nodes": 50,
        "bias_top_k_edges": 500,
        "bias_top_k_region_pairs": 5,
        "bias_stable_frequency": 0.8,
        "top_n_rank_changes": 50,
        "top_n_links": 100,
        "top_share_fractions": [0.01, 0.05, 0.10],
        "stable_edge_thresholds": [0.0, 1e-8, 1e-6],
        "stable_frequency_rules": [0.5, 0.8, 0.95],
        "local_share_thresholds": [0.05, 0.1],
        "grouping_n_clusters": 5,
        "community_min_clusters": 3,
        "community_max_clusters": 8,
        "distance_curve_edges_km": [0, 25, 50, 100, 200, 400, 800, 1600, 2500],
        "distance_class_edges_km": [0, 100, 300, 2500],
        "distance_class_labels": ["short", "medium", "long"],
        "spatial_k_neighbors": 12,
        "compute_quantiles": True,
        "quantiles": [0.025, 0.5, 0.975],
        "cv_epsilon": 1e-12,
        "relative_diff_epsilon": 1e-12,
        "default_change_threshold": 1e-10,
        "scatter_sample_size": 20000,
        "pairwise_sample_size": 50000,
        "spearman_max_points": 200000,
        "heatmap_size": 120,
        "max_component_edges": 2_000_000,
        "ordination_enabled": True,
        "n_selected_times_for_figures": 3,
        "n_selected_times_for_threshold_robustness": 3,
        "include_diagonal": True,
        "random_seed": 42,
        "time_indices": None,
        "skip_plots": False,
        "save_intermediate_npz": True,
    }


def _parse_int_list(value: str | None) -> list[int] | None:
    """Parse a comma-separated integer string."""
    if value is None:
        return None
    items = [item.strip() for item in value.split(",") if item.strip()]
    if not items:
        return None
    return [int(item) for item in items]


def build_runtime_config(args: Any) -> dict[str, Any]:
    """Merge CLI overrides into the default configuration."""
    config = get_default_config()
    if getattr(args, "bootstrap_path", None):
        config["bootstrap_path"] = Path(args.bootstrap_path).expanduser().resolve()
    if getattr(args, "single_path", None):
        config["single_path"] = Path(args.single_path).expanduser().resolve()
    if getattr(args, "output_dir", None):
        config["output_dir"] = Path(args.output_dir).expanduser().resolve()
        config["tables_dir"] = config["output_dir"] / "tables"
        config["figures_dir"] = config["output_dir"] / "figures"
        config["intermediate_dir"] = config["output_dir"] / "intermediate"
    if getattr(args, "block_size", None):
        config["block_size"] = int(args.block_size)
    if getattr(args, "time_indices", None):
        config["time_indices"] = _parse_int_list(args.time_indices)
    if getattr(args, "skip_quantiles", False):
        config["compute_quantiles"] = False
    if getattr(args, "skip_plots", False):
        config["skip_plots"] = True
    return config
