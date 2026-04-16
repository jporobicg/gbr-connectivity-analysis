"""Matplotlib figure helpers for the connectivity comparison workflow."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _save_figure(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()


def _time_axis(df: pd.DataFrame, column: str = "time_label") -> tuple[np.ndarray, list[str]]:
    labels = df.sort_values("time_index")[column].tolist()
    return np.arange(len(labels)), labels


def plot_metric_timeseries(comparison_df: pd.DataFrame, path: Path) -> None:
    data = comparison_df[comparison_df["subset"] == "all"].sort_values("time_index")
    if data.empty:
        return

    x, labels = _time_axis(data)
    metrics = [
        ("pearson_r", "Pearson r"),
        ("rmse", "RMSE"),
        ("mae", "MAE"),
        ("frobenius_norm", "Frobenius norm"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    for ax, (column, title) in zip(axes.ravel(), metrics, strict=True):
        ax.plot(x, data[column].to_numpy(), marker="o", linewidth=1.5)
        ax.set_title(title)
        ax.grid(alpha=0.3)
    for ax in axes[-1]:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right")
    _save_figure(path)


def plot_threshold_sensitivity(
    comparison_df: pd.DataFrame,
    ci_df: pd.DataFrame,
    network_df: pd.DataFrame,
    path: Path,
) -> None:
    threshold_metrics = comparison_df[comparison_df["subset"] == "threshold"].copy()
    if threshold_metrics.empty:
        return

    threshold_metrics["threshold_label"] = threshold_metrics["threshold"].map(lambda value: f"{value:g}")
    grouped = threshold_metrics.groupby("threshold_label", sort=False)
    threshold_labels = list(grouped.groups.keys())
    x = np.arange(len(threshold_labels))

    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
    axes[0].plot(x, grouped["pearson_r"].median().to_numpy(), marker="o")
    axes[0].fill_between(x, grouped["pearson_r"].quantile(0.25), grouped["pearson_r"].quantile(0.75), alpha=0.2)
    axes[0].set_title("Threshold Sensitivity")
    axes[0].set_ylabel("Pearson r")
    axes[0].grid(alpha=0.3)

    axes[1].plot(x, grouped["changed_fraction"].median().to_numpy(), marker="o", color="tab:orange")
    axes[1].fill_between(
        x,
        grouped["changed_fraction"].quantile(0.25),
        grouped["changed_fraction"].quantile(0.75),
        alpha=0.2,
        color="tab:orange",
    )
    axes[1].set_ylabel("Changed fraction")
    axes[1].grid(alpha=0.3)

    if not ci_df.empty:
        ci_threshold = ci_df[ci_df["subset"] == "threshold"].copy()
        ci_threshold["threshold_label"] = ci_threshold["threshold"].map(lambda value: f"{value:g}")
        ci_grouped = ci_threshold.groupby("threshold_label", sort=False)
        axes[2].plot(
            x,
            [ci_grouped["fraction_outside_ci"].median().get(label, np.nan) for label in threshold_labels],
            marker="o",
            label="Outside bootstrap CI",
        )
    if not network_df.empty:
        threshold_pivot = network_df.pivot_table(index=["time_index", "threshold"], columns="dataset", values="connectance").reset_index()
        threshold_pivot["threshold_label"] = threshold_pivot["threshold"].map(lambda value: f"{value:g}")
        connectance_gap = threshold_pivot.groupby("threshold_label").apply(
            lambda frame: np.nanmedian(frame.get("bootstrap_mean", np.nan) - frame.get("single", np.nan)),
            include_groups=False,
        )
        axes[2].plot(
            x,
            [connectance_gap.get(label, np.nan) for label in threshold_labels],
            marker="s",
            label="Connectance gap",
        )
    axes[2].set_ylabel("Fraction / gap")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(threshold_labels, rotation=45, ha="right")
    axes[2].grid(alpha=0.3)
    axes[2].legend()
    _save_figure(path)


def plot_uncertainty_timeseries(edge_ci_df: pd.DataFrame, node_ci_df: pd.DataFrame, similarity_df: pd.DataFrame, path: Path) -> None:
    if edge_ci_df.empty and node_ci_df.empty and similarity_df.empty:
        return
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)

    labels: list[str] = []
    if not edge_ci_df.empty:
        edge_all = edge_ci_df[edge_ci_df["subset"] == "all"].sort_values("time_index")
        x, labels = _time_axis(edge_all)
        axes[0, 0].plot(x, edge_all["fraction_outside_ci"], marker="o")
        axes[0, 0].set_title("Edge Coverage")
        axes[0, 0].set_ylabel("Fraction outside CI")
        axes[0, 0].grid(alpha=0.3)

    if not node_ci_df.empty:
        pivot = node_ci_df.pivot(index="time_index", columns="entity_type", values="fraction_outside_ci").sort_index()
        x = np.arange(pivot.shape[0])
        for column in pivot.columns:
            axes[0, 1].plot(x, pivot[column], marker="o", label=column)
        axes[0, 1].set_title("Node Coverage")
        axes[0, 1].set_ylabel("Fraction outside CI")
        axes[0, 1].legend()
        axes[0, 1].grid(alpha=0.3)

    if not similarity_df.empty:
        data = similarity_df.sort_values("time_index")
        x, labels = _time_axis(data)
        axes[1, 0].plot(x, data["mean_pairwise_bootstrap_r"], marker="o", label="Bootstrap vs bootstrap")
        axes[1, 0].plot(x, data["single_vs_bootstrap_mean_r"], marker="s", label="Single vs mean")
        axes[1, 0].set_title("Sample Similarity")
        axes[1, 0].set_ylabel("Correlation")
        axes[1, 0].legend()
        axes[1, 0].grid(alpha=0.3)

    if not edge_ci_df.empty:
        edge_threshold = edge_ci_df[edge_ci_df["subset"] == "threshold"]
        if not edge_threshold.empty:
            grouped = edge_threshold.groupby("threshold")["fraction_outside_ci"].median()
            axes[1, 1].plot(grouped.index, grouped.values, marker="o")
            axes[1, 1].set_xscale("symlog", linthresh=1e-10)
            axes[1, 1].set_title("Outside-CI Fraction by Threshold")
            axes[1, 1].set_xlabel("Threshold")
            axes[1, 1].set_ylabel("Fraction")
            axes[1, 1].grid(alpha=0.3)

    for ax in axes[-1]:
        if labels:
            ax.set_xticks(np.arange(len(labels)))
            ax.set_xticklabels(labels, rotation=45, ha="right")
    _save_figure(path)


def plot_heatmap_panel(matrices: dict[str, np.ndarray], time_label: str, path: Path) -> None:
    order = ["single", "bootstrap_mean", "difference", "bootstrap_sd", "bootstrap_cv"]
    titles = {
        "single": "Single matrix",
        "bootstrap_mean": "Bootstrap mean",
        "difference": "Single - bootstrap mean",
        "bootstrap_sd": "Bootstrap SD",
        "bootstrap_cv": "Bootstrap CV",
    }
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    for ax, key in zip(axes.ravel(), order, strict=False):
        if key not in matrices:
            ax.axis("off")
            continue
        cmap = "coolwarm" if key == "difference" else "viridis"
        image = ax.imshow(matrices[key], aspect="auto", cmap=cmap)
        ax.set_title(titles[key])
        ax.set_xlabel("Sink block")
        ax.set_ylabel("Source block")
        plt.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    axes.ravel()[-1].axis("off")
    fig.suptitle(f"Block-Averaged Connectivity Structure: {time_label}")
    _save_figure(path)


def plot_scatter(single_values: np.ndarray, mean_values: np.ndarray, time_label: str, path: Path) -> None:
    valid = np.isfinite(single_values) & np.isfinite(mean_values)
    x = single_values[valid]
    y = mean_values[valid]
    if x.size == 0:
        return
    epsilon = 1e-12
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].scatter(x, y, s=5, alpha=0.2)
    limit = max(float(np.nanmax(x)), float(np.nanmax(y)))
    axes[0].plot([0, limit], [0, limit], color="black", linewidth=1)
    axes[0].set_title(f"Raw Edge Weights: {time_label}")
    axes[0].set_xlabel("Single")
    axes[0].set_ylabel("Bootstrap mean")
    axes[0].grid(alpha=0.3)

    xlog = np.log10(x + epsilon)
    ylog = np.log10(y + epsilon)
    axes[1].scatter(xlog, ylog, s=5, alpha=0.2)
    min_log = min(float(np.nanmin(xlog)), float(np.nanmin(ylog)))
    max_log = max(float(np.nanmax(xlog)), float(np.nanmax(ylog)))
    axes[1].plot([min_log, max_log], [min_log, max_log], color="black", linewidth=1)
    axes[1].set_title(f"log10 Edge Weights: {time_label}")
    axes[1].set_xlabel("log10(single + 1e-12)")
    axes[1].set_ylabel("log10(mean + 1e-12)")
    axes[1].grid(alpha=0.3)
    _save_figure(path)


def plot_distribution_panel(single_sampled_edges: np.ndarray, mean_sampled_edges: np.ndarray, row_sums_single: np.ndarray, row_sums_mean: np.ndarray, col_sums_single: np.ndarray, col_sums_mean: np.ndarray, time_label: str, path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    bins = 50
    axes[0].hist(single_sampled_edges[np.isfinite(single_sampled_edges)], bins=bins, alpha=0.6, label="Single")
    axes[0].hist(mean_sampled_edges[np.isfinite(mean_sampled_edges)], bins=bins, alpha=0.6, label="Bootstrap mean")
    axes[0].set_title("Sampled edge weights")
    axes[0].legend()
    axes[1].hist(row_sums_single[np.isfinite(row_sums_single)], bins=bins, alpha=0.6, label="Single")
    axes[1].hist(row_sums_mean[np.isfinite(row_sums_mean)], bins=bins, alpha=0.6, label="Bootstrap mean")
    axes[1].set_title("Outgoing strength")
    axes[1].legend()
    axes[2].hist(col_sums_single[np.isfinite(col_sums_single)], bins=bins, alpha=0.6, label="Single")
    axes[2].hist(col_sums_mean[np.isfinite(col_sums_mean)], bins=bins, alpha=0.6, label="Bootstrap mean")
    axes[2].set_title("Incoming strength")
    axes[2].legend()
    fig.suptitle(f"Distribution Comparison: {time_label}")
    _save_figure(path)


def plot_rank_panel(row_sums_single: np.ndarray, row_sums_mean: np.ndarray, col_sums_single: np.ndarray, col_sums_mean: np.ndarray, time_label: str, path: Path, top_n: int = 25) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    source_candidates = np.unique(np.concatenate([np.argsort(row_sums_single)[-top_n:], np.argsort(row_sums_mean)[-top_n:]]))
    source_order = source_candidates[np.argsort(row_sums_mean[source_candidates])[::-1]]
    axes[0].plot(row_sums_single[source_order], marker="o", label="Single")
    axes[0].plot(row_sums_mean[source_order], marker="s", label="Bootstrap mean")
    axes[0].set_title("Top exporters")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    sink_candidates = np.unique(np.concatenate([np.argsort(col_sums_single)[-top_n:], np.argsort(col_sums_mean)[-top_n:]]))
    sink_order = sink_candidates[np.argsort(col_sums_mean[sink_candidates])[::-1]]
    axes[1].plot(col_sums_single[sink_order], marker="o", label="Single")
    axes[1].plot(col_sums_mean[sink_order], marker="s", label="Bootstrap mean")
    axes[1].set_title("Top importers")
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    fig.suptitle(f"Top Source/Sink Strengths: {time_label}")
    _save_figure(path)


def plot_rank_correlation_timeseries(rank_summary_df: pd.DataFrame, path: Path) -> None:
    if rank_summary_df.empty:
        return
    min_top = rank_summary_df["top_k"].min()
    data = rank_summary_df[rank_summary_df["top_k"] == min_top].sort_values("time_index")
    fig, ax = plt.subplots(figsize=(10, 4.5))
    for entity_type in sorted(data["entity_type"].unique()):
        subset = data[data["entity_type"] == entity_type]
        ax.plot(np.arange(subset.shape[0]), subset["single_vs_mean_spearman"], marker="o", label=entity_type)
    ax.set_title("Source and Sink Rank Correlation")
    ax.set_ylabel("Spearman correlation")
    ax.set_xticks(np.arange(data["time_label"].nunique()))
    ax.set_xticklabels(sorted(data["time_label"].unique()), rotation=45, ha="right")
    ax.grid(alpha=0.3)
    ax.legend()
    _save_figure(path)


def plot_rank_change_bars(source_changes_df: pd.DataFrame, sink_changes_df: pd.DataFrame, time_label: str, path: Path, top_n: int = 15) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, df, title in (
        (axes[0], source_changes_df.head(top_n), "Sources"),
        (axes[1], sink_changes_df.head(top_n), "Sinks"),
    ):
        if df.empty:
            ax.axis("off")
            continue
        ax.bar(np.arange(df.shape[0]), df["abs_rank_shift"])
        ax.set_title(title)
        ax.set_xlabel("Node")
        ax.set_ylabel("Absolute rank shift")
        ax.set_xticks(np.arange(df.shape[0]))
        ax.set_xticklabels(df["node_id"].astype(str), rotation=90)
        ax.grid(alpha=0.3)
    fig.suptitle(f"Largest Rank Changes: {time_label}")
    _save_figure(path)


def plot_profile_divergence_hist(profile_df: pd.DataFrame, time_label: str, path: Path) -> None:
    if profile_df.empty or "time_label" not in profile_df:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, entity_type in zip(axes, ["source", "sink"], strict=True):
        subset = profile_df[profile_df["entity_type"] == entity_type]
        ax.hist(subset["js_divergence"], bins=50, alpha=0.7)
        ax.set_title(f"{entity_type.title()} profiles")
        ax.set_xlabel("Jensen-Shannon divergence")
        ax.set_ylabel("Count")
        ax.grid(alpha=0.3)
    fig.suptitle(f"Profile Divergence: {time_label}")
    _save_figure(path)


def plot_profile_top_changes(source_df: pd.DataFrame, sink_df: pd.DataFrame, time_label: str, path: Path, top_n: int = 15) -> None:
    if (source_df.empty or "time_label" not in source_df) and (sink_df.empty or "time_label" not in sink_df):
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, df, title in (
        (axes[0], source_df.head(top_n), "Sources"),
        (axes[1], sink_df.head(top_n), "Sinks"),
    ):
        if df.empty:
            ax.axis("off")
            continue
        ax.bar(np.arange(df.shape[0]), df["js_divergence"])
        ax.set_title(title)
        ax.set_ylabel("JS divergence")
        ax.set_xticks(np.arange(df.shape[0]))
        ax.set_xticklabels(df["node_id"].astype(str), rotation=90)
        ax.grid(alpha=0.3)
    fig.suptitle(f"Top Profile Changes: {time_label}")
    _save_figure(path)


def plot_grouping_similarity(grouping_df: pd.DataFrame, path: Path) -> None:
    if grouping_df.empty:
        return
    data = grouping_df.sort_values("time_index")
    x, labels = _time_axis(data)
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.plot(x, data["single_vs_bootstrap_mean_ari"], marker="o", label="Single vs mean")
    ax.plot(x, data["bootstrap_mean_vs_samples_ari_median"], marker="s", label="Mean vs bootstrap samples")
    ax.set_title("Grouping Similarity Through Time")
    ax.set_ylabel("Adjusted Rand index")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.grid(alpha=0.3)
    ax.legend()
    _save_figure(path)


def plot_group_sizes(grouping_size_df: pd.DataFrame, time_label: str, path: Path) -> None:
    if grouping_size_df.empty or "time_label" not in grouping_size_df:
        return
    subset = grouping_size_df[grouping_size_df["time_label"] == time_label]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for dataset in subset["dataset"].unique():
        data = subset[subset["dataset"] == dataset]
        ax.plot(data["group_order"], data["group_size"], marker="o", label=dataset)
    ax.set_title(f"Grouping Size Summary: {time_label}")
    ax.set_xlabel("Ordered group")
    ax.set_ylabel("Group size")
    ax.grid(alpha=0.3)
    ax.legend()
    _save_figure(path)


def plot_lorenz_curves(lorenz_df: pd.DataFrame, time_label: str, path: Path) -> None:
    if lorenz_df.empty or "time_label" not in lorenz_df:
        return
    subset = lorenz_df[lorenz_df["time_label"] == time_label]
    if subset.empty:
        return
    fig, ax = plt.subplots(figsize=(6, 5))
    for dataset in subset["dataset"].unique():
        data = subset[subset["dataset"] == dataset]
        ax.plot(data["edge_fraction"], data["weight_fraction"], label=dataset)
    ax.plot([0, 1], [0, 1], color="black", linewidth=1, linestyle="--")
    ax.set_title(f"Lorenz Curves: {time_label}")
    ax.set_xlabel("Cumulative fraction of edges")
    ax.set_ylabel("Cumulative fraction of total connectivity")
    ax.grid(alpha=0.3)
    ax.legend()
    _save_figure(path)


def plot_concentration_timeseries(concentration_df: pd.DataFrame, path: Path) -> None:
    if concentration_df.empty:
        return
    data = concentration_df[concentration_df["dataset"].isin(["single", "bootstrap_mean"])].sort_values("time_index")
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for dataset in ["single", "bootstrap_mean"]:
        subset = data[data["dataset"] == dataset]
        x, labels = _time_axis(subset)
        axes[0].plot(x, subset["gini"], marker="o", label=dataset)
        if "share_top_1pct" in subset:
            axes[1].plot(x, subset["share_top_1pct"], marker="o", label=dataset)
    axes[0].set_title("Concentration Through Time")
    axes[0].set_ylabel("Gini")
    axes[0].grid(alpha=0.3)
    axes[0].legend()
    axes[1].set_ylabel("Top 1% share")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=45, ha="right")
    axes[1].grid(alpha=0.3)
    if axes[1].lines:
        axes[1].legend()
    _save_figure(path)


def plot_tail_timeseries(tail_df: pd.DataFrame, path: Path) -> None:
    if tail_df.empty:
        return
    data = tail_df[tail_df["dataset"].isin(["single", "bootstrap_mean"])].sort_values("time_index")
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for dataset in ["single", "bootstrap_mean"]:
        subset = data[data["dataset"] == dataset]
        x, labels = _time_axis(subset)
        axes[0].plot(x, subset["max_edge"], marker="o", label=dataset)
        if "top_100_mean" in subset:
            axes[1].plot(x, subset["top_100_mean"], marker="o", label=dataset)
    axes[0].set_title("Tail Metrics Through Time")
    axes[0].set_ylabel("Max edge")
    axes[0].grid(alpha=0.3)
    axes[0].legend()
    axes[1].set_ylabel("Top 100 mean")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=45, ha="right")
    axes[1].grid(alpha=0.3)
    axes[1].legend()
    _save_figure(path)


def plot_top_edge_strength_distribution(single_weights: np.ndarray, mean_weights: np.ndarray, time_label: str, path: Path) -> None:
    if single_weights.size == 0 or mean_weights.size == 0:
        return
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(np.sort(single_weights)[::-1], label="Single")
    ax.plot(np.sort(mean_weights)[::-1], label="Bootstrap mean")
    ax.set_title(f"Top Edge Strength Distribution: {time_label}")
    ax.set_xlabel("Ordered top edges")
    ax.set_ylabel("Edge weight")
    ax.grid(alpha=0.3)
    ax.legend()
    _save_figure(path)


def plot_node_strength_uncertainty(node_df: pd.DataFrame, time_label: str, path: Path, top_n: int = 10) -> None:
    if node_df.empty or "time_label" not in node_df:
        return
    subset = node_df[node_df["time_label"] == time_label]
    if subset.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, entity_type in zip(axes, ["source", "sink"], strict=True):
        data = subset[subset["entity_type"] == entity_type].sort_values("bootstrap_mean_value", ascending=False).head(top_n)
        x = np.arange(data.shape[0])
        lower = np.clip(data["bootstrap_sample_mean"] - data["bootstrap_q_low"], a_min=0.0, a_max=None)
        upper = np.clip(data["bootstrap_q_high"] - data["bootstrap_sample_mean"], a_min=0.0, a_max=None)
        ax.errorbar(
            x,
            data["bootstrap_sample_mean"],
            yerr=[lower, upper],
            fmt="o",
            label="Bootstrap interval",
        )
        ax.scatter(x, data["single_value"], marker="s", label="Single")
        ax.set_title(entity_type.title())
        ax.set_xticks(x)
        ax.set_xticklabels(data["node_id"].astype(str), rotation=90)
        ax.set_ylabel("Strength")
        ax.grid(alpha=0.3)
        ax.legend()
    fig.suptitle(f"Node Strength Uncertainty: {time_label}")
    _save_figure(path)


def plot_node_cv_hist(node_df: pd.DataFrame, time_label: str, path: Path) -> None:
    if node_df.empty or "time_label" not in node_df:
        return
    subset = node_df[node_df["time_label"] == time_label]
    if subset.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, entity_type in zip(axes, ["source", "sink"], strict=True):
        data = subset[subset["entity_type"] == entity_type]["bootstrap_sample_cv"]
        ax.hist(data[np.isfinite(data)], bins=50, alpha=0.7)
        ax.set_title(entity_type.title())
        ax.set_xlabel("CV")
        ax.set_ylabel("Count")
        ax.grid(alpha=0.3)
    fig.suptitle(f"Node Strength CV: {time_label}")
    _save_figure(path)


def plot_edge_overlap_timeseries(edge_overlap_df: pd.DataFrame, path: Path) -> None:
    if edge_overlap_df.empty:
        return
    fig, ax = plt.subplots(figsize=(10, 4.5))
    for top_k in sorted(edge_overlap_df["top_k"].unique()):
        subset = edge_overlap_df[edge_overlap_df["top_k"] == top_k].sort_values("time_index")
        x, labels = _time_axis(subset)
        ax.plot(x, subset["overlap_fraction"], marker="o", label=f"top {top_k}")
    ax.set_title("Top-k Edge Overlap Through Time")
    ax.set_ylabel("Overlap fraction")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.grid(alpha=0.3)
    ax.legend()
    _save_figure(path)


def plot_edge_churn(edge_churn_df: pd.DataFrame, time_label: str, path: Path) -> None:
    if edge_churn_df.empty or "time_label" not in edge_churn_df:
        return
    subset = edge_churn_df[edge_churn_df["time_label"] == time_label]
    if subset.empty:
        return
    top_k = subset["top_k"].min()
    subset = subset[subset["top_k"] == top_k]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, status in zip(axes, ["entering", "leaving"], strict=True):
        data = subset[subset["status"] == status].head(15)
        ax.bar(np.arange(data.shape[0]), data["bootstrap_mean_value"] if status == "entering" else data["single_value"])
        ax.set_title(status.title())
        ax.set_xticks(np.arange(data.shape[0]))
        ax.set_xticklabels((data["source"].astype(str) + "→" + data["sink"].astype(str)).tolist(), rotation=90)
        ax.grid(alpha=0.3)
    fig.suptitle(f"Top-link Churn: {time_label}")
    _save_figure(path)


def plot_stable_edge_counts(stable_edge_summary_df: pd.DataFrame, path: Path) -> None:
    if stable_edge_summary_df.empty:
        return
    data = stable_edge_summary_df[stable_edge_summary_df["rule_type"] == "top_k_frequency"]
    fig, ax = plt.subplots(figsize=(10, 4.5))
    for rule in sorted(data["frequency_rule"].unique()):
        subset = data[data["frequency_rule"] == rule].sort_values("reference_value")
        ax.plot(subset["reference_value"], subset["stable_edge_count"], marker="o", label=f"{int(rule * 100)}% rule")
    ax.set_title("Stable-edge Counts")
    ax.set_xlabel("Top-k definition")
    ax.set_ylabel("Stable edge count")
    ax.grid(alpha=0.3)
    ax.legend()
    _save_figure(path)


def plot_stable_edge_overlap(stable_edge_summary_df: pd.DataFrame, path: Path) -> None:
    if stable_edge_summary_df.empty:
        return
    data = stable_edge_summary_df[stable_edge_summary_df["rule_type"] == "top_k_frequency"]
    if data.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for ax, column, title in (
        (axes[0], "single_overlap_fraction", "Single overlap"),
        (axes[1], "bootstrap_mean_overlap_fraction", "Bootstrap mean overlap"),
    ):
        for rule in sorted(data["frequency_rule"].unique()):
            subset = data[data["frequency_rule"] == rule].sort_values("reference_value")
            ax.plot(subset["reference_value"], subset[column], marker="o", label=f"{int(rule * 100)}% rule")
        ax.set_title(title)
        ax.set_xlabel("Top-k definition")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("Overlap fraction")
    if axes[1].lines:
        axes[1].legend()
    _save_figure(path)


def plot_entropy_strength(local_node_metrics_df: pd.DataFrame, time_label: str, path: Path) -> None:
    if local_node_metrics_df.empty or "time_label" not in local_node_metrics_df:
        return
    subset = local_node_metrics_df[local_node_metrics_df["time_label"] == time_label]
    if subset.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, entity_type in zip(axes, ["source", "sink"], strict=True):
        data = subset[subset["entity_type"] == entity_type]
        ax.scatter(data["bootstrap_mean_effective_n"], data["bootstrap_mean_entropy"], s=8, alpha=0.4)
        ax.set_title(entity_type.title())
        ax.set_xlabel("Effective number of partners")
        ax.set_ylabel("Entropy")
        ax.grid(alpha=0.3)
    fig.suptitle(f"Entropy vs Strength Structure: {time_label}")
    _save_figure(path)


def plot_dominance_ratio(local_node_metrics_df: pd.DataFrame, time_label: str, path: Path) -> None:
    if local_node_metrics_df.empty or "time_label" not in local_node_metrics_df:
        return
    subset = local_node_metrics_df[local_node_metrics_df["time_label"] == time_label]
    if subset.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, entity_type in zip(axes, ["source", "sink"], strict=True):
        data = subset[subset["entity_type"] == entity_type]
        ax.hist(data["single_dominance"][np.isfinite(data["single_dominance"])], bins=50, alpha=0.6, label="Single")
        ax.hist(
            data["bootstrap_mean_dominance"][np.isfinite(data["bootstrap_mean_dominance"])],
            bins=50,
            alpha=0.6,
            label="Bootstrap mean",
        )
        ax.set_title(entity_type.title())
        ax.set_xlabel("Dominance ratio")
        ax.grid(alpha=0.3)
        ax.legend()
    fig.suptitle(f"Dominance Ratio: {time_label}")
    _save_figure(path)


def plot_effective_links_distribution(local_node_metrics_df: pd.DataFrame, time_label: str, path: Path) -> None:
    if local_node_metrics_df.empty or "time_label" not in local_node_metrics_df:
        return
    subset = local_node_metrics_df[local_node_metrics_df["time_label"] == time_label]
    if subset.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, entity_type in zip(axes, ["source", "sink"], strict=True):
        data = subset[subset["entity_type"] == entity_type]
        ax.hist(data["single_effective_n"][np.isfinite(data["single_effective_n"])], bins=50, alpha=0.6, label="Single")
        ax.hist(
            data["bootstrap_mean_effective_n"][np.isfinite(data["bootstrap_mean_effective_n"])],
            bins=50,
            alpha=0.6,
            label="Bootstrap mean",
        )
        ax.set_title(entity_type.title())
        ax.set_xlabel("Effective number of partners")
        ax.grid(alpha=0.3)
        ax.legend()
    fig.suptitle(f"Effective Partners: {time_label}")
    _save_figure(path)


def plot_backbone_overlap(backbone_df: pd.DataFrame, path: Path) -> None:
    if backbone_df.empty:
        return
    fig, ax = plt.subplots(figsize=(10, 4.5))
    grouped = backbone_df.groupby("time_index")
    x = np.arange(len(grouped))
    labels = [group["time_label"].iloc[0] for _, group in grouped]
    ax.plot(x, grouped["jaccard_index"].max().to_numpy(), marker="o")
    ax.set_title("Backbone Overlap Through Time")
    ax.set_ylabel("Best Jaccard overlap")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.grid(alpha=0.3)
    _save_figure(path)


def plot_coverage_rates(coverage_df: pd.DataFrame, path: Path) -> None:
    if coverage_df.empty:
        return
    summary = coverage_df.groupby(["time_index", "time_label", "reference"], as_index=False)["outside_ci"].mean()
    fig, ax = plt.subplots(figsize=(10, 4.5))
    for reference in summary["reference"].unique():
        subset = summary[summary["reference"] == reference].sort_values("time_index")
        x, labels = _time_axis(subset)
        ax.plot(x, subset["outside_ci"], marker="o", label=reference)
    ax.set_title("Coverage Rates Through Time")
    ax.set_ylabel("Fraction outside bootstrap CI")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.grid(alpha=0.3)
    ax.legend()
    _save_figure(path)


def plot_coverage_examples(node_df: pd.DataFrame, time_label: str, path: Path, top_n: int = 8) -> None:
    if node_df.empty or "time_label" not in node_df:
        return
    subset = node_df[node_df["time_label"] == time_label].sort_values("bootstrap_sample_cv", ascending=False).head(top_n)
    if subset.empty:
        return
    fig, ax = plt.subplots(figsize=(10, 4.5))
    x = np.arange(subset.shape[0])
    lower = np.clip(subset["bootstrap_sample_mean"] - subset["bootstrap_q_low"], a_min=0.0, a_max=None)
    upper = np.clip(subset["bootstrap_q_high"] - subset["bootstrap_sample_mean"], a_min=0.0, a_max=None)
    ax.errorbar(
        x,
        subset["bootstrap_sample_mean"],
        yerr=[lower, upper],
        fmt="o",
        label="Bootstrap interval",
    )
    ax.scatter(x, subset["single_value"], marker="s", label="Single")
    ax.set_title(f"Coverage Examples: {time_label}")
    ax.set_xticks(x)
    ax.set_xticklabels((subset["entity_type"] + ":" + subset["node_id"].astype(str)).tolist(), rotation=90)
    ax.grid(alpha=0.3)
    ax.legend()
    _save_figure(path)


def plot_discrepancy_vs_variability(discrepancy_df: pd.DataFrame, path: Path) -> None:
    if discrepancy_df.empty:
        return
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(discrepancy_df.shape[0])
    ax.errorbar(
        x,
        discrepancy_df["within_bootstrap_median"],
        yerr=[
            discrepancy_df["within_bootstrap_median"] - discrepancy_df["within_bootstrap_q025"],
            discrepancy_df["within_bootstrap_q975"] - discrepancy_df["within_bootstrap_median"],
        ],
        fmt="o",
        label="Within bootstrap",
    )
    ax.scatter(x, discrepancy_df["single_vs_bootstrap_mean"], marker="s", label="Single vs mean")
    ax.set_xticks(x)
    ax.set_xticklabels(discrepancy_df["metric_name"], rotation=45, ha="right")
    ax.set_title("Single vs Bootstrap Discrepancy vs Within-bootstrap Variability")
    ax.grid(alpha=0.3)
    ax.legend()
    _save_figure(path)


def plot_topk_probability(rank_df: pd.DataFrame, time_label: str, path: Path, top_k: int) -> None:
    if rank_df.empty or "time_label" not in rank_df:
        return
    subset = rank_df[rank_df["time_label"] == time_label]
    if subset.empty or f"prob_top_{top_k}" not in subset:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, entity_type in zip(axes, ["source", "sink"], strict=True):
        data = subset[subset["entity_type"] == entity_type].sort_values(f"prob_top_{top_k}", ascending=False).head(15)
        ax.bar(np.arange(data.shape[0]), data[f"prob_top_{top_k}"])
        ax.set_title(entity_type.title())
        ax.set_xticks(np.arange(data.shape[0]))
        ax.set_xticklabels(data["node_id"].astype(str), rotation=90)
        ax.set_ylabel(f"P(top {top_k})")
        ax.grid(alpha=0.3)
    fig.suptitle(f"Top-k Probability: {time_label}")
    _save_figure(path)


def plot_rank_uncertainty_intervals(rank_df: pd.DataFrame, time_label: str, path: Path) -> None:
    if rank_df.empty or "time_label" not in rank_df:
        return
    subset = rank_df[rank_df["time_label"] == time_label]
    if subset.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, entity_type in zip(axes, ["source", "sink"], strict=True):
        data = subset[subset["entity_type"] == entity_type].sort_values("bootstrap_rank_median").head(12)
        x = np.arange(data.shape[0])
        lower = np.clip(data["bootstrap_rank_median"] - data["bootstrap_rank_q025"], a_min=0.0, a_max=None)
        upper = np.clip(data["bootstrap_rank_q975"] - data["bootstrap_rank_median"], a_min=0.0, a_max=None)
        ax.errorbar(
            x,
            data["bootstrap_rank_median"],
            yerr=[lower, upper],
            fmt="o",
            label="Bootstrap rank interval",
        )
        ax.scatter(x, data["single_rank"], marker="s", label="Single rank")
        ax.set_title(entity_type.title())
        ax.set_xticks(x)
        ax.set_xticklabels(data["node_id"].astype(str), rotation=90)
        ax.set_ylabel("Rank")
        ax.invert_yaxis()
        ax.grid(alpha=0.3)
        ax.legend()
    fig.suptitle(f"Rank Uncertainty: {time_label}")
    _save_figure(path)


def plot_ordination(coords_df: pd.DataFrame, meta_df: pd.DataFrame, time_label: str, path: Path) -> None:
    if coords_df.empty:
        return
    fig, ax = plt.subplots(figsize=(6, 5))
    bootstrap = coords_df[coords_df["sample_type"] == "bootstrap"]
    single = coords_df[coords_df["sample_type"] == "single"]
    mean = coords_df[coords_df["sample_type"] == "bootstrap_mean"]
    ax.scatter(bootstrap["pc1"], bootstrap["pc2"], s=30, alpha=0.4, label="Bootstrap samples")
    if not mean.empty:
        ax.scatter(mean["pc1"], mean["pc2"], s=80, marker="D", color="tab:green", label="Bootstrap mean")
    ax.scatter(single["pc1"], single["pc2"], s=90, marker="*", color="tab:red", label="Single")
    if not meta_df.empty:
        ax.set_xlabel(f"PC1 ({meta_df['pc1_explained_ratio'].iloc[0]:.2%})")
        ax.set_ylabel(f"PC2 ({meta_df['pc2_explained_ratio'].iloc[0]:.2%})")
    ax.set_title(f"Bootstrap Cloud vs Single Matrix: {time_label}")
    ax.grid(alpha=0.3)
    ax.legend()
    _save_figure(path)


def plot_community_similarity_timeseries(community_summary_df: pd.DataFrame, path: Path) -> None:
    if community_summary_df.empty:
        return
    data = community_summary_df.sort_values("time_index")
    x, labels = _time_axis(data)
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    axes[0].plot(x, data["bootstrap_mean_vs_consensus_ari"], marker="o", label="Mean vs consensus")
    axes[0].plot(x, data["sample_vs_consensus_ari_median"], marker="s", label="Bootstrap samples vs consensus")
    axes[0].plot(x, data["single_vs_consensus_ari"], marker="^", label="Single vs consensus")
    axes[0].set_title("Community Stability Through Time")
    axes[0].set_ylabel("Adjusted Rand index")
    axes[0].grid(alpha=0.3)
    axes[0].legend()
    axes[1].plot(x, data["selected_n_communities"], marker="o", label="Selected communities")
    axes[1].plot(x, data["mean_consensus_probability"], marker="s", label="Mean consensus probability")
    axes[1].set_ylabel("Count / probability")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=45, ha="right")
    axes[1].grid(alpha=0.3)
    axes[1].legend()
    _save_figure(path)


def plot_community_consensus_map(assignments_df: pd.DataFrame, time_label: str, path: Path) -> None:
    subset = assignments_df[assignments_df["time_label"] == time_label]
    if subset.empty or subset["LAT"].isna().all() or subset["LON"].isna().all():
        return
    fig, ax = plt.subplots(figsize=(7.5, 9))
    scatter = ax.scatter(
        subset["LON"],
        subset["LAT"],
        c=subset["consensus_community"],
        s=14,
        cmap="tab20",
        alpha=0.85,
    )
    ax.set_title(f"Consensus Connectivity Regions: {time_label}")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(alpha=0.2)
    plt.colorbar(scatter, ax=ax, label="Consensus community")
    _save_figure(path)


def plot_community_stability_map(assignments_df: pd.DataFrame, time_label: str, path: Path) -> None:
    subset = assignments_df[assignments_df["time_label"] == time_label]
    if subset.empty or subset["LAT"].isna().all() or subset["LON"].isna().all():
        return
    fig, ax = plt.subplots(figsize=(7.5, 9))
    scatter = ax.scatter(
        subset["LON"],
        subset["LAT"],
        c=subset["community_instability"],
        s=16,
        cmap="viridis",
        alpha=0.9,
    )
    ax.set_title(f"Community Instability: {time_label}")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(alpha=0.2)
    plt.colorbar(scatter, ax=ax, label="1 - consensus probability")
    _save_figure(path)


def plot_coassignment_heatmap(coassignment_block: np.ndarray, time_label: str, path: Path) -> None:
    if coassignment_block.size == 0:
        return
    fig, ax = plt.subplots(figsize=(6.5, 6))
    image = ax.imshow(coassignment_block, aspect="auto", cmap="magma", vmin=0.0, vmax=1.0)
    ax.set_title(f"Bootstrap Co-assignment Probability: {time_label}")
    ax.set_xlabel("Ordered reef groups")
    ax.set_ylabel("Ordered reef groups")
    plt.colorbar(image, ax=ax, label="Co-assignment probability")
    _save_figure(path)


def plot_bridge_map(bridge_df: pd.DataFrame, time_label: str, path: Path, top_n: int = 150) -> None:
    subset = bridge_df[bridge_df["time_label"] == time_label]
    if subset.empty or subset["LAT"].isna().all() or subset["LON"].isna().all():
        return
    top = subset.sort_values("bootstrap_mean_bridge_score", ascending=False).head(top_n)
    fig, ax = plt.subplots(figsize=(7.5, 9))
    size = 20 + 180 * (top["bootstrap_mean_bridge_score"] / top["bootstrap_mean_bridge_score"].max())
    scatter = ax.scatter(
        top["LON"],
        top["LAT"],
        c=top["bootstrap_mean_bridge_score"],
        s=size,
        cmap="plasma",
        alpha=0.85,
    )
    ax.set_title(f"Top Stepping-stone Reefs: {time_label}")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(alpha=0.2)
    plt.colorbar(scatter, ax=ax, label="Bootstrap-mean bridge score")
    _save_figure(path)


def plot_bridge_uncertainty(bridge_df: pd.DataFrame, time_label: str, path: Path, top_n: int = 20) -> None:
    subset = bridge_df[bridge_df["time_label"] == time_label]
    if subset.empty:
        return
    prob_col = next((column for column in subset.columns if column.startswith("prob_top_bridge_")), None)
    if prob_col is None:
        return
    data = subset.sort_values(prob_col, ascending=False).head(top_n)
    x = np.arange(data.shape[0])
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.bar(x, data[prob_col])
    ax.set_title(f"Stepping-stone Top-k Probability: {time_label}")
    ax.set_ylabel("Probability")
    ax.set_xticks(x)
    ax.set_xticklabels(data["node_id"].astype(str), rotation=90)
    ax.grid(alpha=0.3)
    _save_figure(path)


def plot_bridge_strength_scatter(bridge_df: pd.DataFrame, time_label: str, path: Path) -> None:
    subset = bridge_df[bridge_df["time_label"] == time_label]
    if subset.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].scatter(subset["bootstrap_mean_source_strength"], subset["bootstrap_mean_bridge_score"], s=10, alpha=0.3)
    axes[0].set_xlabel("Bootstrap mean source strength")
    axes[0].set_ylabel("Bootstrap mean bridge score")
    axes[0].set_title("Bridge vs source strength")
    axes[0].grid(alpha=0.3)
    axes[1].scatter(subset["bootstrap_mean_sink_strength"], subset["bootstrap_mean_bridge_score"], s=10, alpha=0.3)
    axes[1].set_xlabel("Bootstrap mean sink strength")
    axes[1].set_ylabel("Bootstrap mean bridge score")
    axes[1].set_title("Bridge vs sink strength")
    axes[1].grid(alpha=0.3)
    fig.suptitle(f"Stepping-stone Strength Comparison: {time_label}")
    _save_figure(path)


def plot_bridge_removal_impact(bridge_df: pd.DataFrame, time_label: str, path: Path, top_n: int = 20) -> None:
    subset = bridge_df[bridge_df["time_label"] == time_label]
    if subset.empty:
        return
    data = subset.sort_values("bootstrap_mean_removal_impact", ascending=False).head(top_n)
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.bar(np.arange(data.shape[0]), data["bootstrap_mean_removal_impact"])
    ax.set_title(f"Removal Impact on Inter-anchor Exchange: {time_label}")
    ax.set_ylabel("Approximate impact fraction")
    ax.set_xticks(np.arange(data.shape[0]))
    ax.set_xticklabels(data["node_id"].astype(str), rotation=90)
    ax.grid(alpha=0.3)
    _save_figure(path)


def plot_distance_decay(distance_curve_df: pd.DataFrame, time_label: str, path: Path) -> None:
    subset = distance_curve_df[distance_curve_df["time_label"] == time_label]
    if subset.empty:
        return
    x = 0.5 * (subset["bin_start_km"] + subset["bin_end_km"])
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(x, subset["bootstrap_sample_density_mean"], marker="o", label="Bootstrap samples")
    ax.fill_between(
        x,
        subset["bootstrap_sample_density_q025"],
        subset["bootstrap_sample_density_q975"],
        alpha=0.2,
    )
    ax.plot(x, subset["single_connectivity_density"], marker="s", label="Single")
    ax.plot(x, subset["bootstrap_mean_connectivity_density"], marker="^", label="Bootstrap mean")
    ax.set_title(f"Distance-decay with Bootstrap Envelope: {time_label}")
    ax.set_xlabel("Distance bin midpoint (km)")
    ax.set_ylabel("Connectivity density")
    ax.grid(alpha=0.3)
    ax.legend()
    _save_figure(path)


def plot_distance_uncertainty(distance_curve_df: pd.DataFrame, time_label: str, path: Path) -> None:
    subset = distance_curve_df[distance_curve_df["time_label"] == time_label]
    if subset.empty:
        return
    x = 0.5 * (subset["bin_start_km"] + subset["bin_end_km"])
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(x, subset["bootstrap_sample_density_cv"], marker="o")
    ax.set_title(f"Distance-dependent Uncertainty: {time_label}")
    ax.set_xlabel("Distance bin midpoint (km)")
    ax.set_ylabel("Bootstrap CV of density")
    ax.grid(alpha=0.3)
    _save_figure(path)


def plot_long_distance_panel(distance_class_df: pd.DataFrame, path: Path) -> None:
    if distance_class_df.empty:
        return
    data = distance_class_df[distance_class_df["distance_class"] == "long"].sort_values("time_index")
    if data.empty:
        return
    x, labels = _time_axis(data)
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.plot(x, data["single_fraction_of_total"], marker="o", label="Single")
    ax.plot(x, data["bootstrap_mean_fraction_of_total"], marker="s", label="Bootstrap mean")
    ax.fill_between(
        x,
        data["bootstrap_sample_fraction_q025"],
        data["bootstrap_sample_fraction_q975"],
        alpha=0.2,
        label="Bootstrap interval",
    )
    ax.set_title("Long-distance Connectivity Fraction")
    ax.set_ylabel("Fraction of total connectivity")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.grid(alpha=0.3)
    ax.legend()
    _save_figure(path)


def plot_regional_exchange_heatmap(exchange_df: pd.DataFrame, time_label: str, path: Path, value_col: str = "bootstrap_mean_value") -> None:
    subset = exchange_df[exchange_df["time_label"] == time_label]
    if subset.empty:
        return
    pivot = subset.pivot(index="source_region", columns="sink_region", values=value_col)
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    image = ax.imshow(pivot.values, aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right")
    ax.set_yticks(np.arange(pivot.shape[0]))
    ax.set_yticklabels(pivot.index)
    ax.set_title(f"Regional Exchange Heatmap: {time_label}")
    plt.colorbar(image, ax=ax, label=value_col.replace("_", " "))
    _save_figure(path)


def plot_regional_roles(role_df: pd.DataFrame, time_label: str, path: Path) -> None:
    subset = role_df[role_df["time_label"] == time_label]
    if subset.empty:
        return
    x = np.arange(subset.shape[0])
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.bar(x - 0.2, subset["bootstrap_mean_export"], width=0.4, label="Export")
    ax.bar(x + 0.2, subset["bootstrap_mean_import"], width=0.4, label="Import")
    ax.set_title(f"Regional Export and Import Roles: {time_label}")
    ax.set_xticks(x)
    ax.set_xticklabels(subset["region"], rotation=45, ha="right")
    ax.grid(alpha=0.3)
    ax.legend()
    _save_figure(path)


def plot_regional_pair_intervals(exchange_df: pd.DataFrame, time_label: str, path: Path, top_n: int = 12) -> None:
    subset = exchange_df[exchange_df["time_label"] == time_label].copy()
    if subset.empty:
        return
    subset["pair_label"] = subset["source_region"] + " -> " + subset["sink_region"]
    data = subset.sort_values("bootstrap_mean_value", ascending=False).head(top_n)
    x = np.arange(data.shape[0])
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.errorbar(
        x,
        data["bootstrap_sample_mean"],
        yerr=[data["bootstrap_sample_mean"] - data["bootstrap_q_low"], data["bootstrap_q_high"] - data["bootstrap_sample_mean"]],
        fmt="o",
        label="Bootstrap interval",
    )
    ax.scatter(x, data["single_value"], marker="s", label="Single")
    ax.scatter(x, data["bootstrap_mean_value"], marker="^", label="Bootstrap mean")
    ax.set_title(f"Regional Pair Exchange Uncertainty: {time_label}")
    ax.set_xticks(x)
    ax.set_xticklabels(data["pair_label"], rotation=45, ha="right")
    ax.grid(alpha=0.3)
    ax.legend()
    _save_figure(path)


def plot_mean_bias_summary(mean_bias_df: pd.DataFrame, path: Path) -> None:
    if mean_bias_df.empty:
        return
    mean_bias_df = mean_bias_df.copy()
    mean_bias_df["label"] = mean_bias_df["analysis_domain"] + "\n" + mean_bias_df["metric_name"]
    pivot = mean_bias_df.pivot_table(index="label", columns="reference", values="reference_percentile", aggfunc="median")
    fig, ax = plt.subplots(figsize=(8, max(4.5, 0.4 * pivot.shape[0])))
    image = ax.imshow(pivot.values, aspect="auto", cmap="coolwarm", vmin=0, vmax=100)
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns)
    ax.set_yticks(np.arange(pivot.shape[0]))
    ax.set_yticklabels(pivot.index)
    ax.set_title("Mean and Single Reference Percentiles Within Bootstrap Distributions")
    plt.colorbar(image, ax=ax, label="Percentile")
    _save_figure(path)


def plot_spatial_hotspots(spatial_hotspot_df: pd.DataFrame, time_label: str, path: Path) -> None:
    subset = spatial_hotspot_df[spatial_hotspot_df["time_label"] == time_label]
    if subset.empty or subset["LAT"].isna().all() or subset["LON"].isna().all():
        return
    metrics = list(subset["metric_name"].unique())[:4]
    fig, axes = plt.subplots(2, 2, figsize=(10, 11))
    for ax, metric_name in zip(axes.ravel(), metrics, strict=False):
        data = subset[subset["metric_name"] == metric_name]
        scatter = ax.scatter(data["LON"], data["LAT"], c=data["hotspot_score"], s=14, cmap="coolwarm", alpha=0.85)
        ax.set_title(metric_name.replace("_", " ").title())
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.grid(alpha=0.2)
        plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
    for ax in axes.ravel()[len(metrics):]:
        ax.axis("off")
    fig.suptitle(f"Spatial Hotspots of Disagreement: {time_label}")
    _save_figure(path)


def plot_spatial_clustering(spatial_summary_df: pd.DataFrame, path: Path) -> None:
    if spatial_summary_df.empty:
        return
    data = spatial_summary_df.groupby("metric_name", as_index=False)["morans_i"].median()
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(np.arange(data.shape[0]), data["morans_i"])
    ax.set_xticks(np.arange(data.shape[0]))
    ax.set_xticklabels(data["metric_name"], rotation=45, ha="right")
    ax.set_ylabel("Median Moran's I")
    ax.set_title("Spatial Clustering of Uncertainty")
    ax.grid(alpha=0.3)
    _save_figure(path)


def plot_variance_decomposition(variance_summary_df: pd.DataFrame, path: Path) -> None:
    if variance_summary_df.empty:
        return
    components = [
        ("within_time_bootstrap_var", "Within bootstrap"),
        ("among_spawning_periods_var", "Spawning periods"),
        ("among_years_var", "Years"),
        ("among_families_var", "Families"),
    ]
    fig, ax = plt.subplots(figsize=(10, max(4.5, 0.45 * variance_summary_df.shape[0])))
    y = np.arange(variance_summary_df.shape[0])
    left = np.zeros(variance_summary_df.shape[0], dtype=np.float64)
    for column, label in components:
        values = variance_summary_df[column].fillna(0.0).to_numpy(dtype=np.float64)
        ax.barh(y, values, left=left, label=label)
        left += values
    ax.set_yticks(y)
    ax.set_yticklabels(variance_summary_df["metric_name"])
    ax.set_xlabel("Approximate variance component")
    ax.set_title("Variance Decomposition by Metric")
    ax.grid(alpha=0.3, axis="x")
    ax.legend()
    _save_figure(path)


def plot_family_similarity_heatmap(family_similarity_df: pd.DataFrame, path: Path, value_col: str = "source_strength_spearman") -> None:
    if family_similarity_df.empty:
        return
    families = sorted(set(family_similarity_df["family_a"]).union(family_similarity_df["family_b"]))
    matrix = pd.DataFrame(np.eye(len(families)), index=families, columns=families)
    for row in family_similarity_df.itertuples(index=False):
        matrix.loc[row.family_a, row.family_b] = getattr(row, value_col)
        matrix.loc[row.family_b, row.family_a] = getattr(row, value_col)
    fig, ax = plt.subplots(figsize=(6, 5))
    image = ax.imshow(matrix.values, cmap="viridis", vmin=np.nanmin(matrix.values), vmax=np.nanmax(matrix.values))
    ax.set_xticks(np.arange(len(families)))
    ax.set_xticklabels(families, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(families)))
    ax.set_yticklabels(families)
    ax.set_title(f"Family Similarity Heatmap: {value_col.replace('_', ' ')}")
    plt.colorbar(image, ax=ax)
    _save_figure(path)


def plot_family_overlap(family_overlap_df: pd.DataFrame, path: Path) -> None:
    if family_overlap_df.empty:
        return
    data = family_overlap_df.copy()
    data["pair_label"] = data["family_a"] + " vs " + data["family_b"]
    pair_labels = sorted(data["pair_label"].unique())
    metric_names = list(data["metric_name"].unique())
    fig, ax = plt.subplots(figsize=(10, 4.5))
    for metric_index, metric_name in enumerate(metric_names):
        subset = data[data["metric_name"] == metric_name].set_index("pair_label").reindex(pair_labels).reset_index()
        ax.bar(
            np.arange(len(pair_labels)) + 0.2 * metric_index,
            subset["overlap_fraction"],
            width=0.2,
            label=metric_name,
        )
    ax.set_xticks(np.arange(len(pair_labels)) + 0.2)
    ax.set_xticklabels(pair_labels, rotation=45, ha="right")
    ax.set_ylabel("Top-k overlap")
    ax.set_title("Family Overlap of Important Reefs")
    ax.grid(alpha=0.3)
    ax.legend()
    _save_figure(path)


def plot_family_ordination(family_ordination_df: pd.DataFrame, path: Path) -> None:
    if family_ordination_df.empty:
        return
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(family_ordination_df["pc1"], family_ordination_df["pc2"], s=80)
    for row in family_ordination_df.itertuples(index=False):
        ax.text(row.pc1, row.pc2, row.family)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("Family Ordination")
    ax.grid(alpha=0.3)
    _save_figure(path)


def plot_family_uncertainty(family_uncertainty_df: pd.DataFrame, path: Path) -> None:
    if family_uncertainty_df.empty:
        return
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    metrics = [
        ("median_edge_cv_nonzero", "Median edge CV"),
        ("median_bridge_cv", "Median bridge CV"),
        ("mean_consensus_probability", "Consensus probability"),
    ]
    x = np.arange(family_uncertainty_df.shape[0])
    for ax, (column, title) in zip(axes, metrics, strict=True):
        ax.bar(x, family_uncertainty_df[column].fillna(0.0))
        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xticklabels(family_uncertainty_df["family"], rotation=45, ha="right")
        ax.grid(alpha=0.3)
    _save_figure(path)


def plot_ecological_threshold_robustness(threshold_df: pd.DataFrame, path: Path) -> None:
    if threshold_df.empty:
        return
    metrics = [
        ("community_ari_vs_baseline", "Community ARI"),
        ("bridge_top_overlap_fraction", "Bridge overlap"),
        ("backbone_jaccard_vs_baseline", "Backbone Jaccard"),
    ]
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    labels = threshold_df["threshold_mode"].tolist()
    x = np.arange(len(labels))
    for ax, (column, title) in zip(axes, metrics, strict=True):
        ax.bar(x, threshold_df[column])
        ax.set_ylabel(title)
        ax.grid(alpha=0.3)
    axes[0].set_title("Ecological Threshold Robustness")
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(labels, rotation=45, ha="right")
    _save_figure(path)
