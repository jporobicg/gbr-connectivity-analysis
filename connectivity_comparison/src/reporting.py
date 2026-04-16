"""Markdown reporting helpers for the connectivity comparison workflow."""

from __future__ import annotations

from datetime import datetime

import numpy as np
import pandas as pd


def _safe_float(value: float | int | np.floating | np.integer | None) -> str:
    if value is None or pd.isna(value):
        return "NA"
    return f"{float(value):.4g}"


def _median_or_nan(df: pd.DataFrame, column: str) -> float:
    if df.empty or column not in df:
        return np.nan
    return float(df[column].median())


def build_report(
    comparison_df: pd.DataFrame,
    overlap_df: pd.DataFrame,
    uncertainty_df: pd.DataFrame,
    node_ci_df: pd.DataFrame,
    similarity_df: pd.DataFrame,
    concentration_df: pd.DataFrame,
    profile_df: pd.DataFrame,
    grouping_df: pd.DataFrame,
    coverage_df: pd.DataFrame,
    discrepancy_df: pd.DataFrame,
    rank_uncertainty_df: pd.DataFrame,
    community_df: pd.DataFrame,
    bridge_df: pd.DataFrame,
    distance_df: pd.DataFrame,
    regional_df: pd.DataFrame,
    mean_bias_df: pd.DataFrame,
    family_df: pd.DataFrame,
    variance_df: pd.DataFrame,
    threshold_df: pd.DataFrame,
    config: dict,
) -> str:
    """Create a Markdown report covering the original and new analyses."""
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")

    comparison_all = comparison_df[comparison_df["subset"] == "all"].sort_values("time_index")
    comparison_nonzero = comparison_df[comparison_df["subset"] == "nonzero"].sort_values("time_index")
    uncertainty_all = uncertainty_df[uncertainty_df["subset"] == "all"].sort_values("time_index") if not uncertainty_df.empty else pd.DataFrame()
    outgoing_ci = node_ci_df[node_ci_df["entity_type"].isin(["outgoing", "source"])] if not node_ci_df.empty else pd.DataFrame()
    incoming_ci = node_ci_df[node_ci_df["entity_type"].isin(["incoming", "sink"])] if not node_ci_df.empty else pd.DataFrame()
    source_rank = rank_uncertainty_df[rank_uncertainty_df["entity_type"] == "source"] if not rank_uncertainty_df.empty else pd.DataFrame()
    sink_rank = rank_uncertainty_df[rank_uncertainty_df["entity_type"] == "sink"] if not rank_uncertainty_df.empty else pd.DataFrame()
    source_profile = profile_df[profile_df["entity_type"] == "source"] if not profile_df.empty else pd.DataFrame()
    sink_profile = profile_df[profile_df["entity_type"] == "sink"] if not profile_df.empty else pd.DataFrame()
    overlap_top_min = overlap_df[overlap_df["top_k"] == overlap_df["top_k"].min()] if not overlap_df.empty else pd.DataFrame()
    single_concentration = concentration_df[concentration_df["dataset"] == "single"] if not concentration_df.empty else pd.DataFrame()
    mean_concentration = concentration_df[concentration_df["dataset"] == "bootstrap_mean"] if not concentration_df.empty else pd.DataFrame()
    long_distance = distance_df[distance_df["distance_class"] == "long"] if not distance_df.empty else pd.DataFrame()
    mean_bias_mean = mean_bias_df[mean_bias_df["reference"] == "bootstrap_mean"] if not mean_bias_df.empty else pd.DataFrame()

    largest_diff_row = comparison_all.loc[comparison_all["frobenius_norm"].idxmax()] if not comparison_all.empty else None
    strongest_profile_row = source_profile.loc[source_profile["mean_js_divergence"].idxmax()] if not source_profile.empty else None

    lines: list[str] = [
        "# Connectivity Comparison Report",
        "",
        f"Generated: {generated}",
        "",
        "## Main Answers",
        "",
        "### 1. Are the single and bootstrap-average matrices broadly similar?",
        "",
    ]

    if not comparison_all.empty:
        lines.extend(
            [
                f"- Across all edges, the median Pearson correlation is {_safe_float(comparison_all['pearson_r'].median())}.",
                f"- Across non-zero edges only, the median Pearson correlation is {_safe_float(comparison_nonzero['pearson_r'].median())}.",
                f"- The median outgoing-strength Spearman correlation is {_safe_float(comparison_all['row_sum_spearman_r'].median())}.",
                f"- The median incoming-strength Spearman correlation is {_safe_float(comparison_all['col_sum_spearman_r'].median())}.",
                f"- The median top-link overlap for the smallest top-k set is {_safe_float(_median_or_nan(overlap_top_min, 'overlap_fraction'))}.",
            ]
        )
    else:
        lines.append("- No matrix-comparison results were available.")
    lines.extend(
        [
            "",
            "This analysis adds a whole-matrix view. It tells you whether the single estimate and the bootstrap mean are broadly describing the same system before focusing on specific reefs or pathways.",
            "It is useful as a starting point, but it should be read together with rank, pathway, and uncertainty analyses because high global similarity can still hide meaningful ecological changes.",
            "",
            "### 2. Where do they differ most?",
            "",
        ]
    )

    if largest_diff_row is not None:
        lines.extend(
            [
                f"- The largest structural difference occurs at `{largest_diff_row['time_label']}`.",
                f"- At that time, the Frobenius norm is {_safe_float(largest_diff_row['frobenius_norm'])}, RMSE is {_safe_float(largest_diff_row['rmse'])}, and MAE is {_safe_float(largest_diff_row['mae'])}.",
                f"- The median source-profile Jensen-Shannon divergence is {_safe_float(_median_or_nan(source_profile, 'median_js_divergence'))}.",
                f"- The median sink-profile Jensen-Shannon divergence is {_safe_float(_median_or_nan(sink_profile, 'median_js_divergence'))}.",
            ]
        )
    else:
        lines.append("- No difference summary was available.")
    lines.extend(
        [
            "",
            "This adds a structural-location view. It shows whether differences are only in raw magnitudes or also in where reefs send and receive connectivity.",
            "That matters because two products can have similar total strength yet imply different pathways, destinations, or source regions.",
            "",
            "### 3. Do bootstraps change the identity of key sources, sinks, or pathways?",
            "",
            f"- The median fraction of outgoing strengths outside bootstrap intervals is {_safe_float(_median_or_nan(outgoing_ci, 'fraction_outside_ci'))}.",
            f"- The median fraction of incoming strengths outside bootstrap intervals is {_safe_float(_median_or_nan(incoming_ci, 'fraction_outside_ci'))}.",
            f"- The median fraction of source ranks outside bootstrap rank intervals is {_safe_float(_median_or_nan(source_rank, 'single_rank_outside_ci_fraction'))}.",
            f"- The median fraction of sink ranks outside bootstrap rank intervals is {_safe_float(_median_or_nan(sink_rank, 'single_rank_outside_ci_fraction'))}.",
            f"- The median overlap of strongest pathways remains {_safe_float(_median_or_nan(overlap_top_min, 'overlap_fraction'))}.",
            "",
            "This adds a management-relevance view. It tests whether the reefs and pathways you would call important stay important once bootstrap uncertainty is considered.",
            "If ranks or top pathways shift materially, then ecological conclusions based on the single estimate may not be robust even when broad matrix agreement looks high.",
            "",
            "### 4. Are differences concentrated in weak links or strong links?",
            "",
            f"- The median all-edge Pearson correlation is {_safe_float(_median_or_nan(comparison_all, 'pearson_r'))}, but the non-zero-edge Pearson correlation is {_safe_float(_median_or_nan(comparison_nonzero, 'pearson_r'))}.",
            f"- The median single Gini coefficient is {_safe_float(_median_or_nan(single_concentration, 'gini'))}, while the bootstrap-mean median is {_safe_float(_median_or_nan(mean_concentration, 'gini'))}.",
            f"- The median share of total connectivity held by the top 1% of edges is {_safe_float(_median_or_nan(single_concentration, 'share_top_1pct'))} for the single matrix and {_safe_float(_median_or_nan(mean_concentration, 'share_top_1pct'))} for the bootstrap mean.",
            "",
            "This adds a strength-concentration view. It shows whether the main differences are mostly in weak background links or in the dominant pathways that carry much of the total connectivity.",
            "That distinction is important because disagreement in weak links may matter less ecologically than disagreement in the strongest exporters, sinks, or corridors.",
            "",
            "### 5. Does bootstraping materially improve representation of uncertainty?",
            "",
            f"- The median fraction of single edges outside bootstrap intervals is {_safe_float(_median_or_nan(uncertainty_all, 'fraction_outside_ci'))}.",
            f"- The median sampled-edge single-vs-bootstrap-mean correlation is {_safe_float(_median_or_nan(similarity_df, 'single_vs_bootstrap_mean_r'))}.",
            f"- The median within-bootstrap sampled-edge correlation is {_safe_float(_median_or_nan(similarity_df, 'mean_pairwise_bootstrap_r'))}.",
            f"- The median system-level coverage failure rate is {_safe_float(np.mean(coverage_df['outside_ci'])) if not coverage_df.empty else np.nan}.",
            f"- The median discrepancy-vs-variability comparison for sampled-edge Pearson is {_safe_float(_median_or_nan(discrepancy_df[discrepancy_df['metric_name'] == 'sampled_edge_pearson'], 'single_vs_bootstrap_mean'))} against a within-bootstrap median of {_safe_float(_median_or_nan(discrepancy_df[discrepancy_df['metric_name'] == 'sampled_edge_pearson'], 'within_bootstrap_median'))}.",
            "",
            "This is the central uncertainty result. It shows whether the single estimate behaves like a typical member of the bootstrap cloud or whether the bootstrap ensemble is adding materially new information.",
            "If the single estimate often falls outside bootstrap intervals or outside normal bootstrap-to-bootstrap variability, then the bootstrap approach is not just smoothing noise, it is providing a better representation of uncertainty.",
            "",
            "## Higher-level Ecological Interpretation",
            "",
            f"- The median bootstrap-mean vs consensus community similarity is {_safe_float(_median_or_nan(community_df, 'bootstrap_mean_vs_consensus_ari'))}, with a median consensus probability of {_safe_float(_median_or_nan(community_df, 'mean_consensus_probability'))}.",
            f"- The median stepping-stone CV is {_safe_float(_median_or_nan(bridge_df, 'median_bridge_cv'))}, showing how stable bridge-reef inference is across bootstrap samples.",
            f"- The median bootstrap-mean long-distance fraction is {_safe_float(_median_or_nan(long_distance, 'bootstrap_mean_fraction_of_total'))}, compared with {_safe_float(_median_or_nan(long_distance, 'single_fraction_of_total'))} for the single estimate.",
            f"- The median bootstrap-mean interregional exchange fraction is {_safe_float(_median_or_nan(regional_df, 'bootstrap_mean_interregional_fraction'))}.",
            f"- Across mean-matrix bias diagnostics, the bootstrap mean falls outside the 95% bootstrap interval in {_safe_float(np.mean(mean_bias_mean['outside_95_interval'])) if not mean_bias_mean.empty else np.nan} of interpretation-focused metrics.",
            "",
            "These additions move the workflow from matrix disagreement toward ecological interpretation.",
            "They test whether uncertainty changes large-scale connectivity regions, bridge reefs, long-distance pathways, and region-to-region roles rather than only changing raw edge values.",
            "",
            "## Cross-family And Variance Structure",
            "",
            f"- The workflow currently detects {int(family_df['family'].nunique()) if not family_df.empty else 0} family products for comparison.",
            f"- The dominant variance source is bootstrap uncertainty for {int((variance_df['dominant_source'] == 'within_time_bootstrap').sum()) if not variance_df.empty and 'dominant_source' in variance_df else 0} tracked metrics.",
            f"- The dominant variance source is year-to-year or spawning-period change for {int(variance_df['dominant_source'].isin(['among_years', 'among_spawning_periods']).sum()) if not variance_df.empty and 'dominant_source' in variance_df else 0} tracked metrics.",
            "",
            "This helps separate sampling uncertainty from biologically meaningful temporal change.",
            "It also makes explicit where family-level biology is already producing different connectivity systems and where more bootstrap family products are still needed.",
            "",
            "## Threshold Robustness",
            "",
            f"- The median community similarity across ecological threshold modes is {_safe_float(_median_or_nan(threshold_df, 'community_ari_vs_baseline'))}.",
            f"- The median bridge-top overlap across threshold modes is {_safe_float(_median_or_nan(threshold_df, 'bridge_top_overlap_fraction'))}.",
            f"- The median backbone Jaccard overlap across threshold modes is {_safe_float(_median_or_nan(threshold_df, 'backbone_jaccard_vs_baseline'))}.",
            "",
            "This extends numerical threshold sensitivity into ecological threshold sensitivity.",
            "It tests whether community structure, bridge reefs, and backbone conclusions survive when the definition of an ecologically meaningful link changes.",
            "",
            "## What Each New Analysis Adds",
            "",
            "- Rank uncertainty adds a direct test of whether top sources and sinks are stable or only appear stable in one estimate. This is often closer to ecological interpretation than comparing every edge equally.",
            "- Profile divergence adds a destination/origin view and can reveal changes in where reefs connect even when total strength stays similar. It helps detect rewiring rather than just strengthening or weakening.",
            "- Grouping similarity adds a broad system-organization view and tests whether the same large-scale structure appears under both approaches. It is useful when the question is about system pattern rather than individual links.",
            "- Community structure uncertainty adds a stronger connectivity-region analysis, consensus assignments, and co-assignment probabilities so broad structure can be interpreted with explicit bootstrap uncertainty.",
            "- Stepping-stone analysis adds connector scores and removal-impact proxies so reefs that hold different regions together can be separated from reefs that are simply strong sources or sinks.",
            "- Concentration and tail analyses add a way to test whether the single estimate exaggerates a few strong pathways. They are especially useful for checking whether bootstrap averaging smooths extreme routes.",
            "- Stable-edge and backbone summaries add a robustness view for dominant pathways across bootstrap samples. They help distinguish strong-looking pathways from truly stable ones.",
            "- Distance-dependent uncertainty adds a direct test of whether averaging changes short, medium, and long-distance dispersal interpretation, especially in the long-distance tail.",
            "- Regional exchange adds a management-scale view by translating reef-level links into region-to-region export and import roles with uncertainty intervals.",
            "- Mean-bias diagnostics add a compact answer to the main scientific question: does the bootstrap mean look more stable or more typical than real bootstrap realizations?",
            "- Spatial clustering adds a check of whether disagreement is geographically concentrated, which helps identify hotspots of unstable interpretation rather than treating all reefs as equally uncertain.",
            "- Family comparison and variance decomposition add the wider biological context needed to judge whether bootstrap uncertainty is small compared with family-specific or temporal change.",
            "- Coverage and discrepancy analyses add a direct uncertainty test for whether the single estimate is consistent with the bootstrap ensemble. These are the clearest summaries of what extra information the bootstrap approach provides.",
            "",
            "## Final Synthesis",
            "",
            "The key question is not exact matrix equality. The key question is whether the same ecological story emerges from the single estimate and the bootstrap ensemble.",
            "If ranks, dominant pathways, grouping, and concentration metrics stay stable, then the single matrix may be adequate for broad interpretation even if many weak edges differ.",
            "If the single matrix falls outside bootstrap intervals, changes top sources or sinks, or looks more concentrated than the bootstrap mean, then bootstraping is materially improving uncertainty representation and should affect interpretation.",
            "",
            "## Possible Ecological Interpretation",
            "",
            "- High agreement in source/sink ranks suggests that major exporters and importers are robust to sampling uncertainty. In that case, management conclusions about key reefs are more likely to be stable.",
            "- Large profile divergence with stable total strength suggests reefs may still be important, but they may connect to different places under bootstrap averaging. That can matter if the ecological question depends on who connects to whom, not just how much.",
            "- Lower Gini or weaker tail metrics in the bootstrap mean suggests bootstrap averaging is smoothing sharp pathways and spreading connectivity more broadly. This would imply the single estimate may overemphasize a small number of routes.",
            "- Stable-edge frequency highlights pathways that remain important across bootstrap samples and are therefore stronger candidates for robust ecological conclusions. These links are often more defensible than pathways that are strong in only one estimate.",
            "",
            "## Limitations And Cautions",
            "",
            "- Sparse, zero-heavy matrices can inflate some whole-matrix similarity statistics. This is why non-zero, thresholded, and rank-based summaries are included alongside all-edge metrics.",
            "- Bootstrap sample concentration metrics use tractable summaries and sampled-edge approximations for some tail properties. They are intended to be interpretable and efficient rather than exhaustive.",
            "- Grouping is intentionally simple and should be interpreted as broad structure, not as definitive ecological communities. It is a stability check, not a final classification of reef network organization.",
            "- Threshold and top-k choices matter, so threshold sensitivity should be checked before drawing strong conclusions. If conclusions shift under reasonable settings, uncertainty should be emphasized explicitly.",
        ]
    )
    return "\n".join(lines) + "\n"


def build_analysis_catalog(selected_time_indices: list[int]) -> str:
    """Create a short catalog of analyses, figures, and interpretation notes."""
    time_tag = f"{selected_time_indices[0]:02d}" if selected_time_indices else "XX"
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    entries = [
        (
            "Matrix comparison through time",
            "Tracks whole-matrix agreement between the single estimate and the bootstrap mean across dates.\nIt combines correlation and error-size summaries so you can see both similarity and mismatch.",
            "time_series_matrix_metrics.png",
            "Higher correlation means the overall structure is more similar.\nLower RMSE, MAE, and Frobenius norm mean the two products differ less in absolute terms.",
        ),
        (
            "Threshold sensitivity",
            "Repeats key comparisons after redefining what counts as an active edge.\nThis checks whether the main conclusions depend on weak links or remain stable under stricter cutoffs.",
            "threshold_sensitivity.png",
            "If the curves stay similar across thresholds, the interpretation is robust.\nIf they change sharply, conclusions depend strongly on how small links are treated.",
        ),
        (
            "Edge uncertainty through time",
            "Summarizes how often the single estimate falls outside bootstrap uncertainty intervals.\nIt also compares node-level coverage and sample-to-sample similarity through time.",
            "uncertainty_time_series.png",
            "Higher outside-CI fractions mean the single estimate is less consistent with bootstrap uncertainty.\nIf single-vs-mean similarity is close to within-bootstrap similarity, the single matrix behaves more like a typical bootstrap realization.",
        ),
        (
            "Block heatmap panel",
            "Shows large-scale spatial structure in the single matrix, bootstrap mean, their difference, and bootstrap uncertainty.\nThe matrices are block-averaged so broad patterns can be seen without plotting all edges directly.",
            f"heatmap_panel_time_{time_tag}.png",
            "Large coherent patches mean differences are spatially structured, not just random noise.\nIf difference, SD, or CV hotspots line up with important regions, those regions deserve closer ecological interpretation.",
        ),
        (
            "Edge scatter",
            "Compares sampled edge weights from the single matrix against the bootstrap mean.\nThe raw and log views help reveal both strong-link agreement and behavior among weaker edges.",
            f"scatter_time_{time_tag}.png",
            "Points close to the diagonal mean good agreement.\nSystematic spread away from the diagonal means one product is consistently stronger, weaker, or more variable than the other.",
        ),
        (
            "Edge and node distributions",
            "Compares sampled edge weights plus outgoing and incoming node strengths.\nThis shows whether the single estimate and bootstrap mean differ in overall spread, skew, or central tendency.",
            f"distributions_time_{time_tag}.png",
            "Similar histograms suggest broadly similar strength structure.\nClear shifts or tail differences suggest the bootstrap mean is smoothing or redistributing connectivity.",
        ),
        (
            "Top sources and sinks",
            "Compares the strongest exporters and importers under both products.\nIt focuses on the reefs most likely to matter for ecological interpretation and management.",
            f"rank_panel_time_{time_tag}.png",
            "If the lines track closely, the most important reefs are stable.\nLarge mismatches suggest bootstraping changes which reefs appear most influential.",
        ),
        (
            "Rank correlation through time",
            "Tracks source and sink rank stability through time using the bootstrap distribution.\nIt summarizes whether important reefs stay important across methods and dates.",
            "rank_correlation_time_series.png",
            "Higher values mean the identity of important reefs is more stable.\nDrops through time suggest that bootstraping changes ecological priorities more strongly at some dates than others.",
        ),
        (
            "Rank-change bars",
            "Highlights reefs with the biggest changes in source or sink rank.\nThis makes the main movers easy to inspect without scanning the full rank table.",
            f"rank_change_bars_time_{time_tag}.png",
            "Large bars indicate reefs whose apparent importance changes the most.\nThose reefs are good candidates for follow-up ecological interpretation or mapping.",
        ),
        (
            "Profile divergence histogram",
            "Measures how much each reef’s normalized destination or origin profile changes.\nThis asks whether reefs connect to different places even when total strength stays similar.",
            f"profile_divergence_hist_time_{time_tag}.png",
            "Values near zero mean the connectivity profile is similar between products.\nA long right tail means a subset of reefs changes where it connects even if broad totals are stable.",
        ),
        (
            "Profile divergence top changes",
            "Shows the reefs with the strongest destination-profile or origin-profile change.\nIt focuses attention on the nodes where ecological routing changes most.",
            f"profile_top_changes_time_{time_tag}.png",
            "Higher bars mean stronger rewiring of connectivity patterns.\nThese reefs may keep similar total strength while shifting their main partners.",
        ),
        (
            "Grouping similarity",
            "Compares simple higher-level node groupings between the single matrix, bootstrap mean, and bootstrap samples.\nThis gives a broad system-organization view without relying on heavy community-detection methods.",
            "grouping_similarity_time_series.png",
            "Higher similarity means the large-scale system structure is more stable.\nLower similarity suggests bootstraping changes how the network separates into major functional groupings.",
        ),
        (
            "Lorenz curves",
            "Shows how concentrated total connectivity is among edges.\nIt is the same logic used in inequality analysis: do a few links carry most of the total connectivity?",
            f"lorenz_curves_time_{time_tag}.png",
            "Curves farther from the diagonal mean stronger concentration in fewer links.\nIf the single curve is more bowed than the bootstrap mean, the single estimate is more dominated by a few pathways.",
        ),
        (
            "Concentration through time",
            "Tracks Gini and top-share metrics through time.\nThis shows whether connectivity is consistently diffuse or concentrated into a small set of strong pathways.",
            "concentration_time_series.png",
            "Higher values mean connectivity is more concentrated in a small number of links.\nIf the single estimate stays higher than the bootstrap mean, it suggests sharper or less smoothed pathways.",
        ),
        (
            "Tail behavior through time",
            "Tracks the strongest-edge and top-k edge summaries through time.\nIt focuses on the upper tail where ecologically dominant pathways are most likely to appear.",
            "tail_metric_time_series.png",
            "Higher tail metrics mean stronger emphasis on dominant pathways.\nIf only the tail changes while broad correlations stay high, differences are concentrated in the strongest links rather than everywhere.",
        ),
        (
            "Top-edge strength distribution",
            "Compares the ranked strongest edge weights in the single matrix and bootstrap mean.\nIt makes it easy to see whether one product has a sharper high-end tail.",
            f"top_edge_strength_distribution_time_{time_tag}.png",
            "A heavier single tail suggests sharper pathways in the single estimate.\nCloser curves suggest the strongest pathways keep similar magnitude under bootstrap averaging.",
        ),
        (
            "Node strength uncertainty",
            "Shows bootstrap intervals for selected source and sink strengths.\nIt directly compares the single node value with the bootstrap range for important reefs.",
            f"node_strength_uncertainty_time_{time_tag}.png",
            "If the single point falls outside the interval, that node is less robust.\nWide intervals mean the reef’s apparent importance is uncertain even if its mean value is high.",
        ),
        (
            "Node CV distribution",
            "Shows how uncertain node strengths are across reefs using bootstrap coefficients of variation.\nIt helps separate consistently important reefs from highly unstable ones.",
            f"node_strength_cv_hist_time_{time_tag}.png",
            "Higher CV means the node’s importance is less stable.\nA broad distribution means some reefs are much more uncertain than others.",
        ),
        (
            "Top-link overlap through time",
            "Tracks how much the strongest pathways overlap between the single matrix and bootstrap mean.\nThis is often more ecologically meaningful than comparing all edges equally.",
            "topk_overlap_time_series.png",
            "Higher overlap means the dominant pathways are more robust.\nLower overlap means bootstraping changes which links you would identify as most important.",
        ),
        (
            "Top-link churn",
            "Shows which edges enter or leave the strongest-link set.\nIt identifies specific pathways that are gained or lost when moving from the single matrix to the bootstrap mean.",
            f"edge_churn_time_{time_tag}.png",
            "Large churn means the identity of dominant pathways is changing.\nIf the entering and leaving edges have similar weights, conclusions about the strongest pathways are less stable.",
        ),
        (
            "Stable-edge counts",
            "Counts edges that remain important across many bootstrap samples under different stability rules.\nThis summarizes how large the robust bootstrap-supported backbone is.",
            "stable_edge_counts.png",
            "More stable edges mean the bootstrap ensemble supports a clearer robust backbone.\nSharp drops under stricter frequency rules mean the strongest pathways are sensitive to sampling variation.",
        ),
        (
            "Stable-edge overlap",
            "Compares stable bootstrap edges with the single and bootstrap-mean strongest links.\nIt tests whether the pathways highlighted by each product are also stable across the bootstrap ensemble.",
            "stable_edge_overlap.png",
            "Higher overlap means the main pathways are robust across uncertainty rules.\nLow overlap means the visually strongest links may not be the most stable ones.",
        ),
        (
            "Entropy vs strength",
            "Shows whether strong nodes are diffuse or focused in where they connect.\nIt links total importance to how broadly or narrowly that importance is distributed across partners.",
            f"entropy_strength_time_{time_tag}.png",
            "High entropy means connections are spread across many partners.\nLow entropy at high strength suggests a reef is important because of a small number of dominant routes.",
        ),
        (
            "Dominance ratio comparison",
            "Compares how much each node is dominated by one main partner.\nThis helps distinguish broad exporters/importers from reefs driven by one dominant connection.",
            f"dominance_ratio_time_{time_tag}.png",
            "Higher dominance means one connection dominates the node’s pattern.\nIf dominance falls in the bootstrap mean, bootstraping is making connectivity appear more diffuse.",
        ),
        (
            "Effective-link distribution",
            "Shows the effective number of destinations or sources for nodes.\nIt translates entropy into an easier ecological interpretation: how many meaningful partners does a reef effectively have?",
            f"effective_links_distribution_time_{time_tag}.png",
            "Higher values mean more diffuse connectivity patterns.\nLower values mean connectivity is concentrated into a smaller set of meaningful partners.",
        ),
        (
            "Backbone overlap",
            "Compares a simple backbone proxy between the single and bootstrap-mean networks.\nThe backbone is defined transparently from thresholds and local row-share dominance rather than a complex black-box method.",
            "backbone_overlap_time_series.png",
            "Higher overlap means both products emphasize similar core links.\nLower overlap means the structural backbone you would report depends on the uncertainty treatment.",
        ),
        (
            "Coverage rates",
            "Summarizes whether single and bootstrap-mean system metrics fall inside bootstrap intervals.\nIt extends coverage beyond edges to broader quantities such as connectance, Gini, self-recruitment, and top-edge summaries.",
            "coverage_rates_time_series.png",
            "Higher outside-CI rates mean the bootstrap ensemble is adding important uncertainty information.\nIf the single estimate often falls outside the bootstrap range, it should not be treated as a typical realization.",
        ),
        (
            "Coverage examples",
            "Shows example bootstrap intervals for selected node-level metrics.\nIt gives a concrete view of which reefs look stable and which ones move outside the bootstrap-supported range.",
            f"coverage_examples_time_{time_tag}.png",
            "Points outside the interval indicate disagreement with bootstrap uncertainty.\nWide intervals mean high uncertainty even when the single point remains inside the interval.",
        ),
        (
            "Discrepancy vs variability",
            "Compares single-vs-mean differences with the normal variability among bootstrap samples.\nThis is a direct check of whether the single estimate is just another plausible realization or something more distinct.",
            "single_vs_bootstrap_variability.png",
            "If the single-vs-mean value is outside the bootstrap range, the single estimate is not just a typical bootstrap realization.\nIf it sits inside the bootstrap range, the single matrix is more consistent with ensemble variability.",
        ),
        (
            "Top-k probability",
            "Shows the probability that selected reefs are in the top-k set across bootstrap samples.\nThis turns rank stability into a direct probability of being an important source or sink.",
            f"topk_probability_time_{time_tag}.png",
            "Probabilities near one indicate very stable importance.\nIntermediate probabilities mean a reef may be important in some realizations but not others.",
        ),
        (
            "Rank uncertainty intervals",
            "Shows bootstrap rank intervals for important reefs.\nIt helps separate reefs with stable ranking from reefs whose rank changes widely across bootstrap samples.",
            f"rank_uncertainty_intervals_time_{time_tag}.png",
            "Wide intervals mean reef importance is uncertain.\nIf the single rank falls outside the interval, the single estimate is not well aligned with bootstrap rank uncertainty.",
        ),
        (
            "Ordination cloud",
            "Places bootstrap samples, the single estimate, and the bootstrap mean in low-dimensional space.\nThis gives a compact view of whether the single matrix sits inside or outside the bootstrap cloud.",
            f"ordination_time_{time_tag}.png",
            "If the single estimate sits inside the bootstrap cloud, it is more consistent with the ensemble.\nIf it sits apart from the cloud, bootstraping is adding information not captured by the single estimate.",
        ),
        (
            "Community similarity",
            "Tracks consensus connectivity-region structure, community count, and bootstrap agreement through time.\nThis turns broad structural interpretation into an explicit uncertainty analysis rather than a simple grouping proxy.",
            "community_similarity_time_series.png",
            "Higher ARI means the same large-scale regions keep appearing.\nLower consensus probability means the apparent community map is less stable across bootstrap samples.",
        ),
        (
            "Community consensus map",
            "Maps reef-level consensus connectivity regions for a representative time.\nIt is the broadest ecological structure output in the workflow.",
            f"community_consensus_time_{time_tag}.png",
            "Large coherent patches mean the network separates into stable connectivity regions.\nFragmented or mixed patterns mean broad structure is less stable than a mean matrix alone might suggest.",
        ),
        (
            "Community instability map",
            "Maps where community assignment changes across bootstrap samples.\nThis identifies reefs whose broad regional role is unstable.",
            f"community_stability_time_{time_tag}.png",
            "Higher values mean a reef changes connectivity-region membership more often.\nHotspots indicate unstable interpretation of large-scale network structure.",
        ),
        (
            "Community co-assignment heatmap",
            "Shows bootstrap co-assignment probabilities after ordering reefs by consensus community.\nIt condenses a full co-assignment matrix into a readable broad-scale structure plot.",
            f"community_coassignment_time_{time_tag}.png",
            "Bright within-block structure means strong community support.\nDiffuse boundaries mean the apparent regions blend together across bootstrap realizations.",
        ),
        (
            "Stepping-stone map",
            "Maps the reefs with the highest bridge score.\nThese are reefs that connect distinct flow regions rather than only having high raw strength.",
            f"bridge_map_time_{time_tag}.png",
            "Larger, brighter points indicate stronger bridge importance.\nIf top bridge reefs are not the same as top sources or sinks, ecological interpretation changes in a meaningful way.",
        ),
        (
            "Stepping-stone uncertainty",
            "Shows the probability that selected reefs remain in the top bridge set across bootstrap samples.\nIt turns bridge-reef inference into a direct uncertainty statement.",
            f"bridge_uncertainty_time_{time_tag}.png",
            "Probabilities near one indicate robust bridge status.\nIntermediate values mean a reef looks important in some realizations but not others.",
        ),
        (
            "Stepping-stone vs strength",
            "Compares bridge score with source and sink strength.\nThis tests whether bridge reefs are truly distinct from simple strong exporters or importers.",
            f"bridge_strength_scatter_time_{time_tag}.png",
            "A diffuse cloud rather than a tight line means bridge importance is not reducible to total strength alone.\nThat is the main scientific reason to include bridge analysis separately.",
        ),
        (
            "Bridge removal impact",
            "Ranks reefs by approximate loss of inter-anchor exchange if they are removed.\nIt gives an interpretable removal-impact summary without an intractable exact betweenness calculation.",
            f"bridge_removal_impact_time_{time_tag}.png",
            "Higher values mean the reef supports a larger share of cross-region exchange.\nThese are candidate stepping stones even when they are not top sources or sinks.",
        ),
        (
            "Distance-decay envelope",
            "Shows connectivity density by dispersal distance with bootstrap uncertainty bands.\nThis tests whether short, medium, and long-distance structure differ in uncertainty.",
            f"distance_decay_time_{time_tag}.png",
            "Wide uncertainty bands at long distance mean rare long-distance pathways are unstable.\nA large gap between the single estimate and bootstrap mean suggests averaging is changing tail interpretation.",
        ),
        (
            "Distance-dependent uncertainty",
            "Tracks the bootstrap coefficient of variation of connectivity density across distance bins.\nIt highlights where uncertainty is concentrated in distance space.",
            f"distance_uncertainty_time_{time_tag}.png",
            "Higher values mean stronger uncertainty at that distance scale.\nIf the right tail stays high, long-distance interpretation is especially sensitive to uncertainty treatment.",
        ),
        (
            "Long-distance fraction",
            "Tracks the fraction of total connectivity carried by long-distance links through time.\nThis is a compact tail summary for ecological dispersal interpretation.",
            "long_distance_fraction_time_series.png",
            "If the single estimate stays above the bootstrap mean, the single product emphasizes long-distance links more strongly.\nLarge bootstrap intervals mean long-distance conclusions should be treated cautiously.",
        ),
        (
            "Regional exchange heatmap",
            "Aggregates reef-level connectivity into region-to-region exchange.\nIt translates the reef network into a management-scale summary.",
            f"regional_exchange_time_{time_tag}.png",
            "Strong off-diagonal cells indicate large interregional exchange.\nDiagonal dominance means retention within regions remains more important than exchange among them.",
        ),
        (
            "Regional roles",
            "Compares regional export and import magnitudes.\nThis shows which regions behave as broad sources, sinks, or balanced exchangers.",
            f"regional_roles_time_{time_tag}.png",
            "Export bars above import bars indicate net exporters.\nClose bars indicate more balanced regional roles.",
        ),
        (
            "Regional pair intervals",
            "Shows bootstrap intervals for selected region-pair exchanges.\nIt makes regional uncertainty visible without scanning the full matrix table.",
            f"regional_pair_intervals_time_{time_tag}.png",
            "If the single point sits outside the interval, the regional conclusion is less robust.\nWide intervals mean the apparent importance of that regional pathway varies strongly across bootstrap samples.",
        ),
        (
            "Mean-bias summary",
            "Places mean-matrix and single-matrix interpretation metrics inside the bootstrap sample distribution.\nIt is the most direct answer to whether mean-based summaries look more stable than real realizations.",
            "mean_bias_summary.png",
            "Percentiles near 50 mean the reference is typical of bootstrap realizations.\nValues near the tails mean that summary is making the system look unusually stable or unusually extreme.",
        ),
        (
            "Spatial clustering summary",
            "Summarizes Moran-style spatial clustering of disagreement and instability metrics.\nIt tests whether uncertainty is geographically clustered rather than spatially random.",
            "spatial_clustering_summary.png",
            "Higher Moran's I means uncertainty hotspots are spatially clustered.\nThat suggests localized ecological interpretations may be less stable than broad system summaries imply.",
        ),
        (
            "Spatial hotspots",
            "Maps hotspots of rank uncertainty, bridge uncertainty, and community instability.\nIt identifies where disagreement concentrates in geographic space.",
            f"spatial_hotspots_time_{time_tag}.png",
            "Red hotspots indicate concentrated instability among neighboring reefs.\nThese are good candidates for follow-up ecological interpretation or targeted validation.",
        ),
        (
            "Variance decomposition",
            "Partitions interpretation-focused metrics into bootstrap, spawning-period, year, and family components.\nIt helps distinguish sampling uncertainty from biologically meaningful temporal or family change.",
            "variance_decomposition.png",
            "Large bootstrap components mean within-time uncertainty is a dominant driver.\nLarge year or family components mean broader biological differences dominate over bootstrap noise.",
        ),
        (
            "Family similarity",
            "Compares time-mean source, sink, and bridge structure among family products.\nIt tests whether family-specific biology is producing distinct connectivity systems.",
            "family_similarity_heatmap.png",
            "Higher similarity means families share broad structure.\nLower similarity means biology is changing ecological interpretation, not only uncertainty magnitude.",
        ),
        (
            "Family overlap",
            "Compares overlap in top sources, sinks, and stepping-stone reefs among families.\nThis focuses on whether the same reefs stay important across families.",
            "family_overlap_summary.png",
            "Higher overlap means the same reefs remain influential across families.\nLower overlap means family-specific biology shifts which reefs matter most.",
        ),
        (
            "Family ordination",
            "Places family summary vectors into a common low-dimensional space.\nIt provides a compact comparison of overall family structure.",
            "family_ordination.png",
            "Families that plot near each other have similar broad connectivity structure.\nFamilies that separate strongly are better treated as distinct systems.",
        ),
        (
            "Family uncertainty panel",
            "Compares available uncertainty summaries among family products.\nIt keeps single-only and bootstrap-capable families in the same comparison table while showing where uncertainty evidence is still missing.",
            "family_uncertainty_panel.png",
            "Higher edge or bridge CV means more uncertainty.\nHigher consensus probability means community structure is more stable.",
        ),
        (
            "Ecological threshold robustness",
            "Compares community, bridge, and backbone conclusions under several ecologically motivated threshold definitions.\nIt extends threshold sensitivity beyond simple numerical cutoffs.",
            "ecological_threshold_robustness.png",
            "Values near one mean the higher-level conclusion survives the threshold change.\nLower values mean ecological interpretation depends strongly on how links are filtered.",
        ),
    ]

    lines = [
        "# Analysis Catalog",
        "",
        f"Generated: {generated}",
        "",
        "This file is a quick guide to the analyses and figures in the connectivity comparison workflow.",
        "",
    ]
    for name, what_it_does, figure, how_to_read in entries:
        lines.extend(
            [
                f"## {name}",
                "",
                f"**What it does:** {what_it_does}",
                "",
                f"**Figure:** `{figure}`",
                "",
                f"**How to read it:** {how_to_read}",
                "",
            ]
        )

    lines.extend(
        [
            "## Recommended Figures For A Paper Or Report",
            "",
            "- `time_series_matrix_metrics.png` for the broad comparison.",
            "- `rank_correlation_time_series.png` and `rank_change_bars_time_XX.png` for source and sink stability.",
            "- `lorenz_curves_time_XX.png` and `concentration_time_series.png` for concentration of connectivity.",
            "- `stable_edge_overlap.png` and `single_vs_bootstrap_variability.png` for robustness and uncertainty.",
            "- `ordination_time_XX.png` for the bootstrap cloud view.",
            "- `community_similarity_time_series.png` and `community_consensus_time_XX.png` for large-scale structure.",
            "- `bridge_map_time_XX.png` and `bridge_strength_scatter_time_XX.png` for stepping-stone interpretation.",
            "- `regional_exchange_time_XX.png` and `mean_bias_summary.png` for management-scale and mean-bias interpretation.",
            "",
            "## Priority Ranking",
            "",
            "### Core Paper Analyses",
            "",
            "- `time_series_matrix_metrics.png`",
            "- `rank_correlation_time_series.png`",
            "- `community_similarity_time_series.png`",
            "- `bridge_map_time_XX.png`",
            "- `distance_decay_time_XX.png`",
            "- `regional_exchange_time_XX.png`",
            "- `mean_bias_summary.png`",
            "- `variance_decomposition.png`",
            "",
            "### Supporting Analyses",
            "",
            "- `stable_edge_overlap.png`",
            "- `single_vs_bootstrap_variability.png`",
            "- `community_coassignment_time_XX.png`",
            "- `bridge_uncertainty_time_XX.png`",
            "- `regional_pair_intervals_time_XX.png`",
            "- `long_distance_fraction_time_series.png`",
            "- `family_similarity_heatmap.png`",
            "- `ecological_threshold_robustness.png`",
            "",
            "### Optional Analyses",
            "",
            "- `heatmap_panel_time_XX.png`",
            "- `scatter_time_XX.png`",
            "- `entropy_strength_time_XX.png`",
            "- `effective_links_distribution_time_XX.png`",
            "- `spatial_hotspots_time_XX.png`",
            "- `family_ordination.png`",
            "",
        ]
    )
    return "\n".join(lines) + "\n"
