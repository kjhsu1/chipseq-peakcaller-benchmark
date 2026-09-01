"""Render the six preregistered Realstudy v2 paper figures and source data."""

"""Imports"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml


"""Constants"""

FORMATS = ["pdf", "svg", "png"]
STUDY_COLORS = {"encode_ceh27_young_adult": "#2364AA", "geo_h3k9me2_fer1_adult": "#D1495B"}


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse paper-figure arguments."""
    parser = argparse.ArgumentParser(description="Render Realstudy v2 paper figures.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--selection-audit", type=Path, required=True)
    parser.add_argument("--files", type=Path, required=True)
    parser.add_argument("--seed-level", type=Path, required=True)
    parser.add_argument("--depth-summary", type=Path, required=True)
    parser.add_argument("--adjacent-gains", type=Path, required=True)
    parser.add_argument("--seed-pairs", type=Path, required=True)
    parser.add_argument("--overlap-sensitivity", type=Path, required=True)
    parser.add_argument("--policy-sensitivity", type=Path, required=True)
    parser.add_argument("--stratified-metrics", type=Path, required=True)
    parser.add_argument("--locus-tracks", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def save_figure(figure: plt.Figure, output_dir: Path, figure_id: str, source: pd.DataFrame, caption: str) -> None:
    """Save vector/raster figures, panel source CSV, and caption stub."""
    figure_dir = output_dir / figure_id
    figure_dir.mkdir(parents=True, exist_ok=True)
    for extension in FORMATS:
        kwargs = {"dpi": 300} if extension == "png" else {}
        figure.savefig(figure_dir / f"{figure_id}.{extension}", bbox_inches="tight", **kwargs)
    source.to_csv(figure_dir / f"{figure_id}_source_data.csv", index=False)
    (figure_dir / f"{figure_id}_caption.md").write_text(caption.strip() + "\n")
    plt.close(figure)


def figure_one(config: dict, audit: pd.DataFrame, files: pd.DataFrame, output_dir: Path) -> None:
    """Render design, dataset selection, eligibility, and nested-sampling overview."""
    source = files.groupby(["study_id", "role"], as_index=False)["raw_reads"].sum()
    figure, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    pivot = source.pivot(index="study_id", columns="role", values="raw_reads") / 1_000_000
    pivot.plot.bar(ax=axes[0], color=["#72B7B2", "#4C78A8"])
    axes[0].set_ylabel("Raw reads (millions)")
    axes[0].set_title("Selected empirical libraries")
    coverages = np.asarray(config["control_coverages_x"], dtype=float)
    targets = np.rint(coverages * config["genome_size_bp"] / config["normalization_fragment_length_bp"]).astype(int)
    axes[1].step(coverages, targets / 1_000_000, where="mid", color="#F58518", linewidth=2)
    axes[1].scatter(coverages, targets / 1_000_000, color="#F58518")
    axes[1].set_xscale("log", base=2)
    axes[1].set_xlabel("Control depth (×)")
    axes[1].set_ylabel("Exact nested fragments (millions)")
    axes[1].set_title("Locked seven-level nested series")
    figure.tight_layout()
    design_source = pd.DataFrame({"control_coverage_x": coverages, "target_fragments": targets})
    design_source["panel"] = "nested_design"
    source["panel"] = "input_depth"
    audit_export = audit.copy()
    audit_export["panel"] = "selection_audit"
    save_figure(
        figure, output_dir, "figure_1_design",
        pd.concat([source, design_source, audit_export], ignore_index=True, sort=False),
        "Figure 1. Selected real C. elegans ChIP-seq datasets and the preregistered exact nested empirical-control design. Raw depth does not guarantee mapped-fragment eligibility; the workflow applies a separate 32× preflight.",
    )


def figure_two(depth: pd.DataFrame, output_dir: Path) -> None:
    """Render global dose-response curves for the four primary stability metrics."""
    metrics = ["anchor_peak_retention", "query_peak_concordance", "genomic_base_jaccard", "peak_count_ratio"]
    figure, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    for axis, metric in zip(axes.flat, metrics):
        for study_id, group in depth.groupby("study_id"):
            axis.plot(group["control_coverage_x"], group[f"{metric}_mean"], marker="o", label=study_id, color=STUDY_COLORS.get(study_id))
            std = group[f"{metric}_std"].fillna(0)
            axis.fill_between(group["control_coverage_x"], group[f"{metric}_mean"] - std, group[f"{metric}_mean"] + std, alpha=0.15, color=STUDY_COLORS.get(study_id))
        axis.set_xscale("log", base=2)
        axis.set_title(metric.replace("_", " ").title())
        axis.set_ylim(bottom=0)
    axes[0, 0].legend(fontsize=8)
    axes[1, 0].set_xlabel("Control depth (×)")
    axes[1, 1].set_xlabel("Control depth (×)")
    figure.tight_layout()
    save_figure(figure, output_dir, "figure_2_dose_response", depth, "Figure 2. Anchor-relative peak stability versus empirical control depth. Lines are seed means and ribbons are one standard deviation across deterministic sampling seeds.")


def figure_three(seed: pd.DataFrame, adjacent: pd.DataFrame, seed_pairs: pd.DataFrame, output_dir: Path) -> None:
    """Render seed variability and adjacent-depth convergence."""
    figure, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    for study_id, group in seed.groupby("study_id"):
        for seed_id, seed_group in group.groupby("seed"):
            axes[0].plot(seed_group["control_coverage_x"], seed_group["anchor_peak_retention"], marker="o", alpha=0.65, label=f"{study_id}: {seed_id}")
    axes[0].set_xscale("log", base=2)
    axes[0].set_title("Seed-specific persistence")
    axes[0].set_xlabel("Control depth (×)")
    axes[0].set_ylabel("Anchor-peak retention")
    for study_id, group in adjacent.groupby("study_id"):
        axes[1].plot(group["control_coverage_x"], group["next_depth_anchor_retention_gain"], marker="s", label=study_id, color=STUDY_COLORS.get(study_id))
    axes[1].axhline(0.02, color="black", linestyle="--", linewidth=1)
    axes[1].set_xscale("log", base=2)
    axes[1].set_title("Adjacent-depth convergence")
    axes[1].set_xlabel("Current control depth (×)")
    axes[1].set_ylabel("Retention gain at next depth")
    axes[1].legend(fontsize=8)
    pair_jaccard = seed_pairs[seed_pairs["metric"] == "genomic_base_jaccard"]
    for study_id, group in pair_jaccard.groupby("study_id"):
        axes[2].scatter(group["control_coverage_x"], group["value"], alpha=0.7, label=study_id, color=STUDY_COLORS.get(study_id))
    axes[2].set_xscale("log", base=2)
    axes[2].set_ylim(0, 1.02)
    axes[2].set_title("Pairwise seed agreement")
    axes[2].set_xlabel("Control depth (×)")
    axes[2].set_ylabel("Peak-base Jaccard")
    figure.tight_layout()
    seed_source = seed.copy(); seed_source["panel"] = "seed_persistence"
    adjacent_source = adjacent.copy(); adjacent_source["panel"] = "adjacent_gain"
    pair_source = seed_pairs.copy(); pair_source["panel"] = "seed_pair_agreement"
    save_figure(figure, output_dir, "figure_3_reproducibility", pd.concat([seed_source, adjacent_source, pair_source], ignore_index=True, sort=False), "Figure 3. Computational sampling uncertainty and convergence. Seeds represent deterministic resamples of the same empirical parent control, not biological replicates.")


def figure_four(seed: pd.DataFrame, output_dir: Path) -> None:
    """Render narrow summit and broad boundary/width/score behavior."""
    figure, axes = plt.subplots(2, 2, figsize=(12, 8))
    panels = [
        ("median_summit_displacement_bp", "Narrow summit displacement"),
        ("median_boundary_displacement_bp", "Broad boundary displacement"),
        ("median_width_displacement_bp", "Broad width displacement"),
        ("score_rank_correlation", "Matched-peak score rank"),
    ]
    plotted = []
    for axis, (metric, title) in zip(axes.flat, panels):
        if metric not in seed:
            axis.set_visible(False)
            continue
        for study_id, group in seed.dropna(subset=[metric]).groupby("study_id"):
            aggregate = group.groupby("control_coverage_x", as_index=False)[metric].mean()
            axis.plot(aggregate["control_coverage_x"], aggregate[metric], marker="o", color=STUDY_COLORS.get(study_id), label=study_id)
            plotted.append(aggregate.assign(study_id=study_id, metric=metric, panel=title))
        axis.set_xscale("log", base=2)
        axis.set_title(title)
        axis.set_xlabel("Control depth (×)")
    axes.flat[-1].legend(fontsize=8)
    figure.tight_layout()
    source = pd.concat(plotted, ignore_index=True, sort=False) if plotted else seed.iloc[0:0]
    save_figure(figure, output_dir, "figure_4_peak_geometry", source, "Figure 4. Peak-geometry and score stability relative to each study's full-control anchor: summit displacement for narrow ceh-27 peaks and boundary behavior for broad H3K9me2 domains.")


def figure_five(stratified: pd.DataFrame, output_dir: Path) -> None:
    """Render ontology/background-class response heatmaps and effect sizes."""
    required = {"stratum_value", "control_coverage_x", "value", "metric"}
    if not required.issubset(stratified.columns):
        raise ValueError(f"Stratified metrics require columns: {sorted(required)}")
    selected_metric = "anchor_peak_retention" if "anchor_peak_retention" in set(stratified["metric"]) else stratified["metric"].iloc[0]
    selected = stratified[stratified["metric"] == selected_metric]
    matrix = selected.pivot_table(index="stratum_value", columns="control_coverage_x", values="value", aggfunc="mean")
    figure, axes = plt.subplots(1, 2, figsize=(13, max(4, 0.45 * len(matrix))))
    image = axes[0].imshow(matrix.to_numpy(), aspect="auto", vmin=0, vmax=1, cmap="viridis")
    axes[0].set_xticks(range(len(matrix.columns)), matrix.columns)
    axes[0].set_yticks(range(len(matrix.index)), matrix.index)
    axes[0].set_xlabel("Control depth (×)")
    axes[0].set_title(f"{selected_metric.replace('_', ' ').title()} by genomic class")
    figure.colorbar(image, ax=axes[0], label="Metric value")
    effects = (matrix.max(axis=1) - matrix.min(axis=1)).sort_values()
    axes[1].barh(effects.index, effects.values, color="#59A14F")
    axes[1].set_xlabel("Max–min depth effect")
    axes[1].set_title("Class-specific response")
    figure.tight_layout()
    save_figure(figure, output_dir, "figure_5_stratified_response", selected, "Figure 5. Genomic ontology/background-class response to empirical control depth. The heatmap exposes heterogeneous convergence and class-specific failure modes.")


def figure_six(tracks: pd.DataFrame, output_dir: Path) -> None:
    """Render representative loci from actual fixed-treatment/control signal tracks and peak calls."""
    required = {"study_id", "locus_id", "track_label", "position", "value", "track_order"}
    if not required.issubset(tracks.columns):
        raise ValueError(f"Locus-track source requires columns: {sorted(required)}")
    loci = list(tracks["locus_id"].drop_duplicates())
    figure, axes = plt.subplots(len(loci), 1, figsize=(12, max(4, 3.2 * len(loci))), squeeze=False)
    for axis, locus_id in zip(axes.flat, loci):
        subset = tracks[tracks["locus_id"] == locus_id]
        for _, group in subset.groupby(["track_order", "track_label"]):
            axis.plot(group["position"], group["value"], label=group["track_label"].iloc[0])
        axis.set_title(locus_id)
        axis.set_xlabel("ce11 coordinate")
        axis.set_ylabel("Normalized signal / peak state")
        axis.legend(fontsize=7, ncol=3)
    figure.tight_layout()
    save_figure(figure, output_dir, "figure_6_representative_loci", tracks, "Figure 6. Representative ceh-27 and H3K9me2 loci showing the fixed treatment, full empirical control, sampled controls, and resulting MACS2 calls. Tracks are generated from the analyzed BAM/peak artifacts, not illustrative placeholders.")


def supplementary_figures(
    audit: pd.DataFrame, overlap: pd.DataFrame, policies: pd.DataFrame, output_dir: Path
) -> None:
    """Render candidate, overlap-threshold, and decision-policy sensitivity supplements."""
    supplement = output_dir / "supplementary"
    audit_plot = audit.copy()
    audit_plot["control_raw_reads"] = pd.to_numeric(audit_plot["control_raw_reads"], errors="coerce")
    figure, axis = plt.subplots(figsize=(9, 4.5))
    colors = ["#59A14F" if value == "included" else "#E15759" for value in audit_plot["inclusion_status"]]
    axis.bar(audit_plot["candidate_id"], audit_plot["control_raw_reads"] / 1_000_000, color=colors)
    axis.set_ylabel("Raw control reads (millions)")
    axis.tick_params(axis="x", rotation=20)
    axis.set_title("Candidate-study eligibility audit")
    figure.tight_layout()
    save_figure(figure, supplement, "supplementary_1_selection", audit, "Supplementary Figure 1. Candidate-study audit, including the shallow excluded SKN-1 control.")
    selected = overlap[overlap["metric"] == "anchor_peak_retention"]
    figure, axis = plt.subplots(figsize=(9, 5))
    for (study_id, threshold), group in selected.groupby(["study_id", "threshold"]):
        axis.plot(group["control_coverage_x"], group["mean"], marker="o", label=f"{study_id}; reciprocal={threshold:g}")
    axis.set_xscale("log", base=2)
    axis.set_xlabel("Control depth (×)")
    axis.set_ylabel("Mean anchor retention")
    axis.set_title("Overlap-definition sensitivity")
    axis.legend(fontsize=7)
    figure.tight_layout()
    save_figure(figure, supplement, "supplementary_2_overlap", overlap, "Supplementary Figure 2. Stability of conclusions under any-overlap, 10%, and 50% reciprocal-overlap definitions.")
    policy_plot = policies.copy()
    policy_plot["selected_coverage_numeric"] = pd.to_numeric(policy_plot["selected_coverage_x"], errors="coerce")
    figure, axis = plt.subplots(figsize=(8, 4.5))
    for index, (study_id, group) in enumerate(policy_plot.groupby("study_id")):
        positions = np.arange(len(group)) + index * 0.25
        axis.bar(positions, group["selected_coverage_numeric"].fillna(0), width=0.25, label=study_id)
    axis.set_xticks(np.arange(policy_plot["policy"].nunique()), list(policy_plot["policy"].drop_duplicates()))
    axis.set_ylabel("Selected depth (×; zero = no selected depth)")
    axis.set_title("Relaxed, registered, and strict policy sensitivity")
    axis.legend(fontsize=8)
    figure.tight_layout()
    save_figure(figure, supplement, "supplementary_3_policy", policies, "Supplementary Figure 3. Sufficiency outcomes under relaxed, registered, and strict policies; zero bars indicate an explicit non-selection outcome.")


def main() -> None:
    """Render all six figures from analysis outputs."""
    args = parse_args()
    design = yaml.safe_load(args.config.read_text())["realstudy_v2"]
    audit = pd.read_csv(args.selection_audit)
    files = pd.read_csv(args.files)
    seed = pd.read_csv(args.seed_level)
    depth = pd.read_csv(args.depth_summary)
    adjacent = pd.read_csv(args.adjacent_gains)
    seed_pairs = pd.read_csv(args.seed_pairs)
    overlap = pd.read_csv(args.overlap_sensitivity)
    policies = pd.read_csv(args.policy_sensitivity)
    stratified = pd.read_csv(args.stratified_metrics)
    tracks = pd.read_csv(args.locus_tracks)
    figure_one(design, audit, files, args.output_dir)
    figure_two(depth, args.output_dir)
    figure_three(seed, adjacent, seed_pairs, args.output_dir)
    figure_four(seed, args.output_dir)
    figure_five(stratified, args.output_dir)
    figure_six(tracks, args.output_dir)
    supplementary_figures(audit, overlap, policies, args.output_dir)


main()
