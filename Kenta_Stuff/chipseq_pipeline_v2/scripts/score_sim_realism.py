"""Score simulated BAMs against local real-data reference ranges."""

"""Imports"""

import argparse
from datetime import datetime
from pathlib import Path
from typing import List

import pandas as pd

from scripts.realism_metrics_lib import (
    interval_lengths,
    load_bed_intervals,
    mean_metaprofile,
    metaprofile_metrics,
    score_distance,
    summarize_depths,
    estimate_depths,
)


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build a BAM-only realism scorecard for simulated and real-reference runs."
    )
    parser.add_argument(
        "--input-manifest",
        type=Path,
        default=Path("analysis_outputs/realism_scorecard_20260513/input_manifest.csv"),
        help="CSV manifest listing simulated and reference BAM/peak inputs",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis_outputs") / f"realism_scorecard_{datetime.now().strftime('%Y%m%d')}",
        help="Output directory for realism scorecard tables",
    )
    parser.add_argument(
        "--half-window",
        type=int,
        default=250,
        help="Half-window size for summit-centered metaprofiles",
    )
    return parser.parse_args()


def compute_row_metrics(row: pd.Series, half_window: int) -> dict:
    """Compute BAM-only realism metrics for one manifest row."""
    peak_intervals = load_bed_intervals(Path(row["peak_bed"]))
    peak_depths = estimate_depths(Path(row["treatment_bam"]), peak_intervals)
    peak_stats = summarize_depths(peak_depths)
    lengths = interval_lengths(peak_intervals)
    profile = mean_metaprofile(Path(row["treatment_bam"]), peak_intervals, half_window)
    profile_stats = metaprofile_metrics(profile)
    control_ratio = 0.0
    if pd.notna(row.get("control_bam")) and str(row.get("control_bam")).strip():
        control_depths = estimate_depths(Path(row["control_bam"]), peak_intervals)
        control_mean = float(control_depths.mean()) if control_depths.size else 0.0
        control_ratio = peak_stats.mean_depth / max(control_mean, 1.0)
    return {
        "entry_id": row["entry_id"],
        "entry_type": row["entry_type"],
        "signal_class": row["signal_class"],
        "selection_status": row["selection_status"],
        "peak_count": int(len(peak_intervals)),
        "peak_width_mean": float(lengths.mean()) if lengths.size else 0.0,
        "peak_width_median": float(pd.Series(lengths).median()) if lengths.size else 0.0,
        "treatment_peak_mean_depth": peak_stats.mean_depth,
        "treatment_peak_cv": peak_stats.cv_depth,
        "treatment_peak_gini": peak_stats.gini_depth,
        "treatment_peak_bumpiness": peak_stats.bumpiness,
        "treatment_control_peak_ratio": control_ratio,
        **profile_stats,
    }


def main() -> None:
    """Run the BAM-only realism scorecard."""
    args = parse_args()
    manifest = pd.read_csv(args.input_manifest)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    available_rows = manifest[
        manifest["selection_status"].isin(["selected", "local_reference_available"])
        & manifest["treatment_bam"].notna()
        & manifest["peak_bed"].notna()
    ].copy()

    metrics_rows: List[dict] = []
    for row in available_rows.to_dict(orient="records"):
        metrics_rows.append(compute_row_metrics(pd.Series(row), args.half_window))

    metrics_df = pd.DataFrame(metrics_rows)
    if metrics_df.empty:
        (args.output_dir / "README.md").write_text(
            "# Realism Scorecard\n\nNo local BAM+peak inputs were available to score.\n",
            encoding="utf-8",
        )
        return

    sim_df = metrics_df[metrics_df["entry_type"] == "simulated"].copy()
    ref_df = metrics_df[metrics_df["entry_type"] == "real_reference"].copy()
    metric_names = [
        "peak_count",
        "peak_width_mean",
        "treatment_peak_mean_depth",
        "treatment_peak_cv",
        "treatment_peak_gini",
        "treatment_peak_bumpiness",
        "summit_height",
        "half_max_width",
        "auc",
        "symmetry",
        "shoulder_score",
    ]

    comparison_rows = []
    for _, sim_row in sim_df.iterrows():
        candidate_refs = ref_df[ref_df["signal_class"] == sim_row["signal_class"]]
        for _, ref_row in candidate_refs.iterrows():
            comparison_rows.append(
                {
                    "sim_entry_id": sim_row["entry_id"],
                    "ref_entry_id": ref_row["entry_id"],
                    "signal_class": sim_row["signal_class"],
                    "distance_to_reference": score_distance(sim_row, ref_row, metric_names),
                }
            )

    metrics_df.to_csv(args.output_dir / "run_metrics.csv", index=False)
    metrics_df.groupby(["entry_type", "signal_class"], as_index=False)[metric_names].mean().to_csv(
        args.output_dir / "group_metric_summary.csv", index=False
    )
    pd.DataFrame(comparison_rows).sort_values(
        ["signal_class", "distance_to_reference", "sim_entry_id"]
    ).to_csv(args.output_dir / "distance_to_reference.csv", index=False)

    lines = [
        "# Realism Scorecard",
        "",
        f"- manifest: `{args.input_manifest}`",
        f"- scored rows: `{len(metrics_df)}`",
        "- output tables separate simulated metrics, reference metrics, and distance-to-reference summaries",
        "- this phase is BAM only and does not use BigWig inputs",
    ]
    (args.output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


main()
