"""Aggregate Realstudy v2 metrics and apply the preregistered sufficiency rule."""

"""Imports"""

import argparse
from pathlib import Path

import pandas as pd
import yaml


"""Constants"""

AGREEMENT_METRICS = [
    "anchor_peak_retention", "query_peak_concordance", "genomic_base_jaccard", "peak_count_ratio"
]


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse analysis arguments."""
    parser = argparse.ArgumentParser(description="Summarize Realstudy v2 control-depth stability.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-table", type=Path, required=True)
    parser.add_argument("--run-metrics", type=Path, required=True)
    parser.add_argument("--parent-libraries", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def metric_wide(run_metrics: pd.DataFrame, threshold: float = 0.0) -> pd.DataFrame:
    """Pivot one overlap-threshold slice into one row per run."""
    selected = run_metrics[pd.to_numeric(run_metrics["threshold"]) == threshold]
    return selected.pivot(index="run_id", columns="metric", values="value").reset_index()


def add_adjacent_gain(frame: pd.DataFrame) -> pd.DataFrame:
    """Add mean anchor-retention gain to the next tested control depth."""
    frame = frame.sort_values(["study_id", "control_coverage_x"]).copy()
    frame["next_depth_anchor_retention_gain"] = (
        frame.groupby("study_id")["mean_anchor_peak_retention"].shift(-1)
        - frame["mean_anchor_peak_retention"]
    )
    return frame


def build_depth_summary(runs: pd.DataFrame, run_metrics: pd.DataFrame, policy: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarize seeds by depth and record per-seed agreement passes."""
    sampled = runs[runs["run_type"] == "control_subsample"].copy()
    sampled["control_coverage_x"] = pd.to_numeric(sampled["control_coverage_x"])
    metrics = metric_wide(run_metrics)
    joined = sampled.merge(metrics, on="run_id", how="left", validate="one_to_one")
    missing = sorted(set(AGREEMENT_METRICS) - set(joined.columns))
    if missing:
        raise ValueError(f"Missing primary metrics: {', '.join(missing)}")
    joined["seed_pass"] = (
        (joined["anchor_peak_retention"] >= policy["minimum_anchor_retention"])
        & (joined["query_peak_concordance"] >= policy["minimum_query_concordance"])
        & (joined["genomic_base_jaccard"] >= policy["minimum_base_jaccard"])
        & joined["peak_count_ratio"].between(
            policy["minimum_peak_count_ratio"], policy["maximum_peak_count_ratio"], inclusive="both"
        )
    )
    aggregations = {metric: ["mean", "std", "min", "max"] for metric in AGREEMENT_METRICS}
    summary = joined.groupby(["study_id", "control_coverage_x"], as_index=False).agg(aggregations)
    summary.columns = [
        "_".join(str(part) for part in column if part).strip("_") if isinstance(column, tuple) else column
        for column in summary.columns
    ]
    seed_counts = joined.groupby(["study_id", "control_coverage_x"], as_index=False).agg(
        passing_seeds=("seed_pass", "sum"), evaluated_seeds=("seed_pass", "size")
    )
    summary = add_adjacent_gain(summary.merge(seed_counts, on=["study_id", "control_coverage_x"]))
    return joined, summary


def decide_sufficiency(summary: pd.DataFrame, parent_libraries: pd.DataFrame, policy: dict, minimum_parent: int) -> pd.DataFrame:
    """Choose the lowest sufficient depth or one explicit non-success outcome."""
    decisions = []
    for study_id, group in summary.groupby("study_id"):
        parent = parent_libraries[
            (parent_libraries["study_id"] == study_id) & (parent_libraries["role"] == "control")
        ]
        eligible = pd.to_numeric(parent["eligible_fragments"], errors="coerce")
        if parent.empty or eligible.isna().all() or eligible.max() < minimum_parent:
            decisions.append(
                {"study_id": study_id, "outcome": "parent_library_insufficient", "selected_coverage_x": "",
                 "passing_seeds": 0, "rationale": f"Eligible parent control fragments are below {minimum_parent}."}
            )
            continue
        group = group.sort_values("control_coverage_x").copy()
        group["mean_pass"] = (
            (group["anchor_peak_retention_mean"] >= policy["minimum_anchor_retention"])
            & (group["query_peak_concordance_mean"] >= policy["minimum_query_concordance"])
            & (group["genomic_base_jaccard_mean"] >= policy["minimum_base_jaccard"])
            & group["peak_count_ratio_mean"].between(
                policy["minimum_peak_count_ratio"], policy["maximum_peak_count_ratio"], inclusive="both"
            )
            & (group["passing_seeds"] >= policy["minimum_passing_seeds"])
        )
        sufficient = group[
            group["mean_pass"]
            & group["next_depth_anchor_retention_gain"].notna()
            & (group["next_depth_anchor_retention_gain"] <= policy["maximum_next_depth_retention_gain"])
        ]
        if not sufficient.empty:
            row = sufficient.iloc[0]
            outcome = "sufficient_depth_identified"
            selected = row["control_coverage_x"]
            rationale = "All agreement, seed, and next-depth plateau criteria passed."
        elif group["mean_pass"].iloc[-1] and pd.isna(group["next_depth_anchor_retention_gain"].iloc[-1]):
            row = group.iloc[-1]
            outcome = "agreement_without_plateau"
            selected = ""
            rationale = "Agreement passed only at the highest tested depth, so the next-depth plateau is untested."
        else:
            row = group.iloc[-1]
            outcome = "no_tested_depth_sufficient"
            selected = ""
            rationale = "No tested depth passed every preregistered agreement and plateau criterion."
        decisions.append(
            {"study_id": study_id, "outcome": outcome, "selected_coverage_x": selected,
             "passing_seeds": int(row["passing_seeds"]), "rationale": rationale}
        )
    return pd.DataFrame(decisions)


def main() -> None:
    """Write seed-level, depth-level, adjacent-gain, and decision tables."""
    args = parse_args()
    design = yaml.safe_load(args.config.read_text())["realstudy_v2"]
    runs = pd.read_csv(args.run_table, dtype=str).fillna("")
    metrics = pd.read_csv(args.run_metrics)
    parents = pd.read_csv(args.parent_libraries)
    seed_level, depth_summary = build_depth_summary(runs, metrics, design["sufficiency"])
    decisions = decide_sufficiency(
        depth_summary, parents, design["sufficiency"], int(design["minimum_eligible_control_fragments"])
    )
    seed_level = seed_level.sort_values(["study_id", "seed", "control_coverage_x"])
    seed_level["adjacent_anchor_retention_gain"] = (
        seed_level.groupby(["study_id", "seed"])["anchor_peak_retention"].shift(-1)
        - seed_level["anchor_peak_retention"]
    )
    adjacent_run_metrics = seed_level[["run_id", "adjacent_anchor_retention_gain"]].rename(
        columns={"adjacent_anchor_retention_gain": "value"}
    )
    adjacent_run_metrics["metric"] = "adjacent_anchor_retention_gain"
    adjacent_run_metrics["threshold"] = 0.0
    adjacent_run_metrics["status"] = "complete"
    augmented_metrics = pd.concat(
        [metrics, adjacent_run_metrics[["run_id", "metric", "threshold", "value", "status"]]],
        ignore_index=True,
    )
    sensitivity = runs[runs["run_type"] == "control_subsample"][["run_id", "study_id", "control_coverage_x", "seed"]].merge(
        metrics, on="run_id", how="inner"
    )
    sensitivity["control_coverage_x"] = pd.to_numeric(sensitivity["control_coverage_x"])
    overlap_summary = sensitivity.groupby(
        ["study_id", "control_coverage_x", "threshold", "metric"]
    )["value"].agg(["mean", "std", "min", "max"]).reset_index()
    policy_rows = []
    variants = {
        "relaxed": {**design["sufficiency"], "minimum_anchor_retention": 0.85, "minimum_query_concordance": 0.85, "minimum_base_jaccard": 0.75, "maximum_next_depth_retention_gain": 0.03},
        "registered": design["sufficiency"],
        "strict": {**design["sufficiency"], "minimum_anchor_retention": 0.95, "minimum_query_concordance": 0.95, "minimum_base_jaccard": 0.90, "minimum_peak_count_ratio": 0.95, "maximum_peak_count_ratio": 1.05, "maximum_next_depth_retention_gain": 0.01, "minimum_passing_seeds": 3},
    }
    for policy_name, policy in variants.items():
        policy_seed, policy_summary = build_depth_summary(runs, metrics, policy)
        policy_decisions = decide_sufficiency(
            policy_summary, parents, policy, int(design["minimum_eligible_control_fragments"])
        )
        policy_decisions["policy"] = policy_name
        policy_rows.append(policy_decisions)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    seed_level.to_csv(args.output_dir / "seed_level_metrics.csv", index=False)
    depth_summary.to_csv(args.output_dir / "depth_summary.csv", index=False)
    depth_summary[["study_id", "control_coverage_x", "next_depth_anchor_retention_gain"]].to_csv(
        args.output_dir / "adjacent_depth_gains.csv", index=False
    )
    decisions.to_csv(args.output_dir / "enough_control_decisions.csv", index=False)
    augmented_metrics.to_csv(args.output_dir / "run_metrics_augmented.csv", index=False)
    overlap_summary.to_csv(args.output_dir / "overlap_sensitivity_summary.csv", index=False)
    pd.concat(policy_rows, ignore_index=True).to_csv(args.output_dir / "policy_sensitivity.csv", index=False)


main()
