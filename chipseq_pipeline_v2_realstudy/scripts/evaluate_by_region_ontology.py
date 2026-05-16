"""Aggregate ontology-based evaluation tables and figures."""

"""Imports"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Evaluate performance by ontology class.")
    parser.add_argument("--input-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Aggregate performance metrics by ontology class and control depth."""
    args = parse_args()
    df = pd.read_csv(args.input_csv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = (
        df.groupby(["ontology_class", "coverage_ctrl"], as_index=False)[["precision", "recall", "f1"]]
        .mean()
        .sort_values(["ontology_class", "coverage_ctrl"])
    )
    per_run = (
        df.groupby(["run_id", "study_id", "coverage_ctrl"], as_index=False)[["precision", "recall", "f1"]]
        .mean()
        .sort_values(["study_id", "run_id"])
    )
    failure_mode_summary = (
        df.groupby(["failure_mode", "coverage_ctrl"], as_index=False)[["precision", "recall", "f1"]]
        .mean()
        .sort_values(["failure_mode", "coverage_ctrl"])
    )
    control_response = (
        summary.groupby("ontology_class", as_index=False)
        .agg(
            min_ctrl=("coverage_ctrl", "min"),
            max_ctrl=("coverage_ctrl", "max"),
            min_f1=("f1", "min"),
            max_f1=("f1", "max"),
        )
    )
    control_response["f1_gain"] = control_response["max_f1"] - control_response["min_f1"]

    per_run.to_csv(args.output_dir / "per_run_overall_metrics.csv", index=False)
    summary.to_csv(args.output_dir / "per_ontology_metrics.csv", index=False)
    control_response.to_csv(args.output_dir / "control_response_by_ontology.csv", index=False)
    failure_mode_summary.to_csv(args.output_dir / "failure_mode_metrics.csv", index=False)

    heatmap = summary.pivot(index="ontology_class", columns="coverage_ctrl", values="f1")
    heatmap.to_csv(args.output_dir / "ontology_f1_heatmap.csv")
    if not heatmap.empty:
        fig, ax = plt.subplots(figsize=(8, max(3, len(heatmap) * 0.5)))
        im = ax.imshow(heatmap.fillna(0.0).values, aspect="auto")
        ax.set_xticks(range(len(heatmap.columns)))
        ax.set_xticklabels([str(value) for value in heatmap.columns], rotation=45, ha="right")
        ax.set_yticks(range(len(heatmap.index)))
        ax.set_yticklabels(list(heatmap.index))
        ax.set_title("F1 by Ontology Class and Control Depth")
        fig.colorbar(im, ax=ax, label="F1")
        fig.tight_layout()
        fig.savefig(args.output_dir / "ontology_f1_heatmap.png", dpi=150)
        plt.close(fig)

    lines = [
        "# Ontology Evaluation Summary",
        "",
        f"- runs summarized: `{df['run_id'].nunique() if 'run_id' in df.columns else 0}`",
        f"- ontology classes: `{df['ontology_class'].nunique() if 'ontology_class' in df.columns else 0}`",
        "",
        "## Questions Addressed",
        "- weak-peak recall response by control depth is available in `control_response_by_ontology.csv`",
        "- false-positive-prone ontology classes are summarized in `failure_mode_metrics.csv`",
        "- control-depth sensitivity by ontology class is summarized in `per_ontology_metrics.csv` and `ontology_f1_heatmap.png`",
    ]
    (args.output_dir / "failure_mode_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


main()
