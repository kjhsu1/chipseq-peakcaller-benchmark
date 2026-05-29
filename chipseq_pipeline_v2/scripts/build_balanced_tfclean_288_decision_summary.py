"""
Build a decision-oriented summary from the local balanced TF-clean 288 progress table.
"""

"""Imports"""

import argparse
from pathlib import Path
from typing import List

import pandas as pd


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build a heartbeat-friendly decision summary for the local balanced TF-clean 288 wave."
    )
    parser.add_argument(
        "--progress-csv",
        type=Path,
        default=Path("analysis_outputs/tfclean_balanced_288_local_progress_current/progress_summary.csv"),
        help="Progress summary CSV from summarize_balanced_tfclean_288_progress.py",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis_outputs/tfclean_balanced_288_local_decision_current"),
        help="Directory to write the decision summary outputs",
    )
    return parser.parse_args()


def next_action(data: pd.DataFrame) -> str:
    """Choose the highest-value next operator action from current state."""
    if (data["launch_state"] == "score_ready").any():
        return "run_score_refresh"
    if (data["launch_state"] == "in_progress").any():
        if (data["launch_state"] == "not_started").any():
            return "defer_until_later_heartbeat_before_more_launches"
        return "defer_until_later_heartbeat"
    if (data["launch_state"] == "not_started").any():
        return "launch_when_capacity_allows"
    return "complete_or_compare"


def recommendation_lines(data: pd.DataFrame) -> List[str]:
    """Write a compact recommendation block from current states."""
    lines: List[str] = []
    launch_counts = data["launch_state"].value_counts().to_dict()
    phase_counts = data["phase"].value_counts().to_dict()
    lines.append(f"- next action: `{next_action(data)}`")
    lines.append(
        "- launch-state counts: "
        + ", ".join(f"`{key}`={value}" for key, value in sorted(launch_counts.items()))
    )
    lines.append(
        "- phase counts: "
        + ", ".join(f"`{key}`={value}" for key, value in sorted(phase_counts.items()))
    )

    active = data[data["launch_state"] == "in_progress"]
    if not active.empty:
        active_labels = ", ".join(f"`{label}`" for label in active["canonical_selector"].tolist())
        lines.append(f"- active roots: {active_labels}")

    waiting = data[data["launch_state"] == "not_started"]
    if not waiting.empty:
        waiting_labels = ", ".join(f"`{label}`" for label in waiting["canonical_selector"].tolist())
        lines.append(f"- not-started roots: {waiting_labels}")

    score_ready = data[data["launch_state"] == "score_ready"]
    if not score_ready.empty:
        ready_labels = ", ".join(f"`{label}`" for label in score_ready["canonical_selector"].tolist())
        lines.append(f"- score-ready roots: {ready_labels}")

    return lines


def write_readme(output_dir: Path, data: pd.DataFrame) -> None:
    """Write the decision summary markdown."""
    lines = [
        "# Local Balanced TF-Clean 288 Decision Summary",
        "",
        "This summary converts the current full-wave local `288` progress table",
        "into a single operator recommendation for heartbeat resumes.",
        "",
        "## Recommendation",
    ]
    lines.extend(recommendation_lines(data))
    lines.extend(
        [
            "",
            "## Per-Root Snapshot",
            "",
            "| canonical selector | launch state | phase | run dirs | treat BAMs | control BAMs | peak BEDs | recommended action |",
            "| --- | --- | --- | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in data.itertuples(index=False):
        lines.append(
            f"| {row.canonical_selector} | {row.launch_state} | {row.phase} | "
            f"{row.run_dirs}/{row.expected_runs} | {row.treat_bams}/{row.expected_runs} | "
            f"{row.control_bams}/{row.expected_runs} | {row.peak_beds}/{row.expected_runs} | "
            f"{row.recommended_action} |"
        )
    (output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Build the decision summary outputs."""
    args = parse_args()
    data = pd.read_csv(args.progress_csv)
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    data.to_csv(output_dir / "progress_snapshot.csv", index=False)
    write_readme(output_dir, data)
    print("\n".join(recommendation_lines(data)))


main()
