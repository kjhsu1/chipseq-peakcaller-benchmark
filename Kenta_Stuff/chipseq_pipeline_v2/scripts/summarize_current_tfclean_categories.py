"""Summarize the current TF-clean six-category controlled baseline."""

"""Imports"""

import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd

from scripts.category_summary_lib import (
    best_control_by_mean_f1,
    canonical_label,
    expected_run_count,
    observed_parameter_values,
    parse_simple_yaml_map,
    precision_dropoff_points,
    recall_saturation_points,
    summarize_counts_by_control,
    warning_lines,
    zero_metric_row_count,
)


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Summarize the current TF-clean six-category balanced report."
    )
    parser.add_argument(
        "--report-root",
        type=Path,
        default=Path("analysis_outputs/tfclean_balanced_288_config_report_20260505"),
        help="Directory containing one subdirectory per TF-clean category",
    )
    parser.add_argument(
        "--name-map",
        type=Path,
        default=Path("configs/category_name_map.yaml"),
        help="Canonical category mapping file",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis_outputs") / f"controlled_simulator_cleanup_{datetime.now().strftime('%Y%m%d')}",
        help="Output directory for cleanup summary tables",
    )
    return parser.parse_args()


def main() -> None:
    """Build cleanup summary tables for the current controlled baseline."""
    args = parse_args()
    mapping = parse_simple_yaml_map(args.name_map)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    category_rows = []
    count_tables = []
    saturation_tables = []
    dropoff_tables = []

    for category_dir in sorted(args.report_root.glob("balanced_tfclean_*_288")):
        per_run_df = pd.read_csv(category_dir / "per_run_stats.csv")
        plot_df = pd.read_csv(category_dir / "plot_point_summary.csv")
        label = canonical_label(category_dir.name, mapping)
        params = observed_parameter_values(per_run_df)
        category_rows.append(
            {
                "config_name": category_dir.name,
                "canonical_name": label,
                "observed_runs": len(per_run_df),
                "expected_runs": expected_run_count(per_run_df),
                "zero_metric_rows": zero_metric_row_count(per_run_df),
                "best_control_mean_f1": best_control_by_mean_f1(per_run_df),
                "swept_parameters": "; ".join(
                    f"{key}={','.join(values)}" for key, values in sorted(params.items())
                ),
                "warnings": "; ".join(warning_lines()),
            }
        )

        counts_df = summarize_counts_by_control(per_run_df)
        counts_df.insert(0, "canonical_name", label)
        count_tables.append(counts_df)

        saturation_df = recall_saturation_points(plot_df)
        saturation_df.insert(0, "canonical_name", label)
        saturation_tables.append(saturation_df)

        dropoff_df = precision_dropoff_points(plot_df)
        dropoff_df.insert(0, "canonical_name", label)
        dropoff_tables.append(dropoff_df)

    summary_df = pd.DataFrame(category_rows).sort_values("canonical_name")
    summary_df.to_csv(args.output_dir / "category_summary.csv", index=False)
    pd.concat(count_tables, ignore_index=True).to_csv(
        args.output_dir / "mean_counts_by_control.csv", index=False
    )
    pd.concat(saturation_tables, ignore_index=True).to_csv(
        args.output_dir / "recall_saturation_points.csv", index=False
    )
    pd.concat(dropoff_tables, ignore_index=True).to_csv(
        args.output_dir / "precision_dropoff_points.csv", index=False
    )

    lines = [
        "# Controlled Simulator Cleanup",
        "",
        f"- source report root: `{args.report_root}`",
        f"- category count: `{len(summary_df)}`",
        "",
        "## Standard Warnings",
    ]
    lines.extend(f"- {line}" for line in warning_lines())
    lines.extend(
        [
            "",
            "## Output Tables",
            "- `category_summary.csv`",
            "- `mean_counts_by_control.csv`",
            "- `recall_saturation_points.csv`",
            "- `precision_dropoff_points.csv`",
        ]
    )
    (args.output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


main()
