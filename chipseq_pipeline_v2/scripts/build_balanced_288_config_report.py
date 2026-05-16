"""
Build per-config PR/recall/F1 vs control-coverage plots and metadata summaries.
"""

"""Imports"""

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build a six-config balanced sweep report from peak_pr_stats outputs."
    )
    parser.add_argument(
        "--input-dirs",
        nargs="+",
        type=Path,
        required=True,
        help="Category report directories containing per_run_stats.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory to write the combined report overview files",
    )
    parser.add_argument(
        "--attempt-log",
        type=Path,
        default=None,
        help="Optional run-attempt log to copy into the report root",
    )
    return parser.parse_args()


def parse_manifest(path: Path) -> Dict[str, str]:
    """Parse a simple key: value manifest file."""
    data: Dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        data[key.strip()] = value.strip()
    return data


def summarize_parameter_values(data: pd.DataFrame) -> Tuple[List[str], List[str]]:
    """Return fixed and swept parameter summaries for info output."""
    ignore_cols = {
        "run_id",
        "total_called",
        "total_planted",
        "tp_called",
        "tp_planted",
        "precision",
        "recall",
        "f1",
    }
    fixed_lines: List[str] = []
    swept_lines: List[str] = []

    for col in sorted(data.columns):
        if col in ignore_cols:
            continue
        unique_values = sorted(data[col].dropna().astype(str).unique().tolist())
        if not unique_values:
            continue
        joined = ", ".join(unique_values)
        if len(unique_values) == 1:
            fixed_lines.append(f"- `{col}`: `{joined}`")
        else:
            swept_lines.append(f"- `{col}`: `{joined}`")

    return fixed_lines, swept_lines


def aggregate_plot_points(data: pd.DataFrame) -> pd.DataFrame:
    """Aggregate per-run stats into one plot point per control/treatment pair."""
    summary = (
        data.groupby(["coverage_treat", "coverage_ctrl"], as_index=False)
        .agg(
            tp_called=("tp_called", "sum"),
            total_called=("total_called", "sum"),
            tp_planted=("tp_planted", "sum"),
            total_planted=("total_planted", "sum"),
            n_runs=("run_id", "count"),
        )
        .sort_values(["coverage_treat", "coverage_ctrl"])
    )
    summary["coverage_treat"] = summary["coverage_treat"].astype(float)
    summary["coverage_ctrl"] = summary["coverage_ctrl"].astype(float)
    summary["precision"] = summary["tp_called"] / summary["total_called"]
    summary["precision"] = summary["precision"].fillna(0.0)
    summary["recall"] = summary["tp_planted"] / summary["total_planted"]
    summary["recall"] = summary["recall"].fillna(0.0)
    denom = summary["precision"] + summary["recall"]
    summary["f1"] = 0.0
    valid = denom > 0
    summary.loc[valid, "f1"] = (
        2 * summary.loc[valid, "precision"] * summary.loc[valid, "recall"] / denom.loc[valid]
    )
    return summary


def plot_summary(
    summary: pd.DataFrame,
    title: str,
    output_png: Path,
    log_x: bool,
) -> None:
    """Write a 3-panel precision/recall/F1 figure."""
    metrics = [("precision", "Precision"), ("recall", "Recall"), ("f1", "F1")]
    styles = {
        5.0: {"color": "tab:blue", "marker": "o", "linestyle": "-", "zorder": 3},
        10.0: {"color": "tab:orange", "marker": "s", "linestyle": "--", "zorder": 2},
        20.0: {"color": "tab:red", "marker": "^", "linestyle": "-.", "zorder": 1},
    }

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharex=True)
    for ax, (metric, label) in zip(axes, metrics):
        for treat_cov in sorted(summary["coverage_treat"].unique(), reverse=True):
            line = summary[summary["coverage_treat"] == treat_cov].sort_values("coverage_ctrl")
            style = styles.get(float(treat_cov), {"color": "black", "marker": "o", "linestyle": "-", "zorder": 1})
            ax.plot(
                line["coverage_ctrl"],
                line[metric],
                marker=style["marker"],
                linewidth=1.8,
                markersize=6.5,
                markerfacecolor="white",
                markeredgewidth=1.2,
                color=style["color"],
                linestyle=style["linestyle"],
                label=f"treat={int(treat_cov)}",
                zorder=style["zorder"],
            )
            ax.scatter(
                line["coverage_ctrl"],
                line[metric],
                s=20,
                color=style["color"],
                zorder=style["zorder"] + 0.1,
            )
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.25)
        ax.set_title(label)
        ax.set_xlabel("control coverage")
        if log_x:
            ax.set_xscale("log")

    axes[0].set_ylabel("metric value")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.99),
        ncol=max(1, len(labels)),
        frameon=False,
    )
    suffix = " (log x-axis)" if log_x else " (linear x-axis)"
    fig.suptitle(f"{title}{suffix}", y=1.06, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_info_file(
    output_path: Path,
    category_name: str,
    source_results_dir: str,
    manifest_data: Dict[str, str],
    fixed_lines: List[str],
    swept_lines: List[str],
    summary: pd.DataFrame,
) -> None:
    """Write a markdown overview of data sources and plotted point counts."""
    lines = [
        f"# {category_name}",
        "",
        "## Source Data",
        f"- archived results dir: `{source_results_dir}`",
        f"- params csv: `{manifest_data.get('params_csv', 'unknown')}`",
        f"- filter used by `peak_pr_stats.py`: `{manifest_data.get('filter', 'unknown')}`",
        f"- included per-run rows: `{manifest_data.get('included_runs', 'unknown')}`",
        "",
        "## Fixed Parameters",
    ]
    lines.extend(fixed_lines or ["- none"])
    lines.extend(
        [
            "",
            "## Swept Parameters",
        ]
    )
    lines.extend(swept_lines or ["- none"])
    lines.extend(
        [
            "",
            "## Plot Point Counts",
            "",
            "| coverage_treat | coverage_ctrl | n_runs | total_called | total_planted | precision | recall | f1 |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in summary.itertuples(index=False):
        lines.append(
            f"| {row.coverage_treat:g} | {row.coverage_ctrl:g} | {row.n_runs} | "
            f"{row.total_called} | {row.total_planted} | {row.precision:.4f} | "
            f"{row.recall:.4f} | {row.f1:.4f} |"
        )

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_root_readme(output_dir: Path, categories: List[Tuple[str, str]]) -> None:
    """Write a root README listing category folders and their contents."""
    lines = [
        "# Balanced 288 Config Report",
        "",
        "This directory contains one subdirectory per completed balanced 288 config.",
        "Each category folder includes:",
        "- `per_run_stats.csv` produced by `scripts/peak_pr_stats.py`",
        "- `group_summary.csv` produced by `scripts/peak_pr_stats.py`",
        "- `plot_point_summary.csv` with one aggregated plot point per treatment/control pair",
        "- `pr_recall_f1_vs_ctrl_coverage.png` with 3 panels and one curve per treatment coverage",
        "- `data_info.md` summarizing the source data, swept parameters, and point counts",
        "",
        "## Included Categories",
    ]
    for category_name, source_dir in categories:
        lines.append(f"- `{category_name}` from `{source_dir}`")
    lines.append("")
    lines.append("If present, `attempt_history.log` records the known Slurm/log history for these six config runs.")
    lines.append("")
    (output_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Build the combined balanced 288 config report."""
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    categories: List[Tuple[str, str]] = []
    for input_dir in args.input_dirs:
        per_run_path = input_dir / "per_run_stats.csv"
        manifest_path = input_dir / "run_filter_manifest.txt"
        if not per_run_path.exists():
            raise FileNotFoundError(f"Missing per_run_stats.csv in {input_dir}")
        if not manifest_path.exists():
            raise FileNotFoundError(f"Missing run_filter_manifest.txt in {input_dir}")

        data = pd.read_csv(per_run_path, dtype={"run_id": str})
        manifest_data = parse_manifest(manifest_path)
        summary = aggregate_plot_points(data)
        summary.to_csv(input_dir / "plot_point_summary.csv", index=False)

        category_name = input_dir.name
        source_results_dir = manifest_data.get("results_dir", "unknown")
        fixed_lines, swept_lines = summarize_parameter_values(data)
        plot_summary(
            summary,
            title=f"{category_name}\nPrecision/Recall/F1 vs control coverage",
            output_png=input_dir / "pr_recall_f1_vs_ctrl_coverage.png",
            log_x=False,
        )
        plot_summary(
            summary,
            title=f"{category_name}\nPrecision/Recall/F1 vs control coverage",
            output_png=input_dir / "pr_recall_f1_vs_ctrl_coverage_logx.png",
            log_x=True,
        )
        write_info_file(
            output_path=input_dir / "data_info.md",
            category_name=category_name,
            source_results_dir=source_results_dir,
            manifest_data=manifest_data,
            fixed_lines=fixed_lines,
            swept_lines=swept_lines,
            summary=summary,
        )
        categories.append((category_name, source_results_dir))

    if args.attempt_log is not None and args.attempt_log.exists():
        (args.output_dir / "attempt_history.log").write_text(
            args.attempt_log.read_text(encoding="utf-8"),
            encoding="utf-8",
        )

    write_root_readme(args.output_dir, categories)


main()
