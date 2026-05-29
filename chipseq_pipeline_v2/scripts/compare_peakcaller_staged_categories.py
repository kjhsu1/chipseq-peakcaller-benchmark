"""Compare staged HOMER and MACS2 category summaries."""

"""Imports"""

import argparse
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

from scripts.category_summary_lib import parse_simple_yaml_map


"""Constants"""

EXPECTED_STAGE_PAIRS: List[Tuple[str, str, str]] = [
    (
        "flatearth_peak_narrow",
        "homer_tfclean_flatearth_peak_narrow_128_stats",
        "macs2_control_tfclean_flatearth_peak_narrow_128_stats",
    ),
    (
        "flatearth_plateau_broad",
        "homer_tfclean_flatearth_plateau_broad_128_stats",
        "macs2_control_tfclean_flatearth_plateau_broad_128_stats",
    ),
    (
        "wavy_peak_narrow",
        "homer_tfclean_realistic_peaks_wavy_narrow_128_stats",
        "macs2_control_tfclean_realistic_peaks_wavy_narrow_128_stats",
    ),
    (
        "hilly_peak_narrow",
        "homer_tfclean_realistic_peaks_hilly_narrow_128_stats",
        "macs2_control_tfclean_realistic_peaks_hilly_narrow_128_stats",
    ),
    (
        "wavy_plateau_broad",
        "homer_tfclean_realistic_plateaus_wavy_broad_128_stats",
        "macs2_control_tfclean_realistic_plateaus_wavy_broad_128_stats",
    ),
    (
        "hilly_plateau_broad",
        "homer_tfclean_realistic_plateaus_hilly_broad_128_stats",
        "macs2_control_tfclean_realistic_plateaus_hilly_broad_128_stats",
    ),
]


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Compare staged HOMER and MACS2 category stats."
    )
    parser.add_argument(
        "--analysis-root",
        type=Path,
        default=Path("analysis_outputs"),
        help="Directory containing staged stats roots.",
    )
    parser.add_argument(
        "--name-map",
        type=Path,
        default=Path("configs/category_name_map.yaml"),
        help="Canonical category mapping file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis_outputs")
        / f"peakcaller_comparison_{datetime.now().strftime('%Y%m%d')}",
        help="Directory to write comparison outputs.",
    )
    return parser.parse_args()


def pick_best_f1_row(summary_df: pd.DataFrame) -> pd.Series:
    """Return the best F1 row, tie-broken toward lower control coverage."""
    ordered = summary_df.sort_values(["f1", "coverage_ctrl"], ascending=[False, True])
    return ordered.iloc[0]


def category_mode(label: str) -> str:
    """Return broad/narrow mode from the canonical label."""
    return "broad" if "broad" in label else "narrow"


def best_row_from_stats(stats_dir: Path) -> Dict[str, object]:
    """Read one stats root and return its best-row headline metrics."""
    summary_df = pd.read_csv(stats_dir / "group_summary_counts_based.csv")
    best_row = pick_best_f1_row(summary_df)
    return {
        "stats_dir": stats_dir.name,
        "best_f1": float(best_row["f1"]),
        "best_precision": float(best_row["precision"]),
        "best_recall": float(best_row["recall"]),
        "best_control_coverage": float(best_row["coverage_ctrl"]),
    }


def build_comparison_rows(
    analysis_root: Path, mapping: Dict[str, str]
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Build completed and missing comparison rows."""
    completed_rows: List[Dict[str, object]] = []
    missing_rows: List[Dict[str, object]] = []

    for expected_label, homer_dir_name, macs2_dir_name in EXPECTED_STAGE_PAIRS:
        canonical_name = mapping.get(homer_dir_name, expected_label)
        homer_dir = analysis_root / homer_dir_name
        macs2_dir = analysis_root / macs2_dir_name
        if homer_dir.exists() and macs2_dir.exists():
            homer_best = best_row_from_stats(homer_dir)
            macs2_best = best_row_from_stats(macs2_dir)
            f1_gain = macs2_best["best_f1"] - homer_best["best_f1"]
            completed_rows.append(
                {
                    "canonical_name": canonical_name,
                    "mode": category_mode(canonical_name),
                    "homer_stats_dir": homer_best["stats_dir"],
                    "macs2_stats_dir": macs2_best["stats_dir"],
                    "homer_best_f1": homer_best["best_f1"],
                    "homer_best_precision": homer_best["best_precision"],
                    "homer_best_recall": homer_best["best_recall"],
                    "homer_best_control_coverage": homer_best["best_control_coverage"],
                    "macs2_best_f1": macs2_best["best_f1"],
                    "macs2_best_precision": macs2_best["best_precision"],
                    "macs2_best_recall": macs2_best["best_recall"],
                    "macs2_best_control_coverage": macs2_best["best_control_coverage"],
                    "macs2_minus_homer_f1": f1_gain,
                    "winner": "macs2" if f1_gain > 0 else "homer" if f1_gain < 0 else "tie",
                }
            )
        else:
            missing_rows.append(
                {
                    "canonical_name": canonical_name,
                    "missing_homer": not homer_dir.exists(),
                    "missing_macs2": not macs2_dir.exists(),
                    "expected_homer_stats_dir": homer_dir_name,
                    "expected_macs2_stats_dir": macs2_dir_name,
                }
            )

    completed_df = pd.DataFrame(completed_rows).sort_values("canonical_name")
    missing_df = pd.DataFrame(
        missing_rows,
        columns=[
            "canonical_name",
            "missing_homer",
            "missing_macs2",
            "expected_homer_stats_dir",
            "expected_macs2_stats_dir",
        ],
    ).sort_values("canonical_name")
    return completed_df, missing_df


def aggregate_mode_rows(comparison_df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate comparison behavior by broad/narrow mode."""
    rows: List[Dict[str, object]] = []
    for mode, group in comparison_df.groupby("mode"):
        rows.append(
            {
                "mode": mode,
                "completed_categories": int(len(group)),
                "mean_homer_best_f1": float(group["homer_best_f1"].mean()),
                "mean_macs2_best_f1": float(group["macs2_best_f1"].mean()),
                "mean_macs2_minus_homer_f1": float(group["macs2_minus_homer_f1"].mean()),
                "macs2_wins": int((group["winner"] == "macs2").sum()),
                "homer_wins": int((group["winner"] == "homer").sum()),
            }
        )
    return pd.DataFrame(rows).sort_values("mode")


def write_stage_status(
    output_dir: Path,
    analysis_root: Path,
    comparison_df: pd.DataFrame,
    missing_df: pd.DataFrame,
) -> None:
    """Write a compact comparison-completeness manifest."""
    lines = [
        f"analysis_root: {analysis_root}",
        f"completed_category_count: {len(comparison_df)}",
        f"missing_category_count: {len(missing_df)}",
        "completed_categories:",
    ]
    lines.extend(f"- {name}" for name in comparison_df["canonical_name"].tolist())
    lines.append("missing_categories:")
    lines.extend(f"- {name}" for name in missing_df["canonical_name"].tolist())
    lines.append("")
    lines.append("This manifest is generated by compare_peakcaller_staged_categories.py.")
    (output_dir / "stage_status.txt").write_text("\n".join(lines), encoding="utf-8")


def write_decision_summary(
    output_dir: Path,
    comparison_df: pd.DataFrame,
    mode_df: pd.DataFrame,
    missing_df: pd.DataFrame,
) -> None:
    """Write a compact staged HOMER-vs-MACS2 decision note."""
    lines = [
        "# Peakcaller Staged Decision",
        "",
        f"- completed paired categories: `{len(comparison_df)}`",
        f"- missing paired categories: `{len(missing_df)}`",
        "",
    ]
    if not comparison_df.empty:
        total_gain = float(comparison_df["macs2_minus_homer_f1"].mean())
        macs2_wins = int((comparison_df["winner"] == "macs2").sum())
        homer_wins = int((comparison_df["winner"] == "homer").sum())
        lines.extend(
            [
                "## Current Signal",
                f"- mean MACS2-minus-HOMER best F1 gain across completed categories: `{total_gain:.4f}`",
                f"- MACS2 category wins: `{macs2_wins}`",
                f"- HOMER category wins: `{homer_wins}`",
                "",
            ]
        )
        if not mode_df.empty:
            for row in mode_df.itertuples(index=False):
                lines.append(
                    f"- `{row.mode}`: HOMER mean best F1 `{row.mean_homer_best_f1:.4f}` vs "
                    f"MACS2 `{row.mean_macs2_best_f1:.4f}` "
                    f"(gain `{row.mean_macs2_minus_homer_f1:.4f}` across `{row.completed_categories}` categories)"
                )
        lines.extend(
            [
                "",
                "## Decision Read",
                "- matched staged controls currently favor MACS2 over HOMER on completed categories",
            ]
        )
        if not missing_df.empty:
            lines.append(
                "- the full six-category caller comparison should be refreshed after the remaining hilly stats finish"
            )
        else:
            lines.append(
                "- the six-category paired caller comparison is now complete for the current staged wave"
            )
        if not missing_df.empty:
            lines.extend(
                [
                    "",
                    "## Remaining Uncertainty",
                    "- this is still a partial paired verdict because one or more category pairs are unfinished",
                ]
            )
    else:
        lines.extend(
            [
                "## Current Signal",
                "- no paired HOMER-vs-MACS2 categories are complete yet",
            ]
        )
    (output_dir / "decision_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_readme(
    output_dir: Path,
    analysis_root: Path,
    comparison_df: pd.DataFrame,
    mode_df: pd.DataFrame,
    missing_df: pd.DataFrame,
) -> None:
    """Write a compact README for the comparison outputs."""
    lines = [
        "# Peakcaller Staged Comparison",
        "",
        f"- analysis root: `{analysis_root}`",
        f"- completed paired categories: `{len(comparison_df)}`",
        f"- missing paired categories: `{len(missing_df)}`",
        "",
        "## Completed Paired Categories",
    ]
    if comparison_df.empty:
        lines.append("- none")
    else:
        for row in comparison_df.itertuples(index=False):
            lines.append(
                f"- `{row.canonical_name}`: HOMER best F1 `{row.homer_best_f1:.4f}` vs "
                f"MACS2 `{row.macs2_best_f1:.4f}` "
                f"(gain `{row.macs2_minus_homer_f1:.4f}`; winner `{row.winner}`)"
            )
    lines.extend(["", "## Missing Paired Categories"])
    if missing_df.empty:
        lines.append("- none")
    else:
        for row in missing_df.itertuples(index=False):
            lines.append(
                f"- `{row.canonical_name}`: missing HOMER=`{row.missing_homer}`, missing MACS2=`{row.missing_macs2}`"
            )
    if not mode_df.empty:
        lines.extend(["", "## Mode Headlines"])
        for row in mode_df.itertuples(index=False):
            lines.append(
                f"- `{row.mode}`: HOMER mean best F1 `{row.mean_homer_best_f1:.4f}` vs "
                f"MACS2 `{row.mean_macs2_best_f1:.4f}` "
                f"(gain `{row.mean_macs2_minus_homer_f1:.4f}`)"
            )
    lines.extend(
        [
            "",
            "## Output Tables",
            "- `category_comparison.csv`",
            "- `missing_categories.csv`",
            "- `mode_comparison.csv`",
            "- `decision_summary.md`",
        ]
    )
    (output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Write staged HOMER-vs-MACS2 comparison outputs."""
    args = parse_args()
    mapping = parse_simple_yaml_map(args.name_map)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    comparison_df, missing_df = build_comparison_rows(args.analysis_root, mapping)
    mode_df = aggregate_mode_rows(comparison_df) if not comparison_df.empty else pd.DataFrame()

    comparison_df.to_csv(args.output_dir / "category_comparison.csv", index=False)
    missing_df.to_csv(args.output_dir / "missing_categories.csv", index=False)
    if not mode_df.empty:
        mode_df.to_csv(args.output_dir / "mode_comparison.csv", index=False)
    write_stage_status(args.output_dir, args.analysis_root, comparison_df, missing_df)
    write_decision_summary(args.output_dir, comparison_df, mode_df, missing_df)
    write_readme(args.output_dir, args.analysis_root, comparison_df, mode_df, missing_df)


main()
