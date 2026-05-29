"""
Build a staged HOMER category report from completed 128-run stats directories.
"""

"""Imports"""

import argparse
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple


"""
Constants
"""

EXPECTED_HOMER_STATS: List[Tuple[str, str]] = [
    ("flatearth_peak_narrow", "homer_tfclean_flatearth_peak_narrow_128_stats"),
    ("flatearth_plateau_broad", "homer_tfclean_flatearth_plateau_broad_128_stats"),
    ("wavy_peak_narrow", "homer_tfclean_realistic_peaks_wavy_narrow_128_stats"),
    ("hilly_peak_narrow", "homer_tfclean_realistic_peaks_hilly_narrow_128_stats"),
    ("wavy_plateau_broad", "homer_tfclean_realistic_plateaus_wavy_broad_128_stats"),
    ("hilly_plateau_broad", "homer_tfclean_realistic_plateaus_hilly_broad_128_stats"),
]


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build the staged HOMER 128-run combined report from completed category stats."
    )
    parser.add_argument(
        "--analysis-root",
        type=Path,
        default=Path("analysis_outputs"),
        help="Directory containing HOMER stats output directories.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory to write the combined staged report.",
    )
    parser.add_argument(
        "--category-map",
        type=Path,
        default=Path("configs/category_name_map.yaml"),
        help="Canonical category label mapping file.",
    )
    parser.add_argument(
        "--require-all",
        action="store_true",
        help="Fail unless all six expected HOMER stats directories are present.",
    )
    return parser.parse_args()


def collect_available_inputs(analysis_root: Path) -> Tuple[List[Path], List[str]]:
    """Return available stats directories and missing canonical labels."""
    available: List[Path] = []
    missing: List[str] = []
    for canonical_label, stats_dir_name in EXPECTED_HOMER_STATS:
        stats_dir = analysis_root / stats_dir_name
        if stats_dir.exists():
            available.append(stats_dir)
        else:
            missing.append(canonical_label)
    return available, missing


def write_stage_manifest(
    output_dir: Path,
    analysis_root: Path,
    available_inputs: List[Path],
    missing_labels: List[str],
) -> None:
    """Write a short manifest describing staged-report completeness."""
    lines = [
        f"analysis_root: {analysis_root}",
        f"completed_category_count: {len(available_inputs)}",
        f"missing_category_count: {len(missing_labels)}",
        "completed_stats_dirs:",
    ]
    lines.extend(f"- {path.name}" for path in available_inputs)
    lines.append("missing_canonical_categories:")
    lines.extend(f"- {label}" for label in missing_labels)
    lines.append("")
    lines.append(
        "This manifest records whether the staged HOMER report is partial or complete."
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "stage_status.txt").write_text("\n".join(lines), encoding="utf-8")


def run_balanced_report_builder(
    input_dirs: List[Path], output_dir: Path, category_map: Path
) -> None:
    """Invoke the existing balanced report builder on the selected inputs."""
    command = [
        sys.executable,
        "scripts/build_balanced_288_config_report.py",
        "--output-dir",
        str(output_dir),
        "--category-map",
        str(category_map),
        "--input-dirs",
    ]
    command.extend(str(path) for path in input_dirs)
    subprocess.run(command, check=True)


def main() -> None:
    """Build the staged HOMER report root from currently completed category stats."""
    args = parse_args()
    available_inputs, missing_labels = collect_available_inputs(args.analysis_root)

    if args.require_all and missing_labels:
        missing_text = ", ".join(missing_labels)
        raise FileNotFoundError(
            f"Missing staged HOMER stats for: {missing_text}"
        )
    if not available_inputs:
        raise FileNotFoundError(
            f"No staged HOMER stats directories found under {args.analysis_root}"
        )

    run_balanced_report_builder(available_inputs, args.output_dir, args.category_map)
    write_stage_manifest(args.output_dir, args.analysis_root, available_inputs, missing_labels)


main()
