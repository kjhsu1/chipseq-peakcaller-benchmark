"""
Summarize local balanced TF-clean 288-run MACS2 progress by config root.
"""

"""Imports"""

import argparse
from pathlib import Path
from typing import Dict, List

import pandas as pd


"""Constants"""

SELECTOR_MAP: Dict[str, List[str]] = {
    "all_six": [
        "balanced_tfclean_flatearth_peaks_broad_integrated_288",
        "balanced_tfclean_flatearth_plateaus_broad_integrated_288",
        "balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288",
        "balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288",
        "balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288",
        "balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288",
    ],
    "flatearth_pair": [
        "balanced_tfclean_flatearth_peaks_broad_integrated_288",
        "balanced_tfclean_flatearth_plateaus_broad_integrated_288",
    ],
    "wavy_pair": [
        "balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288",
        "balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288",
    ],
    "hilly_pair": [
        "balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288",
        "balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288",
    ],
    "flatearth_peak_narrow": ["balanced_tfclean_flatearth_peaks_broad_integrated_288"],
    "flatearth_plateau_broad": ["balanced_tfclean_flatearth_plateaus_broad_integrated_288"],
    "wavy_peak_narrow": ["balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288"],
    "hilly_peak_narrow": ["balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288"],
    "wavy_plateau_broad": ["balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288"],
    "hilly_plateau_broad": ["balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288"],
}

CANONICAL_SELECTOR_MAP: Dict[str, str] = {
    "balanced_tfclean_flatearth_peaks_broad_integrated_288": "flatearth_peak_narrow",
    "balanced_tfclean_flatearth_plateaus_broad_integrated_288": "flatearth_plateau_broad",
    "balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288": "wavy_peak_narrow",
    "balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288": "hilly_peak_narrow",
    "balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288": "wavy_plateau_broad",
    "balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288": "hilly_plateau_broad",
}


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Summarize local 288-run MACS2 progress for selected balanced TF-clean configs."
    )
    parser.add_argument(
        "selectors",
        nargs="*",
        default=["all_six"],
        help="Selectors or raw config-root names to summarize.",
    )
    parser.add_argument(
        "--results-prefix",
        type=Path,
        default=Path("."),
        help="Directory containing results_<config_name> roots.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory to write summary outputs.",
    )
    parser.add_argument(
        "--expected-runs",
        type=int,
        default=288,
        help="Expected run count per config root.",
    )
    return parser.parse_args()


def resolve_config_names(selectors: List[str]) -> List[str]:
    """Expand selectors into ordered config-root names."""
    resolved: List[str] = []
    for selector in selectors:
        expanded = SELECTOR_MAP.get(selector, [selector])
        for config_name in expanded:
            if config_name not in resolved:
                resolved.append(config_name)
    return resolved


def default_output_dir(selectors: List[str]) -> Path:
    """Choose a stable output directory for the requested selector set."""
    if selectors == ["all_six"]:
        return Path("analysis_outputs/tfclean_balanced_288_local_progress_current")
    label = "__".join(selectors)
    safe_label = "".join(ch if ch.isalnum() or ch in {"_", "-"} else "_" for ch in label)
    return Path(f"analysis_outputs/tfclean_balanced_288_local_progress_{safe_label}")


def count_matches(root: Path, pattern: str) -> int:
    """Count files matching a glob pattern under one results root."""
    return sum(1 for path in root.glob(pattern) if path.is_file())


def determine_phase(row: Dict[str, object]) -> str:
    """Assign a coarse workflow phase from milestone counts."""
    if not row["results_dir_exists"]:
        return "not_started"
    if row["peak_beds"] >= row["expected_runs"]:
        return "score_ready"
    if row["peak_beds"] > 0:
        return "peak_calling"
    if row["treat_bams"] > 0 or row["control_bams"] > 0:
        return "alignment"
    if row["treat_reads"] > 0 or row["control_reads"] > 0:
        return "read_generation"
    if row["run_dirs"] > 0 or row["params_csv_exists"]:
        return "initialized"
    return "empty_root"


def determine_launch_state(row: Dict[str, object]) -> str:
    """Assign an operator-level state for the local 1728-wave handoff."""
    if row["stats_dir_exists"]:
        return "scored"
    if row["score_ready"]:
        return "score_ready"
    if row["results_dir_exists"]:
        return "in_progress"
    return "not_started"


def determine_recommended_action(row: Dict[str, object]) -> str:
    """Suggest the next sensible operator action from current state."""
    launch_state = row["launch_state"]
    if launch_state == "not_started":
        return "launch_when_capacity_allows"
    if launch_state == "in_progress":
        return "defer_until_later_heartbeat"
    if launch_state == "score_ready":
        return "run_score_refresh"
    return "complete_or_compare"


def collect_progress(config_name: str, results_prefix: Path, expected_runs: int) -> Dict[str, object]:
    """Collect milestone counts for one local results root."""
    results_dir = results_prefix / f"results_{config_name}"
    params_csv = results_dir / "params" / f"{config_name}_run_params.csv"
    stats_dir = results_prefix / "analysis_outputs" / config_name
    row: Dict[str, object] = {
        "config_name": config_name,
        "canonical_selector": CANONICAL_SELECTOR_MAP.get(config_name, config_name),
        "results_dir": str(results_dir),
        "results_dir_exists": results_dir.is_dir(),
        "params_csv": str(params_csv),
        "params_csv_exists": params_csv.is_file(),
        "stats_dir": str(stats_dir),
        "stats_dir_exists": stats_dir.is_dir(),
        "expected_runs": expected_runs,
        "run_dirs": 0,
        "treat_reads": 0,
        "control_reads": 0,
        "treat_bams": 0,
        "control_bams": 0,
        "peak_beds": 0,
    }
    if results_dir.is_dir():
        row["run_dirs"] = sum(
            1
            for path in results_dir.iterdir()
            if path.is_dir() and path.name.isdigit() and len(path.name) == 4
        )
        row["treat_reads"] = count_matches(results_dir, "*/treat/reads_R1.fasta")
        row["control_reads"] = count_matches(results_dir, "*/con/reads_R1.fasta")
        row["treat_bams"] = count_matches(results_dir, "*/bowtie2/treat/aligned.sorted.bam")
        row["control_bams"] = count_matches(results_dir, "*/bowtie2/con/aligned.sorted.bam")
        row["peak_beds"] = count_matches(results_dir, "*/peaks/macs2/*_peaks.bed")
    row["phase"] = determine_phase(row)
    row["score_ready"] = bool(row["peak_beds"] >= expected_runs)
    row["peak_completion_fraction"] = float(row["peak_beds"]) / float(expected_runs)
    row["launch_state"] = determine_launch_state(row)
    row["recommended_action"] = determine_recommended_action(row)
    return row


def write_readme(output_dir: Path, summary_df: pd.DataFrame) -> None:
    """Write a compact markdown summary of selected config-root progress."""
    lines = [
        "# Local Balanced TF-Clean 288 Progress",
        "",
        "This directory summarizes milestone progress for selected local",
        "`results_balanced_tfclean_*_288` MACS2 roots.",
        "",
        "## Selected Roots",
    ]
    for row in summary_df.itertuples(index=False):
        lines.extend(
            [
                f"### `{row.config_name}`",
                f"- canonical selector: `{row.canonical_selector}`",
                f"- launch state: `{row.launch_state}`",
                f"- recommended action: `{row.recommended_action}`",
                f"- phase: `{row.phase}`",
                f"- results dir exists: `{row.results_dir_exists}`",
                f"- params csv exists: `{row.params_csv_exists}`",
                f"- stats dir exists: `{row.stats_dir_exists}`",
                f"- run dirs: `{row.run_dirs}/{row.expected_runs}`",
                f"- treat reads: `{row.treat_reads}/{row.expected_runs}`",
                f"- control reads: `{row.control_reads}/{row.expected_runs}`",
                f"- treat BAMs: `{row.treat_bams}/{row.expected_runs}`",
                f"- control BAMs: `{row.control_bams}/{row.expected_runs}`",
                f"- MACS2 peak BEDs: `{row.peak_beds}/{row.expected_runs}`",
                f"- score-ready: `{row.score_ready}`",
                "",
            ]
        )
    (output_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Run the local 288 progress summarizer."""
    args = parse_args()
    config_names = resolve_config_names(args.selectors)
    rows = [
        collect_progress(config_name, args.results_prefix, args.expected_runs)
        for config_name in config_names
    ]
    summary_df = pd.DataFrame(rows)
    output_dir = args.output_dir or default_output_dir(args.selectors)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(output_dir / "progress_summary.csv", index=False)
    write_readme(output_dir, summary_df)
    print(summary_df.to_string(index=False))


main()
