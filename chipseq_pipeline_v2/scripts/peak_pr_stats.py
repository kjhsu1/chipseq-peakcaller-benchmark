"""
Compute precision/recall/F1 for peak calls vs planted peak centers.
"""

"""Imports"""

import argparse
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

from scripts.eval_helpers import (
    aggregate_counts_summary,
    aggregate_mean_of_runs_summary,
    counts_to_metrics,
    metric_definition_lines,
    resolve_peak_path,
)

try:
    import pysam
except ImportError:  # pragma: no cover - optional runtime dependency
    pysam = None


"""
Constants
"""

DEFAULT_RESULTS_DIR = Path("archived_results/results_tf_gcacc_ctrl_sweep32_clean64")


"""
Data Structures
"""


"""
Functions
"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Compute precision/recall/F1 for peak calls vs planted centers."
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help="Root results directory containing run folders and params/run_params.csv",
    )
    parser.add_argument(
        "--params-csv",
        type=Path,
        default=None,
        help="Path to run_params.csv (defaults to {results_dir}/params/run_params.csv)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory (defaults to archived_results/pr_stats_YYYYMMDD_HHMMSS)",
    )
    parser.add_argument(
        "--group-by",
        type=str,
        default="coverage_ctrl",
        help="Column to group by for summary statistics",
    )
    return parser.parse_args()


def load_planted_centers(path: Path) -> Dict[str, Set[int]]:
    """Load planted peak centers from a 1-bp BED file."""
    centers: Dict[str, Set[int]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, start, end, *_ = line.rstrip().split("\t")
            centers.setdefault(chrom, set()).add(int(start))
    return centers


def load_called_intervals(path: Path) -> Dict[str, List[Tuple[int, int]]]:
    """Load called peak intervals from a narrowPeak BED-like file."""
    intervals: Dict[str, List[Tuple[int, int]]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *_ = line.rstrip().split("\t")
            if not start.isdigit() or not end.isdigit():
                continue
            intervals.setdefault(chrom, []).append((int(start), int(end)))
    return intervals


def interval_overlaps_any(center: int, intervals: List[Tuple[int, int]]) -> bool:
    """Return True if a center lies within any interval (inclusive start, exclusive end)."""
    for start, end in intervals:
        if start <= center < end:
            return True
    return False


def compute_overlap_stats(
    planted: Dict[str, Set[int]],
    called: Dict[str, List[Tuple[int, int]]],
) -> Tuple[int, int, int, int]:
    """Compute TP counts and totals for precision/recall."""
    total_called = sum(len(v) for v in called.values())
    total_planted = sum(len(v) for v in planted.values())
    tp_called = 0
    tp_planted = 0

    for chrom, intervals in called.items():
        centers = planted.get(chrom, set())
        for start, end in intervals:
            if any(start <= center < end for center in centers):
                tp_called += 1

    for chrom, centers in planted.items():
        intervals = called.get(chrom, [])
        for center in centers:
            if interval_overlaps_any(center, intervals):
                tp_planted += 1

    return tp_called, total_called, tp_planted, total_planted


def precision_recall_f1(
    tp_called: int, total_called: int, tp_planted: int, total_planted: int
) -> Tuple[float, float, float]:
    """Return precision, recall, and F1 given counts."""
    precision = tp_called / total_called if total_called else 0.0
    recall = tp_planted / total_planted if total_planted else 0.0
    if precision + recall == 0:
        f1 = 0.0
    else:
        f1 = 2 * (precision * recall) / (precision + recall)
    return precision, recall, f1


def resolve_params_path(results_dir: Path, params_csv: Path) -> Path:
    """Resolve params CSV path if none provided."""
    if params_csv is not None:
        return params_csv
    return results_dir / "params" / "run_params.csv"


def resolve_output_dir(output_dir: Path) -> Path:
    """Resolve output directory if none provided."""
    if output_dir is not None:
        return output_dir
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return Path("archived_results") / f"pr_stats_{timestamp}"


def filter_runs(params: pd.DataFrame) -> pd.DataFrame:
    """Filter to unbiased (gc+acc == 0) or both biases present."""
    unbiased = (params["gc_exp"] == 0) & (params["acc_exp"] == 0)
    both_biased = (params["gc_exp"] > 0) & (params["acc_exp"] > 0)
    return params[unbiased | both_biased].copy()


def truth_mode_for_row(row: pd.Series) -> str:
    """Return the truth mode used for a run."""
    macs2_mode = str(row.get("macs2_mode", "narrow"))
    peakcaller = str(row.get("peakcaller", "macs2"))
    if peakcaller == "epic2" or macs2_mode == "broad":
        return "interval"
    return "center"


def planted_truth_path(results_dir: Path, row: pd.Series) -> Path:
    """Return the planted truth path for a run."""
    run_id = row["run_id"]
    if truth_mode_for_row(row) == "interval":
        return results_dir / run_id / "treat" / "planted_peak_intervals.bed"
    return results_dir / run_id / "treat" / "planted_peaks.bed"


def bam_path(results_dir: Path, row: pd.Series, cond: str) -> Path:
    """Return BAM path for one run and condition."""
    return results_dir / row["run_id"] / row["aligner"] / cond / "aligned.sorted.bam"


def file_size_or_none(path: Path) -> Optional[int]:
    """Return file size when present."""
    if path.exists():
        return path.stat().st_size
    return None


def bam_mapped_reads_or_none(path: Path) -> Optional[int]:
    """Return mapped BAM read count when possible."""
    if not path.exists() or pysam is None:
        return None
    with pysam.AlignmentFile(path, "rb") as handle:
        return int(handle.mapped)


def write_manifest(
    output_dir: Path,
    params_csv: Path,
    results_dir: Path,
    group_by: str,
    run_count: int,
) -> None:
    """Write a short manifest describing the run."""
    manifest = [
        f"results_dir: {results_dir}",
        f"params_csv: {params_csv}",
        "filter: (gc_exp == 0 and acc_exp == 0) OR (gc_exp > 0 and acc_exp > 0)",
        f"group_by: {group_by}",
        "headline_aggregation: ratio_of_summed_counts",
        f"included_runs: {run_count}",
    ]
    (output_dir / "run_filter_manifest.txt").write_text("\n".join(manifest), encoding="utf-8")


def main() -> None:
    """Run precision/recall/F1 computations and write outputs."""
    args = parse_args()
    results_dir = args.results_dir
    params_csv = resolve_params_path(results_dir, args.params_csv)
    output_dir = resolve_output_dir(args.output_dir)

    params = pd.read_csv(params_csv, dtype={"run_id": str})
    filtered = filter_runs(params)

    per_run_rows = []
    for _, row in filtered.iterrows():
        run_id = row["run_id"]
        planted_path = planted_truth_path(results_dir, row)
        called_peak_path = resolve_peak_path(results_dir, run_id, str(row.get("peakcaller", "macs2")))
        planted = load_planted_centers(planted_path)
        called = load_called_intervals(called_peak_path)
        tp_called, total_called, tp_planted, total_planted = compute_overlap_stats(
            planted, called
        )
        precision, recall, f1 = counts_to_metrics(
            tp_called, total_called, tp_planted, total_planted
        )
        control_bam = bam_path(results_dir, row, "con")
        treatment_bam = bam_path(results_dir, row, "treat")
        per_run_rows.append(
            {
                "run_id": run_id,
                "total_called": total_called,
                "total_planted": total_planted,
                "tp_called": tp_called,
                "tp_planted": tp_planted,
                "precision": precision,
                "recall": recall,
                "f1": f1,
                "truth_mode": truth_mode_for_row(row),
                "planted_truth_path": str(planted_path),
                "called_peak_path": str(called_peak_path),
                "called_peak_size_bytes": file_size_or_none(called_peak_path),
                "control_bam_path": str(control_bam),
                "control_bam_size_bytes": file_size_or_none(control_bam),
                "control_mapped_reads": bam_mapped_reads_or_none(control_bam),
                "treat_bam_path": str(treatment_bam),
                "treat_bam_size_bytes": file_size_or_none(treatment_bam),
                "treat_mapped_reads": bam_mapped_reads_or_none(treatment_bam),
                **{col: row[col] for col in row.index if col not in {"run_id"}},
            }
        )

    per_run_df = pd.DataFrame(per_run_rows)
    output_dir.mkdir(parents=True, exist_ok=True)
    per_run_df.to_csv(output_dir / "per_run_stats.csv", index=False)

    group_summary = aggregate_counts_summary(per_run_df, [args.group_by])
    group_summary.to_csv(output_dir / "group_summary.csv", index=False)
    group_summary.to_csv(output_dir / "group_summary_counts_based.csv", index=False)
    mean_of_runs_summary = aggregate_mean_of_runs_summary(per_run_df, [args.group_by])
    mean_of_runs_summary.to_csv(output_dir / "group_summary_mean_of_runs.csv", index=False)

    write_manifest(output_dir, params_csv, results_dir, args.group_by, len(per_run_df))
    metric_definition = metric_definition_lines(
        truth_mode="center for narrow runs; interval for broad or epic2 runs",
        aggregation_rule="ratio_of_summed_counts for headline summaries; mean_of_runs emitted separately",
        caller_settings="caller-specific paths resolved from run metadata",
        evaluation_scope=f"grouped by `{args.group_by}`",
    )
    (output_dir / "metric_definition.md").write_text(
        "\n".join(metric_definition) + "\n",
        encoding="utf-8",
    )


main()
