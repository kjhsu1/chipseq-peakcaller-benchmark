"""Score ChIPs simulated peak recovery over template real-study regions."""

"""Imports"""

import argparse
import csv
import math
from pathlib import Path

import pysam


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build per-region signal summaries for downstream ontology classification."
    )
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--run-table-csv", type=Path, required=True)
    parser.add_argument("--truth-bed", type=Path, required=True)
    parser.add_argument("--called-bed", type=Path, required=True)
    parser.add_argument("--treat-bam", type=Path, required=True)
    parser.add_argument("--control-bam", type=Path, required=True)
    parser.add_argument("--score-column", type=int, default=5)
    parser.add_argument("--output-csv", type=Path, required=True)
    return parser.parse_args()


def load_run_row(run_table_csv: Path, run_id: str) -> dict:
    """Return the run table row matching the requested run id."""
    with run_table_csv.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("run_id") == run_id:
                return row
    raise ValueError(f"Run id not found in run table: {run_id}")


def parse_bed(path: Path, score_column: int = 5) -> list[dict]:
    """Load BED-like intervals with an optional one-based score column."""
    intervals = []
    score_index = max(0, score_column - 1)
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) < 3:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            if end <= start:
                continue
            score = 0.0
            if len(fields) > score_index:
                try:
                    score = float(fields[score_index])
                except ValueError:
                    score = 0.0
            intervals.append(
                {
                    "chrom": fields[0],
                    "start": start,
                    "end": end,
                    "score": score,
                }
            )
    return intervals


def overlaps_any(interval: dict, intervals_by_chrom: dict[str, list[dict]]) -> bool:
    """Return whether one interval overlaps any interval on the same chromosome."""
    chrom_intervals = intervals_by_chrom.get(interval["chrom"], [])
    for other in chrom_intervals:
        if other["end"] <= interval["start"]:
            continue
        if other["start"] >= interval["end"]:
            return False
        return True
    return False


def group_by_chrom(intervals: list[dict]) -> dict[str, list[dict]]:
    """Group sorted intervals by chromosome."""
    grouped: dict[str, list[dict]] = {}
    for interval in sorted(intervals, key=lambda row: (row["chrom"], row["start"], row["end"])):
        grouped.setdefault(interval["chrom"], []).append(interval)
    return grouped


def mean(values: list[float]) -> float:
    """Return the arithmetic mean of a list, or zero for empty input."""
    if not values:
        return 0.0
    return sum(values) / len(values)


def sample_std(values: list[float]) -> float:
    """Return the population standard deviation of a list, or zero for empty input."""
    if not values:
        return 0.0
    avg = mean(values)
    return math.sqrt(sum((value - avg) ** 2 for value in values) / len(values))


def region_counts(bam: pysam.AlignmentFile, chrom: str, start: int, end: int) -> tuple[int, list[int]]:
    """Count overlapping reads and simple read-center bins for one genomic region."""
    width = max(1, end - start)
    bins = max(1, min(50, math.ceil(width / 100)))
    binned = [0 for _ in range(bins)]
    count = 0
    try:
        reads = bam.fetch(chrom, start, end)
    except ValueError:
        return 0, binned
    for read in reads:
        if read.is_unmapped or read.reference_start is None or read.reference_end is None:
            continue
        if read.reference_end <= start or read.reference_start >= end:
            continue
        count += 1
        center = max(start, min(end - 1, (read.reference_start + read.reference_end) // 2))
        bin_index = min(bins - 1, int(((center - start) / width) * bins))
        binned[bin_index] += 1
    return count, binned


def z_from_bins(values: list[int]) -> float:
    """Return a max-bin z-score from binned read counts."""
    avg = mean([float(value) for value in values])
    std = sample_std([float(value) for value in values])
    if std == 0.0:
        return 0.0
    return (max(values) - avg) / std


def cv_from_bins(values: list[int]) -> float:
    """Return coefficient of variation from binned read counts."""
    avg = mean([float(value) for value in values])
    if avg == 0.0:
        return 0.0
    return sample_std([float(value) for value in values]) / avg


def numeric(row: dict, name: str, default: str = "0") -> str:
    """Return a run-table field as a stable string value."""
    value = str(row.get(name, default)).strip()
    return value if value else default


def write_scores(args: argparse.Namespace) -> None:
    """Write per-template-region ontology scoring input table."""
    run_row = load_run_row(args.run_table_csv, args.run_id)
    truth_regions = parse_bed(args.truth_bed, args.score_column)
    called_by_chrom = group_by_chrom(parse_bed(args.called_bed, args.score_column))
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "run_id",
        "study_id",
        "coverage_treat",
        "coverage_ctrl",
        "seed",
        "macs2_mode",
        "chrom",
        "start",
        "end",
        "width_bp",
        "template_score",
        "called_overlap",
        "treat_reads",
        "control_reads",
        "treat_rpk",
        "control_rpk",
        "ctrl_cv",
        "ctrl_bump_z",
        "treat_peak_z",
        "log2_enrichment",
        "plateau_width_bp",
        "precision",
        "recall",
        "f1",
        "metric_basis",
    ]
    with pysam.AlignmentFile(str(args.treat_bam), "rb") as treat_bam:
        with pysam.AlignmentFile(str(args.control_bam), "rb") as control_bam:
            with args.output_csv.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                for region in truth_regions:
                    chrom = region["chrom"]
                    start = int(region["start"])
                    end = int(region["end"])
                    width = max(1, end - start)
                    treat_reads, treat_bins = region_counts(treat_bam, chrom, start, end)
                    control_reads, control_bins = region_counts(control_bam, chrom, start, end)
                    called_overlap = 1.0 if overlaps_any(region, called_by_chrom) else 0.0
                    treat_rpk = treat_reads / (width / 1000.0)
                    control_rpk = control_reads / (width / 1000.0)
                    log2_enrichment = math.log2((treat_reads + 1.0) / (control_reads + 1.0))
                    writer.writerow(
                        {
                            "run_id": args.run_id,
                            "study_id": numeric(run_row, "study_id", "unknown"),
                            "coverage_treat": numeric(run_row, "coverage_treat"),
                            "coverage_ctrl": numeric(run_row, "coverage_ctrl"),
                            "seed": numeric(run_row, "seed"),
                            "macs2_mode": numeric(run_row, "macs2_mode", "unknown"),
                            "chrom": chrom,
                            "start": start,
                            "end": end,
                            "width_bp": width,
                            "template_score": region["score"],
                            "called_overlap": called_overlap,
                            "treat_reads": treat_reads,
                            "control_reads": control_reads,
                            "treat_rpk": round(treat_rpk, 6),
                            "control_rpk": round(control_rpk, 6),
                            "ctrl_cv": round(cv_from_bins(control_bins), 6),
                            "ctrl_bump_z": round(z_from_bins(control_bins), 6),
                            "treat_peak_z": round(z_from_bins(treat_bins), 6),
                            "log2_enrichment": round(log2_enrichment, 6),
                            "plateau_width_bp": width,
                            "precision": called_overlap,
                            "recall": called_overlap,
                            "f1": called_overlap,
                            "metric_basis": "truth_template_region_recovery_proxy",
                        }
                    )


def main() -> None:
    """Run the scoring workflow."""
    write_scores(parse_args())


main()
