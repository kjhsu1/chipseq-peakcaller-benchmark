"""Measure peak-set reproducibility between deterministic control-sampling seeds."""

"""Imports"""

import argparse
import itertools
from pathlib import Path

import pandas as pd


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse seed-pair comparison arguments."""
    parser = argparse.ArgumentParser(description="Compare Realstudy v2 peak calls between seeds.")
    parser.add_argument("--run-table", type=Path, required=True)
    parser.add_argument("--peak-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def read_intervals(path: Path) -> dict[str, list[tuple[int, int]]]:
    """Read and merge a BED peak set by chromosome."""
    grouped = {}
    if path.stat().st_size == 0:
        return grouped
    frame = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2])
    for chromosome, group in frame.groupby(0):
        merged = []
        for start, end in group[[1, 2]].sort_values(1).itertuples(index=False, name=None):
            start, end = int(start), int(end)
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        grouped[str(chromosome)] = [(start, end) for start, end in merged]
    return grouped


def total_bp(intervals: dict[str, list[tuple[int, int]]]) -> int:
    """Return merged interval width."""
    return sum(end - start for values in intervals.values() for start, end in values)


def intersection_bp(left: dict, right: dict) -> int:
    """Return intersection width between two merged interval sets."""
    total = 0
    for chromosome in set(left) & set(right):
        a, b = left[chromosome], right[chromosome]
        i = j = 0
        while i < len(a) and j < len(b):
            total += max(0, min(a[i][1], b[j][1]) - max(a[i][0], b[j][0]))
            if a[i][1] <= b[j][1]:
                i += 1
            else:
                j += 1
    return total


def main() -> None:
    """Write pairwise base-Jaccard and peak-count-ratio metrics for all seed pairs."""
    args = parse_args()
    runs = pd.read_csv(args.run_table, dtype=str).fillna("")
    runs = runs[runs["run_type"] == "control_subsample"]
    rows = []
    for (study_id, coverage), group in runs.groupby(["study_id", "control_coverage_x"]):
        records = list(group.to_dict(orient="records"))
        for left, right in itertools.combinations(records, 2):
            left_intervals = read_intervals(args.peak_root / left["run_id"] / "peaks.bed")
            right_intervals = read_intervals(args.peak_root / right["run_id"] / "peaks.bed")
            left_bp, right_bp = total_bp(left_intervals), total_bp(right_intervals)
            shared = intersection_bp(left_intervals, right_intervals)
            union = left_bp + right_bp - shared
            values = {
                "genomic_base_jaccard": shared / union if union else 1.0,
                "peak_count_ratio": min(sum(map(len, left_intervals.values())), sum(map(len, right_intervals.values())))
                / max(1, max(sum(map(len, left_intervals.values())), sum(map(len, right_intervals.values())))),
            }
            for metric, value in values.items():
                rows.append(
                    {"study_id": study_id, "control_coverage_x": coverage, "seed_a": int(left["seed"]),
                     "seed_b": int(right["seed"]), "metric": metric, "value": value}
                )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False)


main()
