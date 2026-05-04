"""
Investigate which planted peaks were missed and which called peaks were false positives.
"""

"""Imports"""

import argparse
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Trace per-run recall and precision failures back to planted/called peak locations."
    )
    parser.add_argument(
        "--input-dirs",
        nargs="+",
        type=Path,
        required=True,
        help="Archived result directories containing run folders and params/run_params.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory to write investigation CSVs and README",
    )
    return parser.parse_args()


def resolve_output_dir(output_dir: Optional[Path]) -> Path:
    """Return explicit output directory or a timestamped default."""
    if output_dir is not None:
        return output_dir
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return Path("analysis_outputs") / f"peak_recovery_investigation_{timestamp}"


def category_name_from_path(path: Path) -> str:
    """Strip timestamp suffix from archived result directory name."""
    parts = path.name.split("_")
    if len(parts) >= 3 and parts[-2].isdigit() and parts[-1].isdigit():
        return "_".join(parts[:-2])
    return path.name


def filter_runs(params: pd.DataFrame) -> pd.DataFrame:
    """Use the same bias-family filter as peak_pr_stats.py."""
    unbiased = (params["gc_exp"] == 0) & (params["acc_exp"] == 0)
    both_biased = (params["gc_exp"] > 0) & (params["acc_exp"] > 0)
    return params[unbiased | both_biased].copy()


def resolve_peak_path(results_dir: Path, run_id: str, peakcaller: str) -> Path:
    """Resolve the called peak file for a run."""
    if peakcaller == "epic2":
        return results_dir / run_id / "peaks" / "epic2" / f"{run_id}_domains.bed"
    macs2_dir = results_dir / run_id / "peaks" / "macs2"
    candidate_paths = [
        macs2_dir / f"{run_id}_peaks.bed",
        macs2_dir / f"{run_id}_peaks.narrowPeak",
        macs2_dir / f"{run_id}_peaks.broadPeak",
    ]
    return next((path for path in candidate_paths if path.exists()), candidate_paths[0])


def load_planted_centers(path: Path) -> List[Dict]:
    """Load planted 1-bp peak centers from BED."""
    rows: List[Dict] = []
    with path.open("r", encoding="utf-8") as handle:
        for idx, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *_ = line.rstrip().split("\t")
            rows.append(
                {
                    "planted_index": idx,
                    "planted_chrom": chrom,
                    "planted_start": int(start),
                    "planted_end": int(end),
                    "planted_center": int(start),
                }
            )
    return rows


def load_called_intervals(path: Path) -> List[Dict]:
    """Load called peak intervals from BED-like output."""
    rows: List[Dict] = []
    with path.open("r", encoding="utf-8") as handle:
        for idx, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *rest = line.rstrip().split("\t")
            rows.append(
                {
                    "called_index": idx,
                    "called_chrom": chrom,
                    "called_start": int(start),
                    "called_end": int(end),
                    "called_width": int(end) - int(start),
                    "called_extra_fields": "\t".join(rest),
                }
            )
    return rows


def interval_contains_center(interval: Dict, planted: Dict) -> bool:
    """Return whether one called interval contains one planted center."""
    if interval["called_chrom"] != planted["planted_chrom"]:
        return False
    return interval["called_start"] <= planted["planted_center"] < interval["called_end"]


def distance_center_to_interval(center: int, start: int, end: int) -> int:
    """Return distance from a point to an interval, with zero for overlap."""
    if start <= center < end:
        return 0
    if center < start:
        return start - center
    return center - end + 1


def nearest_called_interval(planted: Dict, called: List[Dict]) -> Tuple[Optional[Dict], Optional[int]]:
    """Return nearest called interval on the same chromosome and its distance."""
    candidates = [row for row in called if row["called_chrom"] == planted["planted_chrom"]]
    if not candidates:
        return None, None
    nearest = min(
        candidates,
        key=lambda row: distance_center_to_interval(
            planted["planted_center"], row["called_start"], row["called_end"]
        ),
    )
    distance = distance_center_to_interval(
        planted["planted_center"], nearest["called_start"], nearest["called_end"]
    )
    return nearest, distance


def nearest_planted_center(interval: Dict, planted: List[Dict]) -> Tuple[Optional[Dict], Optional[int]]:
    """Return nearest planted center on the same chromosome and its distance."""
    candidates = [row for row in planted if row["planted_chrom"] == interval["called_chrom"]]
    if not candidates:
        return None, None
    nearest = min(
        candidates,
        key=lambda row: distance_center_to_interval(
            row["planted_center"], interval["called_start"], interval["called_end"]
        ),
    )
    distance = distance_center_to_interval(
        nearest["planted_center"], interval["called_start"], interval["called_end"]
    )
    return nearest, distance


def metric_values(tp_called: int, total_called: int, tp_planted: int, total_planted: int) -> Tuple[float, float, float]:
    """Return precision, recall, and F1."""
    precision = tp_called / total_called if total_called else 0.0
    recall = tp_planted / total_planted if total_planted else 0.0
    if precision + recall == 0:
        return precision, recall, 0.0
    return precision, recall, 2 * precision * recall / (precision + recall)


def parameter_payload(row: pd.Series) -> Dict:
    """Return run parameters to attach to detail rows."""
    keys = [
        "genome",
        "acc_key",
        "gc_key",
        "fragment_length",
        "read_length",
        "nb_k",
        "aligner",
        "peakcaller",
        "macs2_mode",
        "tf_exp",
        "seed",
        "tf_seed",
        "map_seed",
        "gc_exp",
        "acc_exp",
        "map_coverage_pct",
        "map_sigma",
        "map_enrich",
        "map_exp",
        "use_control",
        "coverage_ctrl",
        "coverage_treat",
        "tf_peak_count_treat",
        "tf_sigma",
        "tf_enrich",
    ]
    return {key: row[key] for key in keys if key in row.index}


def analyze_run(results_dir: Path, category: str, row: pd.Series) -> Tuple[Dict, List[Dict], List[Dict]]:
    """Analyze one run and return summary, missed planted rows, and false-positive rows."""
    run_id = row["run_id"]
    params = parameter_payload(row)
    planted_path = results_dir / run_id / "treat" / "planted_peaks.bed"
    called_path = resolve_peak_path(results_dir, run_id, str(row.get("peakcaller", "macs2")))
    planted = load_planted_centers(planted_path)
    called = load_called_intervals(called_path)

    matched_planted_indexes = set()
    matched_called_indexes = set()
    for called_row in called:
        for planted_row in planted:
            if interval_contains_center(called_row, planted_row):
                matched_called_indexes.add(called_row["called_index"])
                matched_planted_indexes.add(planted_row["planted_index"])

    tp_called = len(matched_called_indexes)
    tp_planted = len(matched_planted_indexes)
    precision, recall, f1 = metric_values(tp_called, len(called), tp_planted, len(planted))
    base = {
        "category": category,
        "source_dir": str(results_dir),
        "run_id": run_id,
        "planted_bed": str(planted_path),
        "called_peak_path": str(called_path),
        "total_called": len(called),
        "total_planted": len(planted),
        "tp_called": tp_called,
        "tp_planted": tp_planted,
        "false_positive_called": len(called) - tp_called,
        "missed_planted": len(planted) - tp_planted,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "has_perfect_precision": precision == 1.0,
        "has_perfect_recall": recall == 1.0,
    }
    summary = {**base, **params}

    missed_rows: List[Dict] = []
    for planted_row in planted:
        if planted_row["planted_index"] in matched_planted_indexes:
            continue
        nearest, distance = nearest_called_interval(planted_row, called)
        nearest_payload = {}
        if nearest is not None:
            nearest_payload = {
                "nearest_called_index": nearest["called_index"],
                "nearest_called_chrom": nearest["called_chrom"],
                "nearest_called_start": nearest["called_start"],
                "nearest_called_end": nearest["called_end"],
                "nearest_called_distance_bp": distance,
            }
        missed_rows.append({**summary, **planted_row, **nearest_payload})

    false_positive_rows: List[Dict] = []
    for called_row in called:
        if called_row["called_index"] in matched_called_indexes:
            continue
        nearest, distance = nearest_planted_center(called_row, planted)
        nearest_payload = {}
        if nearest is not None:
            nearest_payload = {
                "nearest_planted_index": nearest["planted_index"],
                "nearest_planted_chrom": nearest["planted_chrom"],
                "nearest_planted_center": nearest["planted_center"],
                "nearest_planted_distance_bp": distance,
            }
        false_positive_rows.append({**summary, **called_row, **nearest_payload})

    return summary, missed_rows, false_positive_rows


def write_readme(output_dir: Path, input_dirs: List[Path], category_summary: pd.DataFrame) -> None:
    """Write a compact README explaining outputs and matching criteria."""
    lines = [
        "# Peak Recovery Investigation",
        "",
        "Matching criterion: a called peak matches a planted peak when the planted 1-bp center lies inside the called interval (`called_start <= planted_center < called_end`).",
        "There is no 50 bp or 100 bp padding, and broad/narrow peaks use the same criterion.",
        "",
        "## Input Directories",
    ]
    for path in input_dirs:
        lines.append(f"- `{path}`")
    lines.extend(
        [
            "",
            "## Output Files",
            "- `per_run_recovery_summary.csv`: one row per run with precision/recall/F1 and run parameters.",
            "- `category_summary.csv`: counts of perfect/non-perfect precision and recall by category.",
            "- `missed_planted_peaks.csv`: one row per planted center not covered by any called peak, with seed/parameters and nearest called interval.",
            "- `false_positive_called_peaks.csv`: one row per called peak containing no planted center, with seed/parameters and nearest planted center.",
            "",
            "## Category Summary",
            "",
            "| category | runs | perfect_recall_runs | imperfect_recall_runs | perfect_precision_runs | imperfect_precision_runs |",
            "| --- | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in category_summary.itertuples(index=False):
        lines.append(
            f"| {row.category} | {row.runs} | {row.perfect_recall_runs} | "
            f"{row.imperfect_recall_runs} | {row.perfect_precision_runs} | "
            f"{row.imperfect_precision_runs} |"
        )
    (output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Run the peak recovery investigation."""
    args = parse_args()
    output_dir = resolve_output_dir(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    summaries: List[Dict] = []
    missed_planted_rows: List[Dict] = []
    false_positive_rows: List[Dict] = []

    for results_dir in args.input_dirs:
        category = category_name_from_path(results_dir)
        params = pd.read_csv(results_dir / "params" / "run_params.csv", dtype={"run_id": str})
        for _, row in filter_runs(params).iterrows():
            summary, missed, false_positives = analyze_run(results_dir, category, row)
            summaries.append(summary)
            missed_planted_rows.extend(missed)
            false_positive_rows.extend(false_positives)

    per_run = pd.DataFrame(summaries)
    per_run.to_csv(output_dir / "per_run_recovery_summary.csv", index=False)
    pd.DataFrame(missed_planted_rows).to_csv(output_dir / "missed_planted_peaks.csv", index=False)
    pd.DataFrame(false_positive_rows).to_csv(output_dir / "false_positive_called_peaks.csv", index=False)

    category_summary = (
        per_run.groupby("category", as_index=False)
        .agg(
            runs=("run_id", "count"),
            perfect_recall_runs=("has_perfect_recall", "sum"),
            perfect_precision_runs=("has_perfect_precision", "sum"),
            max_recall=("recall", "max"),
            max_precision=("precision", "max"),
            mean_recall=("recall", "mean"),
            mean_precision=("precision", "mean"),
            total_missed_planted=("missed_planted", "sum"),
            total_false_positive_called=("false_positive_called", "sum"),
        )
        .sort_values("category")
    )
    category_summary["imperfect_recall_runs"] = category_summary["runs"] - category_summary["perfect_recall_runs"]
    category_summary["imperfect_precision_runs"] = category_summary["runs"] - category_summary["perfect_precision_runs"]
    category_summary.to_csv(output_dir / "category_summary.csv", index=False)
    write_readme(output_dir, args.input_dirs, category_summary)


main()
