"""One-time bias calibration helper for GC/accessibility settings."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Calibrate GC/ACC bias parameters")
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--params-csv", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--accessibility-bed",
        type=Path,
        default=Path("data/accessibility/ce11_1pct/bedA.bed"),
    )
    parser.add_argument("--target-gc-q90-q10", type=float, default=5.5)
    parser.add_argument("--target-open-closed", type=float, default=2.0)
    return parser.parse_args()


def load_accessibility_intervals(path: Path) -> Dict[str, List[Tuple[int, int]]]:
    intervals: Dict[str, List[Tuple[int, int]]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *_ = line.rstrip().split("\t")
            intervals.setdefault(chrom, []).append((int(start), int(end)))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def in_open_region(pos: int, intervals: List[Tuple[int, int]]) -> bool:
    for start, end in intervals:
        if start <= pos < end:
            return True
    return False


def compute_open_closed_ratio(pmf_csv: Path, acc_bed: Path) -> float:
    pmf = pd.read_csv(pmf_csv)
    intervals = load_accessibility_intervals(acc_bed)
    open_vals: List[float] = []
    closed_vals: List[float] = []

    for row in pmf.itertuples(index=False):
        chrom = str(row.chrom).split()[0]
        pos = int(row.bin_idx)
        val = float(row.pmf)
        if in_open_region(pos, intervals.get(chrom, [])):
            open_vals.append(val)
        else:
            closed_vals.append(val)

    if not open_vals or not closed_vals:
        return float("nan")
    return float(np.mean(open_vals) / np.mean(closed_vals))


def main() -> None:
    args = parse_args()
    params_csv = args.params_csv or args.results_dir / "params" / "run_params.csv"
    params = pd.read_csv(params_csv, dtype={"run_id": str})

    tables = args.output_dir / "tables"
    figures = args.output_dir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    rows = []
    for row in params.itertuples(index=False):
        run_id = row.run_id
        pmf_csv = args.results_dir / run_id / "con" / "pmf.csv"
        if not pmf_csv.exists():
            continue
        pmf = pd.read_csv(pmf_csv)
        q90 = float(pmf["pmf"].quantile(0.9))
        q10 = float(pmf["pmf"].quantile(0.1))
        gc_shape = q90 / q10 if q10 > 0 else float("inf")
        rows.append(
            {
                "run_id": run_id,
                "gc_exp": row.gc_exp,
                "acc_exp": row.acc_exp,
                "acc_weight": getattr(row, "acc_weight", 1.0),
                "gc_q90_q10": gc_shape,
                "open_closed_ratio": compute_open_closed_ratio(pmf_csv, args.accessibility_bed),
            }
        )

    if not rows:
        raise RuntimeError("No calibration rows were found. Check results directory.")

    per_run = pd.DataFrame(rows)
    summary = (
        per_run.groupby(["gc_exp", "acc_exp", "acc_weight"], as_index=False)
        .agg(
            gc_q90_q10=("gc_q90_q10", "mean"),
            open_closed_ratio=("open_closed_ratio", "mean"),
            n_runs=("run_id", "count"),
        )
        .sort_values(["gc_exp", "acc_exp", "acc_weight"])
    )
    summary["score"] = (
        (summary["gc_q90_q10"] - args.target_gc_q90_q10).abs()
        + (summary["open_closed_ratio"] - args.target_open_closed).abs()
    )

    per_run.to_csv(tables / "calibration_per_run.csv", index=False)
    summary.to_csv(tables / "calibration_summary.csv", index=False)
    summary.sort_values("score").head(1).to_csv(
        tables / "calibration_recommendation.csv", index=False
    )

    fig, ax = plt.subplots(figsize=(7, 5))
    sc = ax.scatter(
        summary["gc_q90_q10"],
        summary["open_closed_ratio"],
        c=summary["score"],
        cmap="viridis_r",
        s=60,
    )
    for row in summary.itertuples(index=False):
        ax.text(row.gc_q90_q10, row.open_closed_ratio, f"g{row.gc_exp}/a{row.acc_exp}", fontsize=7)
    ax.axvline(args.target_gc_q90_q10, linestyle="--", linewidth=1)
    ax.axhline(args.target_open_closed, linestyle="--", linewidth=1)
    ax.set_xlabel("GC PMF q90/q10")
    ax.set_ylabel("Open/Closed PMF mean ratio")
    ax.set_title("Bias calibration candidates")
    fig.colorbar(sc, label="distance-to-target score")
    fig.tight_layout()
    fig.savefig(figures / "calibration_tradeoff.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
