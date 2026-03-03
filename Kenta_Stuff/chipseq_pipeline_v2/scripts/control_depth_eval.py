"""Cumulative group-level evaluation for control-depth sweeps."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


KNOWN_GROUPS = [
    "proxy_uniform",
    "proxy_bumpy",
    "ctrltreat_uniform",
    "ctrltreat_bumpy",
]


@dataclass(frozen=True)
class OverlapCounts:
    tp_called: int
    total_called: int
    tp_planted: int
    total_planted: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Aggregate control-depth sweep results")
    parser.add_argument("--input-dirs", nargs="+", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--group-labels", nargs="*", default=[])
    return parser.parse_args()


def parse_group_from_dir(path: Path) -> str:
    name = path.name
    for group in KNOWN_GROUPS:
        if group in name:
            return group
    return name


def resolve_groups(input_dirs: List[Path], labels: List[str]) -> Dict[Path, str]:
    if labels:
        if len(labels) != len(input_dirs):
            raise ValueError("--group-labels must match --input-dirs length")
        return {p: labels[i] for i, p in enumerate(input_dirs)}
    return {p: parse_group_from_dir(p) for p in input_dirs}


def resolve_peak_path(results_dir: Path, run_id: str, peakcaller: str) -> Path:
    if peakcaller == "epic2":
        return results_dir / run_id / "peaks" / "epic2" / f"{run_id}_domains.bed"
    normalized = results_dir / run_id / "peaks" / "macs2" / f"{run_id}_peaks.bed"
    if normalized.exists():
        return normalized
    return results_dir / run_id / "peaks" / "macs2" / f"{run_id}_peaks.narrowPeak"


def load_planted_centers(path: Path) -> Dict[str, set[int]]:
    centers: Dict[str, set[int]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, *_ = line.rstrip().split("\t")
            centers.setdefault(chrom, set()).add(int(start))
    return centers


def load_called_intervals(path: Path) -> Dict[str, List[Tuple[int, int]]]:
    intervals: Dict[str, List[Tuple[int, int]]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *_ = line.rstrip().split("\t")
            intervals.setdefault(chrom, []).append((int(start), int(end)))
    return intervals


def interval_overlaps_any(center: int, intervals: List[Tuple[int, int]]) -> bool:
    for start, end in intervals:
        if start <= center < end:
            return True
    return False


def compute_overlap_stats(
    planted: Dict[str, set[int]],
    called: Dict[str, List[Tuple[int, int]]],
) -> OverlapCounts:
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

    return OverlapCounts(tp_called, total_called, tp_planted, total_planted)


def compute_summary(per_run: pd.DataFrame) -> pd.DataFrame:
    summary = (
        per_run.groupby(["group", "ratio"], as_index=False)
        .agg(
            tp_called=("tp_called", "sum"),
            total_called=("total_called", "sum"),
            tp_planted=("tp_planted", "sum"),
            total_planted=("total_planted", "sum"),
            n_runs=("run_id", "count"),
            run_ids=("run_id", lambda v: ",".join(sorted(v))),
            source_dirs=("source_dir", lambda v: ",".join(sorted(set(v)))),
        )
        .sort_values(["group", "ratio"])
    )

    summary["precision"] = np.where(summary["total_called"] > 0, summary["tp_called"] / summary["total_called"], 0.0)
    summary["recall"] = np.where(summary["total_planted"] > 0, summary["tp_planted"] / summary["total_planted"], 0.0)
    denom = summary["precision"] + summary["recall"]
    summary["f1"] = np.where(denom > 0, 2 * (summary["precision"] * summary["recall"]) / denom, 0.0)
    summary["fdr"] = 1.0 - summary["precision"]

    summary["inflation"] = np.nan
    for group, grp in summary.groupby("group"):
        baseline_ratio = grp["ratio"].max()
        baseline_called = float(grp.loc[grp["ratio"] == baseline_ratio, "total_called"].iloc[0])
        mask = summary["group"] == group
        if baseline_called > 0:
            summary.loc[mask, "inflation"] = summary.loc[mask, "total_called"] / baseline_called

    return summary


def write_outputs(summary: pd.DataFrame, output_dir: Path) -> None:
    tables_dir = output_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    summary[
        ["group", "ratio", "precision", "recall", "f1", "fdr", "inflation", "n_runs"]
    ].to_csv(tables_dir / "group_ratio_summary.csv", index=False)

    summary[["group", "ratio", "n_runs", "run_ids", "source_dirs"]].to_csv(
        tables_dir / "figure_table_manifest.csv", index=False
    )

    for group, grp in summary.groupby("group"):
        grp = grp.sort_values("ratio")
        fig_dir = output_dir / f"group_{group}" / "figures"
        fig_dir.mkdir(parents=True, exist_ok=True)
        ratios = grp["ratio"].to_numpy()

        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot(ratios, grp["precision"], marker="o", label="precision")
        ax.plot(ratios, grp["recall"], marker="o", label="recall")
        ax.plot(ratios, grp["f1"], marker="o", label="f1")
        ax.set_xlabel("control:treat ratio")
        ax.set_ylabel("metric value")
        ax.set_ylim(0, 1.05)
        ax.set_title(f"PR/F1 vs ratio ({group})")
        ax.legend()
        fig.tight_layout()
        fig.savefig(fig_dir / "pr_f1_vs_ratio.png", dpi=200)
        plt.close(fig)

        fig, ax1 = plt.subplots(figsize=(7, 4.5))
        l1 = ax1.plot(ratios, grp["fdr"], marker="o", color="tab:red", label="fdr")
        ax1.set_xlabel("control:treat ratio")
        ax1.set_ylabel("fdr", color="tab:red")
        ax1.tick_params(axis="y", labelcolor="tab:red")
        ax1.set_ylim(0, 1.05)
        ax2 = ax1.twinx()
        l2 = ax2.plot(ratios, grp["inflation"], marker="s", color="tab:blue", label="inflation")
        ax2.set_ylabel("peak inflation", color="tab:blue")
        ax2.tick_params(axis="y", labelcolor="tab:blue")
        ax2.set_ylim(bottom=0)
        lines = l1 + l2
        labels = [l.get_label() for l in lines]
        ax1.legend(lines, labels, loc="best")
        ax1.set_title(f"FDR/Inflation vs ratio ({group})")
        fig.tight_layout()
        fig.savefig(fig_dir / "fdr_inflation_vs_ratio.png", dpi=200)
        plt.close(fig)

        metric_names = ["precision", "recall", "f1", "fdr", "inflation"]
        heat = np.array([grp[m].to_numpy() for m in metric_names])

        fig, ax = plt.subplots(figsize=(7, 3.8))
        im = ax.imshow(heat, aspect="auto", cmap="magma")
        ax.set_xticks(range(len(ratios)))
        ax.set_xticklabels([f"{r:.3f}" for r in ratios])
        ax.set_yticks(range(len(metric_names)))
        ax.set_yticklabels(metric_names)
        ax.set_xlabel("control:treat ratio")
        ax.set_title(f"Interaction heatmap ({group})")
        fig.colorbar(im, ax=ax, label="value")
        fig.tight_layout()
        fig.savefig(fig_dir / "interaction_heatmap.png", dpi=200)
        plt.close(fig)


def main() -> None:
    args = parse_args()
    group_map = resolve_groups(args.input_dirs, args.group_labels)

    per_run_rows: List[Dict] = []
    for input_dir in args.input_dirs:
        params_csv = input_dir / "params" / "run_params.csv"
        params = pd.read_csv(params_csv, dtype={"run_id": str})
        group = group_map[input_dir]

        for row in params.itertuples(index=False):
            run_id = row.run_id
            planted_path = input_dir / run_id / "treat" / "planted_peaks.bed"
            peak_path = resolve_peak_path(input_dir, run_id, getattr(row, "peakcaller", "macs2"))
            if not planted_path.exists() or not peak_path.exists():
                continue

            counts = compute_overlap_stats(
                load_planted_centers(planted_path),
                load_called_intervals(peak_path),
            )

            per_run_rows.append(
                {
                    "group": group,
                    "source_dir": str(input_dir),
                    "run_id": run_id,
                    "ratio": float(row.coverage_ctrl) / float(row.coverage_treat),
                    "tp_called": counts.tp_called,
                    "total_called": counts.total_called,
                    "tp_planted": counts.tp_planted,
                    "total_planted": counts.total_planted,
                }
            )

    if not per_run_rows:
        raise RuntimeError("No aggregate rows were found. Check input directories and outputs.")

    write_outputs(compute_summary(pd.DataFrame(per_run_rows)), args.output_dir)


if __name__ == "__main__":
    main()
