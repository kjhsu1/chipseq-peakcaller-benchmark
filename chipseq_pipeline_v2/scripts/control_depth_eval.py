"""Aggregate archived sweep results by category and aligner/peakcaller combo."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts.eval_helpers import counts_to_metrics, resolve_peak_path


class OverlapCounts(NamedTuple):
    """Plain immutable container for overlap count totals."""

    tp_called: int
    total_called: int
    tp_planted: int
    total_planted: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Aggregate archived sweep results")
    parser.add_argument("--input-dirs", nargs="+", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--group-labels", nargs="*", default=[])
    return parser.parse_args()


def resolve_groups(input_dirs: List[Path], labels: List[str]) -> Dict[Path, Optional[str]]:
    if labels:
        if len(labels) != len(input_dirs):
            raise ValueError("--group-labels must match --input-dirs length")
        return {p: labels[i] for i, p in enumerate(input_dirs)}
    return {p: None for p in input_dirs}


def classify_category(row: pd.Series) -> Optional[str]:
    use_control = bool(row.get("use_control", False))
    tf_count = float(row.get("tf_peak_count_treat", 0))
    gc_exp = float(row.get("gc_exp", 0))
    acc_exp = float(row.get("acc_exp", 0))
    map_coverage_pct = float(row.get("map_coverage_pct", 0))
    tf_sigma = float(row.get("tf_sigma", 0))
    macs2_mode = str(row.get("macs2_mode", "narrow"))

    if tf_count == 0 and not use_control and gc_exp == 0 and acc_exp == 0 and map_coverage_pct == 0:
        if macs2_mode == "narrow":
            return "shotgun_flat_narrow_none"
        if macs2_mode == "broad":
            return "shotgun_flat_broad_none"
        return None

    if tf_count <= 0 or not use_control:
        return None

    # Primary benchmark boundary:
    # - TF-like / narrow: tf_sigma <= 5
    # - Histone-like / broad: tf_sigma >= 15
    # - 5 < tf_sigma < 15 is ambiguous and excluded from headline categories
    if gc_exp == 0 and acc_exp == 0 and map_coverage_pct == 0 and macs2_mode == "broad":
        if tf_sigma <= 5:
            return "flatearth_peaks_broad_integrated"
        if tf_sigma >= 15:
            return "flatearth_plateaus_broad_integrated"
        return None

    if gc_exp > 0 and acc_exp > 0:
        background = "wavy" if map_coverage_pct == 0 else "hilly"
        if macs2_mode == "narrow" and tf_sigma <= 5:
            return f"realistic_peaks_{background}_narrow_integrated"
        if macs2_mode == "broad" and tf_sigma >= 15:
            return f"realistic_plateaus_{background}_broad_integrated"

    return None


def infer_method_combo(row: pd.Series) -> str:
    aligner = str(row.get("aligner", "unknown"))
    peakcaller = str(row.get("peakcaller", "unknown"))
    return f"{aligner}+{peakcaller}"


def truth_path_for_row(results_dir: Path, row: pd.Series) -> Path:
    """Return the appropriate planted-truth file for one run."""
    if str(row.get("peakcaller", "macs2")) == "epic2" or str(row.get("macs2_mode", "narrow")) == "broad":
        return results_dir / row["run_id"] / "treat" / "planted_peak_intervals.bed"
    return results_dir / row["run_id"] / "treat" / "planted_peaks.bed"


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


def compute_summary(per_run: pd.DataFrame, group_cols: List[str]) -> pd.DataFrame:
    summary = (
        per_run.groupby(group_cols + ["ratio"], as_index=False)
        .agg(
            tp_called=("tp_called", "sum"),
            total_called=("total_called", "sum"),
            tp_planted=("tp_planted", "sum"),
            total_planted=("total_planted", "sum"),
            n_runs=("run_id", "count"),
            run_ids=("run_id", lambda v: ",".join(sorted(v))),
            source_dirs=("source_dir", lambda v: ",".join(sorted(set(v)))),
        )
        .sort_values(group_cols + ["ratio"])
    )

    metrics = summary.apply(
        lambda row: counts_to_metrics(
            row["tp_called"],
            row["total_called"],
            row["tp_planted"],
            row["total_planted"],
        ),
        axis=1,
        result_type="expand",
    )
    metrics.columns = ["precision", "recall", "f1"]
    summary[["precision", "recall", "f1"]] = metrics
    summary["fdr"] = 1.0 - summary["precision"]

    summary["inflation"] = np.nan
    for keys, grp in summary.groupby(group_cols):
        baseline_ratio = grp["ratio"].max()
        baseline_called = float(grp.loc[grp["ratio"] == baseline_ratio, "total_called"].iloc[0])
        if not isinstance(keys, tuple):
            keys = (keys,)
        mask = np.ones(len(summary), dtype=bool)
        for col, key in zip(group_cols, keys):
            mask &= summary[col] == key
        if baseline_called > 0:
            summary.loc[mask, "inflation"] = summary.loc[mask, "total_called"] / baseline_called

    return summary


def compute_aggregate_table(per_run: pd.DataFrame, group_cols: List[str]) -> pd.DataFrame:
    summary = (
        per_run.groupby(group_cols, as_index=False)
        .agg(
            tp_called=("tp_called", "sum"),
            total_called=("total_called", "sum"),
            tp_planted=("tp_planted", "sum"),
            total_planted=("total_planted", "sum"),
            n_runs=("run_id", "count"),
        )
        .sort_values(group_cols)
    )

    metrics = summary.apply(
        lambda row: counts_to_metrics(
            row["tp_called"],
            row["total_called"],
            row["tp_planted"],
            row["total_planted"],
        ),
        axis=1,
        result_type="expand",
    )
    metrics.columns = ["precision", "recall", "f1"]
    summary[["precision", "recall", "f1"]] = metrics
    summary["fdr"] = 1.0 - summary["precision"]

    if "method_combo" in group_cols:
        baseline = summary.groupby("method_combo", as_index=False)["total_called"].max().rename(
            columns={"total_called": "baseline_called"}
        )
        summary = summary.merge(baseline, on="method_combo", how="left")
    else:
        summary["baseline_called"] = float(summary["total_called"].max())

    summary["inflation"] = np.where(summary["baseline_called"] > 0, summary["total_called"] / summary["baseline_called"], np.nan)
    return summary.drop(columns=["baseline_called"])


def write_outputs(summary: pd.DataFrame, output_dir: Path) -> None:
    tables_dir = output_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    per_run = summary.attrs.get("per_run")
    if per_run is not None:
        per_run.sort_values(["category", "method_combo", "ratio", "run_id"]).to_csv(
            tables_dir / "per_run_stats.csv", index=False
        )

    summary[
        ["category", "method_combo", "ratio", "precision", "recall", "f1", "fdr", "inflation", "n_runs"]
    ].to_csv(tables_dir / "category_method_ratio_summary.csv", index=False)

    summary[["category", "method_combo", "ratio", "n_runs", "run_ids", "source_dirs"]].to_csv(
        tables_dir / "figure_table_manifest.csv", index=False
    )

    category_method_summary = compute_aggregate_table(per_run, ["category", "method_combo"])
    category_method_summary.to_csv(tables_dir / "category_method_summary.csv", index=False)

    category_summary = compute_aggregate_table(per_run, ["category"])
    category_summary.to_csv(tables_dir / "category_summary.csv", index=False)

    for (category, method_combo), grp in summary.groupby(["category", "method_combo"]):
        grp = grp.sort_values("ratio")
        fig_dir = output_dir / category / method_combo.replace("/", "_") / "figures"
        fig_dir.mkdir(parents=True, exist_ok=True)
        ratios = grp["ratio"].to_numpy()

        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot(ratios, grp["precision"], marker="o", label="precision")
        ax.plot(ratios, grp["recall"], marker="o", label="recall")
        ax.plot(ratios, grp["f1"], marker="o", label="f1")
        ax.set_xlabel("control:treat ratio")
        ax.set_ylabel("metric value")
        ax.set_ylim(0, 1.05)
        ax.set_title(f"PR/F1 vs ratio ({category}, {method_combo})")
        ax.legend()
        fig.tight_layout()
        fig.savefig(fig_dir / "pr_f1_vs_ratio.png", dpi=200)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(5.5, 5.5))
        ax.plot(grp["recall"], grp["precision"], marker="o")
        for _, row in grp.iterrows():
            ax.annotate(f"{row['ratio']:.3f}", (row["recall"], row["precision"]), fontsize=7, xytext=(4, 4), textcoords="offset points")
        ax.set_xlabel("recall")
        ax.set_ylabel("precision")
        ax.set_xlim(0, 1.05)
        ax.set_ylim(0, 1.05)
        ax.set_title(f"Precision-Recall ({category}, {method_combo})")
        fig.tight_layout()
        fig.savefig(fig_dir / "precision_recall_curve.png", dpi=200)
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
        ax1.set_title(f"FDR/Inflation vs ratio ({category}, {method_combo})")
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
        ax.set_title(f"Interaction heatmap ({category}, {method_combo})")
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
        for row in params.itertuples(index=False):
            run_id = row.run_id
            row_series = pd.Series(row._asdict())
            planted_path = truth_path_for_row(input_dir, row_series)
            peak_path = resolve_peak_path(input_dir, run_id, getattr(row, "peakcaller", "macs2"))
            if not planted_path.exists() or not peak_path.exists():
                continue

            category = group_map[input_dir] or classify_category(row_series)
            if category is None:
                continue

            counts = compute_overlap_stats(
                load_planted_centers(planted_path),
                load_called_intervals(peak_path),
            )

            per_run_rows.append(
                {
                    "category": category,
                    "method_combo": infer_method_combo(row_series),
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

    per_run_df = pd.DataFrame(per_run_rows)
    summary = compute_summary(per_run_df, ["category", "method_combo"])
    summary.attrs["per_run"] = per_run_df
    write_outputs(summary, args.output_dir)


if __name__ == "__main__":
    main()
