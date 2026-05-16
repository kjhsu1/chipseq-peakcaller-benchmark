"""Audit repeated control-depth outcomes in per-run evaluation outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit repeated high-control outcomes")
    parser.add_argument("--per-run-stats", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = pd.read_csv(args.per_run_stats, dtype={"run_id": str})
    if "control_mapped_reads" not in data.columns:
        data["control_mapped_reads"] = pd.NA
    if "treat_mapped_reads" not in data.columns:
        data["treat_mapped_reads"] = pd.NA
    rows = []
    for keys, group in data.groupby(["coverage_treat", "tf_enrich", "seed"]):
        ordered = group.sort_values("coverage_ctrl")
        high = ordered[ordered["coverage_ctrl"] >= 4].copy()
        if high.empty:
            continue
        unique_metric_states = {
            (
                row["coverage_ctrl"],
                row["tp_called"],
                row["total_called"],
                row["tp_planted"],
                row.get("control_mapped_reads"),
                row.get("treat_mapped_reads"),
            )
            for _, row in high.iterrows()
        }
        rows.append(
            {
                "coverage_treat": keys[0],
                "tf_enrich": keys[1],
                "seed": keys[2],
                "n_high_control_rows": len(high),
                "n_unique_high_control_states": len(unique_metric_states),
                "min_control_coverage": float(high["coverage_ctrl"].min()),
                "max_control_coverage": float(high["coverage_ctrl"].max()),
                "control_mapped_reads_unique": int(high["control_mapped_reads"].nunique(dropna=True)),
                "total_called_unique": int(high["total_called"].nunique(dropna=True)),
            }
        )
    audit = pd.DataFrame(rows).sort_values(["coverage_treat", "tf_enrich", "seed"])
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(args.output_csv, index=False)


if __name__ == "__main__":
    main()
