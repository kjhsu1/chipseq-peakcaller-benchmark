"""Fit smoothed reference-intensity summaries from binned tracks."""

"""Imports"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Fit reference intensity summaries from binned tracks.")
    parser.add_argument("--binned-track-csv", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--study-id", type=str, default="unknown_study")
    parser.add_argument("--clip-percentile", type=float, default=99.9)
    parser.add_argument("--smoothing-window", type=int, default=5)
    return parser.parse_args()


def smooth_by_contig(df: pd.DataFrame, column: str, window: int) -> pd.Series:
    """Apply centered rolling-mean smoothing within each contig."""
    return (
        df.groupby("chrom", group_keys=False)[column]
        .apply(lambda series: series.rolling(window=window, min_periods=1, center=True).mean())
    )


def main() -> None:
    """Write clipped reference-intensity tables and metadata."""
    args = parse_args()
    df = pd.read_csv(args.binned_track_csv)
    clipped = df.copy()
    for column in ["treatment", "control", "log2_treat_control"]:
        upper = clipped[column].quantile(args.clip_percentile / 100.0)
        clipped[column] = clipped[column].clip(upper=upper)
        clipped[f"reference_{column}"] = smooth_by_contig(clipped, column, args.smoothing_window)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    clipped.to_csv(args.output_prefix.with_suffix(".csv"), index=False)
    for column in ["treatment", "control", "log2_treat_control"]:
        clipped[["study_id", "chrom", "start", "end", f"reference_{column}"]].rename(
            columns={f"reference_{column}": "value"}
        ).to_csv(args.output_prefix.with_suffix(f".{column}.csv"), index=False)
    args.output_prefix.with_suffix(".json").write_text(
        json.dumps(
            {
                "study_id": args.study_id,
                "clip_percentile": args.clip_percentile,
                "smoothing": "rolling_mean",
                "smoothing_window": args.smoothing_window,
                "row_count": int(len(clipped)),
                "language": "estimated reference intensity, not true distribution",
            },
            indent=2,
        ),
        encoding="utf-8",
    )


main()
