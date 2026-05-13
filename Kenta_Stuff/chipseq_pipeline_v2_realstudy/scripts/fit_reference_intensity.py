"""Fit simple reference-intensity summaries from binned tracks."""

"""Imports"""

import argparse
import json
from pathlib import Path

import pandas as pd


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Fit reference intensity summaries from binned tracks.")
    parser.add_argument("--binned-track-csv", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--clip-percentile", type=float, default=99.9)
    return parser.parse_args()


def main() -> None:
    """Write clipped reference-intensity tables and metadata."""
    args = parse_args()
    df = pd.read_csv(args.binned_track_csv)
    clipped = df.copy()
    for column in ["treatment", "control", "log2_treat_control"]:
        upper = clipped[column].quantile(args.clip_percentile / 100.0)
        clipped[column] = clipped[column].clip(upper=upper)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    clipped.to_csv(args.output_prefix.with_suffix(".csv"), index=False)
    args.output_prefix.with_suffix(".json").write_text(
        json.dumps(
            {
                "clip_percentile": args.clip_percentile,
                "language": "estimated reference intensity, not true distribution",
            },
            indent=2,
        ),
        encoding="utf-8",
    )


main()
