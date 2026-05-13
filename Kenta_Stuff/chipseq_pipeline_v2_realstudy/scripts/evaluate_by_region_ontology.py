"""Aggregate ontology-based evaluation tables."""

"""Imports"""

import argparse
from pathlib import Path

import pandas as pd


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Evaluate performance by ontology class.")
    parser.add_argument("--input-csv", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Aggregate performance metrics by ontology class and control depth."""
    args = parse_args()
    df = pd.read_csv(args.input_csv)
    summary = (
        df.groupby(["ontology_class", "coverage_ctrl"], as_index=False)[["precision", "recall", "f1"]]
        .mean()
        .sort_values(["ontology_class", "coverage_ctrl"])
    )
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output_csv, index=False)


main()
