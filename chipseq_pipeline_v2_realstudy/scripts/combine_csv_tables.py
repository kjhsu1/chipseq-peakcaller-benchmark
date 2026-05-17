"""Combine CSV tables with matching columns."""

"""Imports"""

import argparse
from pathlib import Path

import pandas as pd


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Combine CSV tables into a single output file.")
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("input_csvs", nargs="+", type=Path)
    return parser.parse_args()


def main() -> None:
    """Concatenate input CSV tables."""
    args = parse_args()
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    tables = [pd.read_csv(path) for path in args.input_csvs]
    pd.concat(tables, ignore_index=True).to_csv(args.output_csv, index=False)


main()
