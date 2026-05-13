"""Assign simple ontology classes from regional signal summaries."""

"""Imports"""

import argparse
from pathlib import Path

import pandas as pd

from scripts.ontology_lib import classify_row


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Classify regions into a simple ontology.")
    parser.add_argument("--input-csv", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Write ontology labels for a regional summary table."""
    args = parse_args()
    df = pd.read_csv(args.input_csv)
    df["ontology_class"] = df.apply(classify_row, axis=1)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_csv, index=False)


main()
