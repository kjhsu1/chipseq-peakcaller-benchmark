"""Assign ontology classes from regional signal summaries."""

"""Imports"""

import argparse
import json
from pathlib import Path

import pandas as pd

from scripts.ontology_lib import DEFAULT_THRESHOLDS, classify_row


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Classify regions into the realstudy ontology.")
    parser.add_argument("--input-csv", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Write ontology labels for a regional summary table."""
    args = parse_args()
    df = pd.read_csv(args.input_csv)
    classified = pd.DataFrame([classify_row(row) for _, row in df.iterrows()])
    result = pd.concat([df.reset_index(drop=True), classified], axis=1)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output_prefix.with_suffix(".csv"), index=False)
    summary = (
        result.groupby(
            ["background_class", "foreground_class", "failure_mode", "ontology_class"], as_index=False
        )
        .size()
        .rename(columns={"size": "region_count"})
        .sort_values(["background_class", "foreground_class", "ontology_class"])
    )
    summary.to_csv(args.output_prefix.with_name(f"{args.output_prefix.name}_summary.csv"), index=False)
    if {"chrom", "start", "end"}.issubset(result.columns):
        result[["chrom", "start", "end", "ontology_class"]].to_csv(
            args.output_prefix.with_suffix(".bed"),
            index=False,
            sep="\t",
            header=False,
        )
    args.output_prefix.with_name(f"{args.output_prefix.name}_definition.json").write_text(
        json.dumps(DEFAULT_THRESHOLDS, indent=2),
        encoding="utf-8",
    )


main()
