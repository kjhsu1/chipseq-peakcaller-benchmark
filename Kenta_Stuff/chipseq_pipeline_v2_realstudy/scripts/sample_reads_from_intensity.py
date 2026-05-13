"""Label planned coverage sweeps from reference intensity models."""

"""Imports"""

import argparse
from pathlib import Path

import pandas as pd

from scripts.realstudy_sampling_lib import sampling_label


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Create prototype sampling run tables.")
    parser.add_argument("--observed-treatment-depth", type=float, required=True)
    parser.add_argument("--coverage-treat", nargs="+", type=float, required=True)
    parser.add_argument("--coverage-ctrl", nargs="+", type=float, required=True)
    parser.add_argument("--seed", nargs="+", type=int, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Write the prototype run table for coverage sweeps."""
    args = parse_args()
    rows = []
    for cov_treat in args.coverage_treat:
        for cov_ctrl in args.coverage_ctrl:
            for seed in args.seed:
                rows.append(
                    {
                        "coverage_treat": cov_treat,
                        "coverage_ctrl": cov_ctrl,
                        "seed": seed,
                        "sampling_label": sampling_label(cov_treat, args.observed_treatment_depth),
                    }
                )
    df = pd.DataFrame(rows)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_csv, index=False)


main()
