"""Label planned coverage sweeps from reference intensity models."""

"""Imports"""

import argparse
import json
from pathlib import Path

import pandas as pd

from scripts.realstudy_sampling_lib import build_run_table_rows, load_observed_depths


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Create prototype sampling run tables.")
    parser.add_argument("--study-ids", nargs="+", required=True)
    parser.add_argument("--observed-treatment-depth", type=float, default=None)
    parser.add_argument("--observed-control-depth", type=float, default=1.0)
    parser.add_argument("--study-depths-csv", type=Path, default=None)
    parser.add_argument("--coverage-treat", nargs="+", type=float, required=True)
    parser.add_argument("--coverage-ctrl", nargs="+", type=float, required=True)
    parser.add_argument("--seed", nargs="+", type=int, required=True)
    parser.add_argument("--fragment-length", nargs="+", type=int, default=[150])
    parser.add_argument("--read-length", nargs="+", type=int, default=[38])
    parser.add_argument("--aligners", nargs="+", default=["bowtie2"])
    parser.add_argument("--peakcaller-list", nargs="+", default=["macs2"])
    parser.add_argument("--macs2-mode", nargs="+", default=["narrow"])
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    """Write the prototype run table for coverage sweeps."""
    args = parse_args()
    observed_by_study = load_observed_depths(args.study_depths_csv)
    rows = build_run_table_rows(
        study_ids=args.study_ids,
        observed_by_study=observed_by_study,
        coverage_treat=args.coverage_treat,
        coverage_ctrl=args.coverage_ctrl,
        seeds=args.seed,
        fragment_lengths=args.fragment_length,
        read_lengths=args.read_length,
        aligners=args.aligners,
        peakcallers=args.peakcaller_list,
        macs2_modes=args.macs2_mode,
        observed_treatment_depth=args.observed_treatment_depth,
        observed_control_depth=args.observed_control_depth,
    )
    df = pd.DataFrame(rows)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_csv, index=False)
    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(
            json.dumps(
                {
                    "study_count": len(args.study_ids),
                    "coverage_treat_count": len(args.coverage_treat),
                    "coverage_ctrl_count": len(args.coverage_ctrl),
                    "seed_count": len(args.seed),
                    "run_count": len(df),
                },
                indent=2,
            ),
            encoding="utf-8",
        )


main()
