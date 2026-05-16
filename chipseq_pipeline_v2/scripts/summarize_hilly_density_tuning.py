"""Summarize hilly-density tuning candidates from existing tuning outputs."""

"""Imports"""

import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd

from scripts.hilly_tuning_lib import build_recommendation_table


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Summarize mappability-density tuning candidates for hilly backgrounds."
    )
    parser.add_argument(
        "--map-qc-csv",
        type=Path,
        default=Path("analysis_outputs/tfclean_map_tuning_20260502_202816/map_bias_depth_qc.csv"),
        help="Existing map tuning QC table",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis_outputs") / f"hilly_density_tuning_{datetime.now().strftime('%Y%m%d')}",
        help="Output directory for the ranked tuning summary",
    )
    return parser.parse_args()


def main() -> None:
    """Build a ranked recommendation table for hilly-density tuning."""
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if not args.map_qc_csv.exists():
        lines = [
            "# Hilly Density Tuning",
            "",
            f"- source QC table missing: `{args.map_qc_csv}`",
            "- rerun this summary after the lite pilot QC table is available",
        ]
        (args.output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
        pd.DataFrame(
            [
                {
                    "recommendation": "blocked_missing_qc",
                    "reason": "source QC table not found",
                    "map_qc_csv": str(args.map_qc_csv),
                }
            ]
        ).to_csv(args.output_dir / "ranked_recommendations.csv", index=False)
        return
    ranked = build_recommendation_table(pd.read_csv(args.map_qc_csv))
    ranked.to_csv(args.output_dir / "ranked_recommendations.csv", index=False)

    lines = [
        "# Hilly Density Tuning",
        "",
        f"- source QC table: `{args.map_qc_csv}`",
        "- ranking rejects candidates with excessive false-positive burden versus wavy when that column is available",
        "- preferred settings are the mildest candidates that remain measurably distinct from wavy in bumpiness proxies",
    ]
    (args.output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


main()
