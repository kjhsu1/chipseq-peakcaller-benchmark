"""Register real-study ingest tasks without executing live external fetches."""

"""Imports"""

import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Prepare resumable ingest manifests for real-study data."
    )
    parser.add_argument(
        "--file-manifest",
        type=Path,
        default=Path("manifests/study_file_manifest.csv"),
        help="Input study file manifest",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis_outputs") / f"realstudy_ingest_prep_{datetime.now().strftime('%Y%m%d')}",
        help="Output directory for ingest planning tables",
    )
    return parser.parse_args()


def main() -> None:
    """Write a resumable ingest plan and stop before network access."""
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.file_manifest)
    manifest["download_status"] = "not_started"
    manifest["ingest_mode"] = manifest["file_format"].apply(
        lambda value: "download_processed_bam" if str(value) == "bam" else "resolve_metadata_first"
    )
    manifest.to_csv(args.output_dir / "ingest_plan.csv", index=False)
    (args.output_dir / "BLOCKED_NETWORK_STEP.txt").write_text(
        "Next blocked step: query live ENCODE or NCBI endpoints for metadata and downloads.\n"
        "Exact capability needed: plain network access for curl/wget/API requests.\n"
        "Plugin-specific access is optional, not required, unless you want connector-mediated provenance.\n",
        encoding="utf-8",
    )


main()
