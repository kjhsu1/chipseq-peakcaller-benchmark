"""Prepare resumable real-study ingest bookkeeping without live network fetches."""

"""Imports"""

import argparse
from datetime import datetime
import json
from pathlib import Path

import pandas as pd

from scripts.realstudy_manifest_lib import (
    build_data_manifest,
    load_file_manifest,
    load_study_selection,
    validate_selection_table,
)


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Prepare resumable ingest manifests for real-study data."
    )
    parser.add_argument(
        "--study-selection",
        type=Path,
        default=Path("manifests/study_selection.yaml"),
        help="Input study selection YAML manifest",
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
    parser.add_argument(
        "--metadata-dir",
        type=Path,
        default=Path("metadata"),
        help="Output directory for normalized metadata tables",
    )
    parser.add_argument(
        "--encode-cache-dir",
        type=Path,
        default=Path("metadata"),
        help="Optional directory containing cached ENCODE experiment JSON files",
    )
    return parser.parse_args()


def build_encode_cache_summary(selection_df: pd.DataFrame, cache_dir: Path) -> pd.DataFrame:
    """Summarize which ENCODE studies have cached live metadata on disk."""
    rows = []
    encode_rows = selection_df[selection_df["source"].astype(str).str.upper() == "ENCODE"]
    for row in encode_rows.to_dict(orient="records"):
        accession = str(row.get("accession", "")).strip()
        control_accession = str(row.get("control_accession", "")).strip()
        cache_path = cache_dir / f"{accession}.json" if accession else None
        control_cache_path = cache_dir / f"{control_accession}.json" if control_accession else None
        rows.append(
            {
                "study_id": row["study_id"],
                "accession": accession,
                "control_accession": control_accession,
                "experiment_json_cached": bool(cache_path and cache_path.exists()),
                "control_json_cached": bool(control_cache_path and control_cache_path.exists()),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    """Write a resumable ingest/data manifest and stop before network access."""
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_dir.mkdir(parents=True, exist_ok=True)
    selection_df = validate_selection_table(load_study_selection(args.study_selection))
    file_df = load_file_manifest(args.file_manifest)
    data_manifest = build_data_manifest(selection_df, file_df, local_search_roots=[Path.cwd()])

    selection_df.to_csv(args.metadata_dir / "study_selection.normalized.csv", index=False)
    data_manifest.to_csv(args.metadata_dir / "data_manifest.csv", index=False)
    data_manifest.to_csv(args.output_dir / "ingest_plan.csv", index=False)
    summary = (
        data_manifest.groupby(["study_id", "metadata_status", "ingest_mode"], as_index=False)
        .size()
        .rename(columns={"size": "file_count"})
        .sort_values(["study_id", "metadata_status", "ingest_mode"])
    )
    summary.to_csv(args.output_dir / "ingest_summary.csv", index=False)
    download_ready = data_manifest[data_manifest["ingest_mode"].isin(["download_processed_bam", "download_fastq"])].copy()
    if not download_ready.empty:
        download_ready.to_csv(args.output_dir / "download_ready_manifest.csv", index=False)
    encode_cache_summary = build_encode_cache_summary(selection_df, args.encode_cache_dir)
    if not encode_cache_summary.empty:
        encode_cache_summary.to_csv(args.output_dir / "encode_cache_summary.csv", index=False)
    (args.output_dir / "BLOCKED_NETWORK_STEP.txt").write_text(
        "Next blocked step: query live ENCODE or NCBI endpoints for metadata and downloads.\n"
        "Exact capability needed: plain network access for curl/wget/API requests.\n"
        "Plugin-specific access is optional, not required, unless you want connector-mediated provenance.\n",
        encoding="utf-8",
    )


main()
