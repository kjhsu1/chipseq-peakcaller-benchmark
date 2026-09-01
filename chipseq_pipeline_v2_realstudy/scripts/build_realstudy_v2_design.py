"""Build and validate the locked Realstudy v2 run table."""

"""Imports"""

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
import yaml

from realstudy_sampling_lib import audit_realstudy_v2_design, build_realstudy_v2_run_rows


"""Constants"""

REQUIRED_FILE_COLUMNS = {
    "file_id", "study_id", "experiment_id", "role", "fastq_url", "raw_reads",
    "compressed_bytes", "md5", "expected_local_path",
}


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse design-builder arguments."""
    parser = argparse.ArgumentParser(description="Build the Realstudy v2 44-call design.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--files-manifest", type=Path, required=True)
    parser.add_argument("--run-table", type=Path, required=True)
    parser.add_argument("--audit-json", type=Path, required=True)
    return parser.parse_args()


def validate_file_registry(path: Path, design: dict) -> dict:
    """Validate input accession, checksum, count, and study coverage."""
    frame = pd.read_csv(path, dtype=str).fillna("")
    missing = sorted(REQUIRED_FILE_COLUMNS - set(frame.columns))
    if missing:
        raise ValueError(f"File registry is missing columns: {', '.join(missing)}")
    if frame["file_id"].duplicated().any() or frame["file_id"].str.strip().eq("").any():
        raise ValueError("File identifiers must be unique and nonempty.")
    for column in ["raw_reads", "compressed_bytes"]:
        values = pd.to_numeric(frame[column], errors="coerce")
        if values.isna().any() or values.le(0).any():
            raise ValueError(f"{column} must contain positive integers.")
    if not frame["md5"].str.fullmatch(r"[0-9a-f]{32}").all():
        raise ValueError("Every input file must have a lowercase 32-character MD5.")
    required_files = {
        file_id
        for study in design["studies"].values()
        for role in ["treatment_files", "control_files"]
        for file_id in study[role]
    }
    observed_files = set(frame["file_id"])
    if required_files != observed_files:
        raise ValueError(
            f"Configured and registered files differ: missing={sorted(required_files - observed_files)}, "
            f"extra={sorted(observed_files - required_files)}"
        )
    return {"input_files": len(frame), "registry_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def main() -> None:
    """Write the deterministic run table and static design audit."""
    args = parse_args()
    loaded = yaml.safe_load(args.config.read_text())
    design = loaded["realstudy_v2"]
    rows = build_realstudy_v2_run_rows(design)
    audit = audit_realstudy_v2_design(rows, design)
    audit.update(validate_file_registry(args.files_manifest, design))
    args.run_table.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.run_table, index=False)
    args.audit_json.parent.mkdir(parents=True, exist_ok=True)
    args.audit_json.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")


main()
