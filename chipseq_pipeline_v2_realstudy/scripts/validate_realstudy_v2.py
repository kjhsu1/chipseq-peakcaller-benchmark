"""Perform final integrity checks for a completed Realstudy v2 output package."""

"""Imports"""

import argparse
import hashlib
import json
import sqlite3
from pathlib import Path

import pandas as pd


"""Constants"""

REQUIRED_TABLES = {
    "schema_metadata", "publications", "publication_claims", "studies", "experiments",
    "samples", "input_files", "experiment_publications", "parent_libraries", "software_versions",
    "sampling_blocks", "control_subsamples", "runs", "artifacts", "validation_events", "peaks",
    "reference_regions", "peak_overlaps", "reference_region_recovery", "run_metrics",
    "seed_pair_metrics", "stratified_metrics", "enough_control_decisions",
}


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse integrity-validator arguments."""
    parser = argparse.ArgumentParser(description="Validate a completed Realstudy v2 package.")
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--run-table", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--export-manifest", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    return parser.parse_args()


def checksum(path: Path) -> str:
    """Return a streaming SHA-256 checksum."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    """Check database constraints, run completeness, exports, and all figure assets."""
    args = parse_args()
    checks = []
    connection = sqlite3.connect(f"file:{args.database}?mode=ro", uri=True)
    try:
        tables = {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        checks.append({"check": "required_database_tables", "status": "pass" if REQUIRED_TABLES <= tables else "fail", "detail": sorted(REQUIRED_TABLES - tables)})
        violations = connection.execute("PRAGMA foreign_key_check").fetchall()
        checks.append({"check": "foreign_keys", "status": "pass" if not violations else "fail", "detail": violations})
        run_counts = dict(connection.execute("SELECT run_type, COUNT(*) FROM runs GROUP BY run_type"))
        checks.append({"check": "44_call_matrix", "status": "pass" if run_counts == {"control_subsample": 42, "full_control_anchor": 2} else "fail", "detail": run_counts})
        decisions = connection.execute("SELECT COUNT(*) FROM enough_control_decisions").fetchone()[0]
        checks.append({"check": "study_decisions", "status": "pass" if decisions == 2 else "fail", "detail": decisions})
    finally:
        connection.close()
    run_table = pd.read_csv(args.run_table)
    checks.append({"check": "run_table_rows", "status": "pass" if len(run_table) == 44 else "fail", "detail": len(run_table)})
    exports = pd.read_csv(args.export_manifest)
    bad_exports = []
    for row in exports.itertuples():
        path = Path(row.csv_path)
        if not path.exists() or checksum(path) != row.sha256:
            bad_exports.append(str(path))
    checks.append({"check": "csv_export_checksums", "status": "pass" if not bad_exports else "fail", "detail": bad_exports})
    missing_figures = []
    for number in range(1, 7):
        candidates = list(args.figure_dir.glob(f"figure_{number}_*/figure_{number}_*"))
        extensions = {path.suffix for path in candidates}
        required = {".pdf", ".svg", ".png", ".csv", ".md"}
        if not required <= extensions:
            missing_figures.append({"figure": number, "missing_extensions": sorted(required - extensions)})
    checks.append({"check": "paper_figure_package", "status": "pass" if not missing_figures else "fail", "detail": missing_figures})
    missing_supplements = []
    for number in range(1, 4):
        candidates = list((args.figure_dir / "supplementary").glob(f"supplementary_{number}_*/supplementary_{number}_*"))
        extensions = {path.suffix for path in candidates}
        required = {".pdf", ".svg", ".png", ".csv", ".md"}
        if not required <= extensions:
            missing_supplements.append({"supplement": number, "missing_extensions": sorted(required - extensions)})
    checks.append({"check": "supplementary_figure_package", "status": "pass" if not missing_supplements else "fail", "detail": missing_supplements})
    args.report.parent.mkdir(parents=True, exist_ok=True)
    report = {"overall_status": "pass" if all(row["status"] == "pass" for row in checks) else "fail", "checks": checks}
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if report["overall_status"] != "pass":
        raise SystemExit("Realstudy v2 integrity validation failed; inspect the JSON report.")


main()
