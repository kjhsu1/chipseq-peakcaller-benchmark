"""Build the normalized transactional Realstudy v2 SQLite database."""

"""Imports"""

import argparse
import csv
import hashlib
import json
import sqlite3
from pathlib import Path

import pandas as pd


"""Constants"""

SCHEMA_VERSION = "1.0.0"
ANALYSIS_TABLES = [
    "parent_libraries", "software_versions", "sampling_blocks", "control_subsamples",
    "artifacts", "validation_events", "peaks", "reference_regions", "peak_overlaps",
    "reference_region_recovery", "run_metrics", "seed_pair_metrics", "stratified_metrics",
    "enough_control_decisions",
]

SCHEMA_SQL = """
PRAGMA foreign_keys = ON;
CREATE TABLE schema_metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE publications (
  citation_id TEXT PRIMARY KEY, title TEXT NOT NULL, authors TEXT NOT NULL, year INTEGER NOT NULL,
  journal TEXT NOT NULL, doi TEXT, pmid TEXT, pmcid TEXT, url TEXT NOT NULL,
  role_in_paper TEXT NOT NULL, classification TEXT NOT NULL
);
CREATE TABLE publication_claims (
  claim_id TEXT PRIMARY KEY, citation_id TEXT NOT NULL REFERENCES publications(citation_id),
  claim_text TEXT NOT NULL, paper_section TEXT NOT NULL
);
CREATE TABLE studies (
  study_id TEXT PRIMARY KEY, organism TEXT NOT NULL, condition TEXT NOT NULL,
  target_or_mark TEXT NOT NULL, assembly TEXT NOT NULL, signal_class TEXT NOT NULL
);
CREATE TABLE experiments (
  experiment_id TEXT PRIMARY KEY, study_id TEXT NOT NULL REFERENCES studies(study_id),
  repository TEXT NOT NULL, role TEXT NOT NULL CHECK(role IN ('treatment','control')),
  matched_experiment_id TEXT REFERENCES experiments(experiment_id) DEFERRABLE INITIALLY DEFERRED,
  antibody_or_tag TEXT, layout TEXT NOT NULL, read_length INTEGER NOT NULL,
  platform TEXT NOT NULL, inclusion_status TEXT NOT NULL
);
CREATE TABLE samples (
  sample_id TEXT PRIMARY KEY, study_id TEXT NOT NULL REFERENCES studies(study_id),
  experiment_id TEXT NOT NULL REFERENCES experiments(experiment_id), role TEXT NOT NULL,
  replicate TEXT NOT NULL
);
CREATE TABLE input_files (
  file_id TEXT PRIMARY KEY, study_id TEXT NOT NULL REFERENCES studies(study_id),
  experiment_id TEXT NOT NULL REFERENCES experiments(experiment_id),
  sample_id TEXT NOT NULL REFERENCES samples(sample_id), role TEXT NOT NULL,
  fastq_url TEXT NOT NULL, raw_reads INTEGER NOT NULL CHECK(raw_reads >= 0),
  raw_bases INTEGER CHECK(raw_bases >= 0), compressed_bytes INTEGER NOT NULL CHECK(compressed_bytes >= 0),
  md5 TEXT NOT NULL CHECK(length(md5) = 32), expected_local_path TEXT NOT NULL,
  staged_checksum TEXT NOT NULL CHECK(length(staged_checksum) IN (32,64)),
  staged_status TEXT NOT NULL CHECK(staged_status = 'verified')
);
CREATE TABLE experiment_publications (
  experiment_id TEXT NOT NULL REFERENCES experiments(experiment_id),
  citation_id TEXT NOT NULL REFERENCES publications(citation_id),
  PRIMARY KEY(experiment_id, citation_id)
);
CREATE TABLE parent_libraries (
  parent_library_id TEXT PRIMARY KEY, study_id TEXT NOT NULL REFERENCES studies(study_id),
  role TEXT NOT NULL CHECK(role IN ('treatment','control')), bam_path TEXT NOT NULL,
  source_checksum TEXT NOT NULL, eligible_fragments INTEGER NOT NULL CHECK(eligible_fragments >= 0),
  aligned_reads INTEGER NOT NULL CHECK(aligned_reads >= 0), duplicate_marked_fragments INTEGER NOT NULL CHECK(duplicate_marked_fragments >= 0),
  status TEXT NOT NULL
);
CREATE TABLE software_versions (software_id TEXT PRIMARY KEY, software_name TEXT NOT NULL, version TEXT NOT NULL, command TEXT NOT NULL);
CREATE TABLE sampling_blocks (
  matched_block_id TEXT PRIMARY KEY, study_id TEXT NOT NULL REFERENCES studies(study_id),
  seed INTEGER,
  reference_run_id TEXT NOT NULL REFERENCES runs(run_id) DEFERRABLE INITIALLY DEFERRED,
  algorithm_version TEXT NOT NULL
);
CREATE TABLE control_subsamples (
  control_subsample_id TEXT PRIMARY KEY,
  control_parent_id TEXT NOT NULL REFERENCES parent_libraries(parent_library_id) DEFERRABLE INITIALLY DEFERRED,
  study_id TEXT NOT NULL REFERENCES studies(study_id), seed INTEGER NOT NULL,
  control_coverage_x REAL NOT NULL CHECK(control_coverage_x > 0),
  requested_fragments INTEGER NOT NULL CHECK(requested_fragments >= 0),
  realized_fragments INTEGER CHECK(realized_fragments >= 0), realized_reads INTEGER CHECK(realized_reads >= 0),
  source_checksum TEXT NOT NULL, output_checksum TEXT, rank_cutoff TEXT,
  UNIQUE(study_id, seed, control_coverage_x)
);
CREATE TABLE runs (
  run_id TEXT PRIMARY KEY,
  matched_block_id TEXT NOT NULL REFERENCES sampling_blocks(matched_block_id) DEFERRABLE INITIALLY DEFERRED,
  study_id TEXT NOT NULL REFERENCES studies(study_id),
  run_type TEXT NOT NULL CHECK(run_type IN ('full_control_anchor','control_subsample')),
  reference_run_id TEXT NOT NULL REFERENCES runs(run_id) DEFERRABLE INITIALLY DEFERRED,
  seed INTEGER, control_coverage_x TEXT NOT NULL,
  requested_control_fragments INTEGER CHECK(requested_control_fragments >= 0),
  treatment_parent_id TEXT NOT NULL REFERENCES parent_libraries(parent_library_id) DEFERRABLE INITIALLY DEFERRED,
  control_parent_id TEXT NOT NULL REFERENCES parent_libraries(parent_library_id) DEFERRABLE INITIALLY DEFERRED,
  control_subsample_id TEXT REFERENCES control_subsamples(control_subsample_id) DEFERRABLE INITIALLY DEFERRED,
  signal_class TEXT NOT NULL CHECK(signal_class IN ('narrow','broad')),
  macs2_mode TEXT NOT NULL CHECK(macs2_mode IN ('narrow','broad')), status TEXT NOT NULL,
  CHECK((run_type='full_control_anchor' AND seed IS NULL AND control_subsample_id IS NULL) OR
        (run_type='control_subsample' AND seed IS NOT NULL AND control_subsample_id IS NOT NULL)),
  UNIQUE(study_id, seed, control_coverage_x)
);
CREATE TABLE artifacts (
  artifact_id TEXT PRIMARY KEY, run_id TEXT REFERENCES runs(run_id), artifact_type TEXT NOT NULL,
  path TEXT NOT NULL, checksum TEXT, bytes INTEGER CHECK(bytes >= 0), status TEXT NOT NULL
);
CREATE TABLE validation_events (
  validation_id TEXT PRIMARY KEY, run_id TEXT REFERENCES runs(run_id), check_name TEXT NOT NULL,
  status TEXT NOT NULL CHECK(status IN ('pass','fail','not_run')), detail TEXT NOT NULL
);
CREATE TABLE peaks (
  peak_id TEXT PRIMARY KEY, run_id TEXT NOT NULL REFERENCES runs(run_id), chromosome TEXT NOT NULL,
  start INTEGER NOT NULL CHECK(start >= 0), end INTEGER NOT NULL CHECK(end > start),
  score REAL, summit INTEGER, width INTEGER NOT NULL CHECK(width > 0)
);
CREATE TABLE reference_regions (
  reference_region_id TEXT PRIMARY KEY, study_id TEXT NOT NULL REFERENCES studies(study_id),
  reference_run_id TEXT NOT NULL REFERENCES runs(run_id), chromosome TEXT NOT NULL,
  start INTEGER NOT NULL CHECK(start >= 0), end INTEGER NOT NULL CHECK(end > start),
  ontology_class TEXT, score REAL, summit INTEGER
);
CREATE TABLE peak_overlaps (
  overlap_id TEXT PRIMARY KEY, run_id TEXT NOT NULL REFERENCES runs(run_id),
  peak_id TEXT NOT NULL REFERENCES peaks(peak_id), reference_region_id TEXT NOT NULL REFERENCES reference_regions(reference_region_id),
  overlap_bp INTEGER NOT NULL CHECK(overlap_bp > 0), query_fraction REAL NOT NULL,
  reference_fraction REAL NOT NULL, jaccard REAL NOT NULL, threshold REAL NOT NULL
);
CREATE TABLE reference_region_recovery (
  run_id TEXT NOT NULL REFERENCES runs(run_id), reference_region_id TEXT NOT NULL REFERENCES reference_regions(reference_region_id),
  threshold REAL NOT NULL, recovered INTEGER NOT NULL CHECK(recovered IN (0,1)),
  PRIMARY KEY(run_id, reference_region_id, threshold)
);
CREATE TABLE run_metrics (
  run_id TEXT NOT NULL REFERENCES runs(run_id), metric TEXT NOT NULL, threshold REAL NOT NULL,
  value REAL, status TEXT NOT NULL, PRIMARY KEY(run_id, metric, threshold)
);
CREATE TABLE seed_pair_metrics (
  study_id TEXT NOT NULL REFERENCES studies(study_id), control_coverage_x REAL NOT NULL,
  seed_a INTEGER NOT NULL, seed_b INTEGER NOT NULL, metric TEXT NOT NULL, value REAL,
  PRIMARY KEY(study_id, control_coverage_x, seed_a, seed_b, metric)
);
CREATE TABLE stratified_metrics (
  run_id TEXT NOT NULL REFERENCES runs(run_id), stratum_type TEXT NOT NULL, stratum_value TEXT NOT NULL,
  metric TEXT NOT NULL, threshold REAL NOT NULL, value REAL,
  PRIMARY KEY(run_id, stratum_type, stratum_value, metric, threshold)
);
CREATE TABLE enough_control_decisions (
  study_id TEXT PRIMARY KEY REFERENCES studies(study_id), outcome TEXT NOT NULL CHECK(outcome IN
  ('sufficient_depth_identified','no_tested_depth_sufficient','agreement_without_plateau','parent_library_insufficient')),
  selected_coverage_x REAL, passing_seeds INTEGER, rationale TEXT NOT NULL
);
CREATE INDEX publications_by_year ON publications(year);
CREATE INDEX input_files_by_study ON input_files(study_id, role);
CREATE INDEX runs_by_design ON runs(study_id, control_coverage_x, seed);
CREATE INDEX peaks_by_interval ON peaks(run_id, chromosome, start, end);
CREATE INDEX reference_by_interval ON reference_regions(study_id, chromosome, start, end);
CREATE INDEX metrics_by_name ON run_metrics(metric, threshold);
CREATE INDEX overlap_by_run ON peak_overlaps(run_id, threshold);
"""


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse database-builder arguments."""
    parser = argparse.ArgumentParser(description="Build the Realstudy v2 SQLite database.")
    parser.add_argument("--publications", type=Path, required=True)
    parser.add_argument("--experiments", type=Path, required=True)
    parser.add_argument("--files", type=Path, required=True)
    parser.add_argument("--run-table", type=Path, required=True)
    parser.add_argument("--table", action="append", default=[], help="table_name=csv_path")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--export-dir", type=Path, required=True)
    parser.add_argument("--export-manifest", type=Path, required=True)
    parser.add_argument("--mark-runs-complete", action="store_true")
    parser.add_argument("--artifact-root", type=Path)
    return parser.parse_args()


def csv_rows(path: Path) -> list[dict]:
    """Read CSV rows with blanks preserved as empty strings."""
    return pd.read_csv(path, dtype=str).fillna("").to_dict(orient="records")


def insert_dicts(connection: sqlite3.Connection, table: str, rows: list[dict]) -> None:
    """Insert dictionaries using the intersection with declared table columns."""
    if not rows:
        return
    allowed = {row[1] for row in connection.execute(f"PRAGMA table_info({table})")}
    columns = [column for column in rows[0] if column in allowed]
    if not columns:
        raise ValueError(f"No CSV columns map to database table {table}.")
    sql = f"INSERT INTO {table} ({','.join(columns)}) VALUES ({','.join('?' for _ in columns)})"
    values = []
    for row in rows:
        values.append(tuple(None if row[column] == "" else row[column] for column in columns))
    connection.executemany(sql, values)


def load_registries(connection: sqlite3.Connection, args: argparse.Namespace) -> None:
    """Normalize the versioned publication, experiment, file, and run registries."""
    publications = csv_rows(args.publications)
    insert_dicts(connection, "publications", publications)
    claim_rows = [
        {
            "claim_id": f"{row['citation_id']}__registered_claim",
            "citation_id": row["citation_id"],
            "claim_text": row["supported_claim_or_method"],
            "paper_section": row["role_in_paper"],
        }
        for row in publications
    ]
    insert_dicts(connection, "publication_claims", claim_rows)
    experiments = csv_rows(args.experiments)
    study_rows = {}
    for row in experiments:
        if row["study_id"] in study_rows and row["role"] != "treatment":
            continue
        signal_class = "broad" if "H3K9" in row["target_or_mark"] else "narrow"
        study_rows[row["study_id"]] = {
            "study_id": row["study_id"], "organism": row["organism"],
            "condition": row["condition"], "target_or_mark": row["target_or_mark"],
            "assembly": "ce11", "signal_class": signal_class,
        }
    insert_dicts(connection, "studies", list(study_rows.values()))
    insert_dicts(connection, "experiments", experiments)
    insert_dicts(
        connection,
        "experiment_publications",
        [{"experiment_id": row["experiment_id"], "citation_id": row["publication_id"]} for row in experiments],
    )
    files = csv_rows(args.files)
    samples = {
        row["sample_id"]: {
            "sample_id": row["sample_id"], "study_id": row["study_id"],
            "experiment_id": row["experiment_id"], "role": row["role"], "replicate": row["replicate"],
        }
        for row in files
    }
    insert_dicts(connection, "samples", list(samples.values()))
    insert_dicts(connection, "input_files", files)
    runs = csv_rows(args.run_table)
    for row in runs:
        row["status"] = "complete" if args.mark_runs_complete else (row.get("status") or "planned")
    insert_dicts(connection, "runs", runs)


def export_database(connection: sqlite3.Connection, export_dir: Path, manifest_path: Path) -> None:
    """Export every table deterministically with row counts and SHA-256 checksums."""
    export_dir.mkdir(parents=True, exist_ok=True)
    tables = [row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")]
    manifest_rows = []
    for table in tables:
        columns = [row[1] for row in connection.execute(f"PRAGMA table_info({table})")]
        records = connection.execute(f"SELECT * FROM {table} ORDER BY rowid").fetchall()
        path = export_dir / f"{table}.csv"
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(columns)
            writer.writerows(records)
        manifest_rows.append(
            {"table_name": table, "row_count": len(records), "csv_path": str(path),
             "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        )
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(manifest_rows).to_csv(manifest_path, index=False)


def load_rendered_artifacts(connection: sqlite3.Connection, root: Path) -> None:
    """Register rendered figure/source/caption artifacts after figure completion."""
    rows = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.name.startswith("."):
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        relative = path.relative_to(root)
        rows.append(
            {"artifact_id": f"figure__{hashlib.sha256(str(relative).encode()).hexdigest()[:20]}",
             "run_id": "", "artifact_type": f"paper_{path.suffix.lstrip('.')}",
             "path": str(path), "checksum": digest, "bytes": path.stat().st_size, "status": "complete"}
        )
    insert_dicts(connection, "artifacts", rows)


def main() -> None:
    """Build into a temporary database and atomically publish it after validation."""
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    if temporary.exists():
        temporary.unlink()
    connection = sqlite3.connect(temporary)
    try:
        connection.executescript(SCHEMA_SQL)
        connection.execute("BEGIN")
        connection.executemany(
            "INSERT INTO schema_metadata VALUES (?, ?)",
            [("schema_version", SCHEMA_VERSION), ("foreign_keys", "enabled")],
        )
        load_registries(connection, args)
        for item in args.table:
            table, separator, path = item.partition("=")
            if not separator or table not in ANALYSIS_TABLES:
                raise ValueError(f"Unsupported --table mapping: {item}")
            insert_dicts(connection, table, csv_rows(Path(path)))
        if args.artifact_root:
            load_rendered_artifacts(connection, args.artifact_root)
        violations = connection.execute("PRAGMA foreign_key_check").fetchall()
        if violations:
            raise ValueError(f"Foreign-key violations: {json.dumps(violations)}")
        connection.commit()
        export_database(connection, args.export_dir, args.export_manifest)
    except Exception:
        connection.rollback()
        connection.close()
        if temporary.exists():
            temporary.unlink()
        raise
    connection.close()
    temporary.replace(args.output)


main()
