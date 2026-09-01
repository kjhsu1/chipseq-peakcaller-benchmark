"""Summarize one fixed treatment or control parent BAM."""

"""Imports"""

import argparse
import hashlib
import sqlite3
from pathlib import Path

import pandas as pd
import pysam


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse parent-library QC arguments."""
    parser = argparse.ArgumentParser(description="Summarize an empirical parent BAM.")
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--parent-library-id", required=True)
    parser.add_argument("--study-id", required=True)
    parser.add_argument("--role", choices=["treatment", "control"], required=True)
    parser.add_argument("--minimum-eligible-fragments", type=int, default=0)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--checksum-output", type=Path, required=True)
    return parser.parse_args()


def eligible(record: pysam.AlignedSegment) -> bool:
    """Match the sampler's mapped-primary-QC-pass eligibility rule."""
    return not (record.is_unmapped or record.is_secondary or record.is_supplementary or record.is_qcfail)


def sha256(path: Path) -> str:
    """Return a streaming BAM checksum."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    """Write unique-fragment eligibility, aligned-read, duplicate, and checksum provenance."""
    args = parse_args()
    ledger_path = args.output_csv.with_suffix(".fragment_ledger.sqlite")
    ledger_path.parent.mkdir(parents=True, exist_ok=True)
    if ledger_path.exists():
        ledger_path.unlink()
    connection = sqlite3.connect(ledger_path)
    connection.execute(
        "CREATE TABLE fragments (fragment_id TEXT PRIMARY KEY, duplicate_marked INTEGER NOT NULL)"
    )
    aligned_reads = 0
    with pysam.AlignmentFile(args.bam, "rb") as source:
        batch = []
        for record in source.fetch(until_eof=True):
            if not eligible(record):
                continue
            aligned_reads += 1
            group = record.get_tag("RG") if record.has_tag("RG") else "NO_RG"
            identifier = f"{group}\x1f{record.query_name}"
            batch.append((identifier, int(record.is_duplicate)))
            if len(batch) >= 50000:
                connection.executemany(
                    "INSERT INTO fragments VALUES (?, ?) ON CONFLICT(fragment_id) DO UPDATE SET duplicate_marked=MAX(duplicate_marked, excluded.duplicate_marked)",
                    batch,
                )
                connection.commit()
                batch = []
        if batch:
            connection.executemany(
                "INSERT INTO fragments VALUES (?, ?) ON CONFLICT(fragment_id) DO UPDATE SET duplicate_marked=MAX(duplicate_marked, excluded.duplicate_marked)",
                batch,
            )
    eligible_fragments, duplicate_fragments = connection.execute(
        "SELECT COUNT(*), COALESCE(SUM(duplicate_marked), 0) FROM fragments"
    ).fetchone()
    connection.close()
    ledger_path.unlink()
    status = "eligible" if eligible_fragments >= args.minimum_eligible_fragments else "insufficient"
    digest = sha256(args.bam)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [{
            "parent_library_id": args.parent_library_id, "study_id": args.study_id, "role": args.role,
            "bam_path": str(args.bam), "source_checksum": digest,
            "eligible_fragments": eligible_fragments, "aligned_reads": aligned_reads,
            "duplicate_marked_fragments": duplicate_fragments, "status": status,
        }]
    ).to_csv(args.output_csv, index=False)
    args.checksum_output.parent.mkdir(parents=True, exist_ok=True)
    args.checksum_output.write_text(f"{digest}  {args.bam}\n")
    if status == "insufficient":
        raise ValueError(
            f"{args.study_id} {args.role} parent has {eligible_fragments} eligible fragments; "
            f"requires {args.minimum_eligible_fragments}."
        )


main()
