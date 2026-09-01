"""Create exact deterministic nested empirical-control BAM subsets."""

"""Imports"""

import argparse
import csv
import hashlib
import sqlite3
from pathlib import Path

import pysam

from realstudy_sampling_lib import eligible_alignment, fragment_identifier, fragment_rank_digest


"""Constants"""

ALGORITHM_VERSION = "sha256-rg-qname-v1"


"""Functions"""


def parse_key_value(text: str, value_type) -> tuple[str, object]:
    """Parse one label=value command-line item."""
    label, separator, value = text.partition("=")
    if not separator or not label or not value:
        raise argparse.ArgumentTypeError(f"Expected label=value, observed {text!r}")
    return label, value_type(value)


def parse_args() -> argparse.Namespace:
    """Parse nested-sampling arguments."""
    parser = argparse.ArgumentParser(description="Build exact nested control BAMs.")
    parser.add_argument("--parent-bam", type=Path, required=True)
    parser.add_argument("--source-checksum")
    parser.add_argument("--source-checksum-file", type=Path)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--study-id", required=True)
    parser.add_argument("--control-parent-id", required=True)
    parser.add_argument("--genome-size-bp", type=int, required=True)
    parser.add_argument("--normalization-fragment-length-bp", type=int, required=True)
    parser.add_argument("--minimum-eligible-fragments", type=int, required=True)
    parser.add_argument("--target", action="append", required=True, help="coverage_label=fragment_count")
    parser.add_argument("--output", action="append", required=True, help="coverage_label=output.bam")
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--algorithm-version", default=ALGORITHM_VERSION)
    return parser.parse_args()


def eligible_record(record: pysam.AlignedSegment) -> bool:
    """Return whether an alignment is eligible; duplicate-marked reads remain eligible."""
    return eligible_alignment(record)


def fragment_id(record: pysam.AlignedSegment) -> str:
    """Return stable read-group plus query-name fragment identity."""
    return fragment_identifier(record)


def rank_digest(version: str, checksum: str, seed: int, identifier: str) -> str:
    """Return the deterministic SHA-256 fragment rank."""
    return fragment_rank_digest(version, checksum, seed, identifier)


def file_checksum(path: Path, algorithm: str = "sha256") -> str:
    """Calculate a streaming file checksum."""
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def initialize_ledger(path: Path) -> sqlite3.Connection:
    """Create a fresh disk-backed fragment ranking ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA synchronous=NORMAL")
    connection.execute(
        "CREATE TABLE fragments (fragment_id TEXT PRIMARY KEY, rank_hex TEXT NOT NULL, duplicate_marked INTEGER NOT NULL)"
    )
    return connection


def populate_ledger(connection: sqlite3.Connection, args: argparse.Namespace) -> tuple[int, int]:
    """Scan the parent once to register eligible fragments and duplicate burden."""
    insert = (
        "INSERT INTO fragments(fragment_id, rank_hex, duplicate_marked) VALUES (?, ?, ?) "
        "ON CONFLICT(fragment_id) DO UPDATE SET duplicate_marked=MAX(duplicate_marked, excluded.duplicate_marked)"
    )
    batch = []
    with pysam.AlignmentFile(args.parent_bam, "rb") as source:
        for record in source.fetch(until_eof=True):
            if not eligible_record(record):
                continue
            identifier = fragment_id(record)
            batch.append(
                (identifier, rank_digest(args.algorithm_version, args.source_checksum, args.seed, identifier), int(record.is_duplicate))
            )
            if len(batch) >= 50000:
                connection.executemany(insert, batch)
                connection.commit()
                batch = []
    if batch:
        connection.executemany(insert, batch)
    connection.execute("CREATE INDEX fragments_by_rank ON fragments(rank_hex, fragment_id)")
    connection.commit()
    eligible, duplicates = connection.execute(
        "SELECT COUNT(*), COALESCE(SUM(duplicate_marked), 0) FROM fragments"
    ).fetchone()
    return int(eligible), int(duplicates)


def mark_selected(connection: sqlite3.Connection, targets: list[tuple[str, int]]) -> dict[str, str]:
    """Assign each selected fragment to its earliest nested depth."""
    connection.execute(
        "CREATE TABLE selected (fragment_id TEXT PRIMARY KEY, first_depth INTEGER NOT NULL, rank_hex TEXT NOT NULL)"
    )
    maximum = max(count for _, count in targets)
    rows = connection.execute(
        "SELECT fragment_id, rank_hex FROM fragments ORDER BY rank_hex, fragment_id LIMIT ?", (maximum,)
    )
    cutoffs = {}
    batch = []
    target_index = 0
    for position, (identifier, rank_hex) in enumerate(rows, start=1):
        while target_index < len(targets) - 1 and position > targets[target_index][1]:
            target_index += 1
        batch.append((identifier, target_index, rank_hex))
        for label, count in targets:
            if position == count:
                cutoffs[label] = rank_hex
        if len(batch) >= 50000:
            connection.executemany("INSERT INTO selected VALUES (?, ?, ?)", batch)
            connection.commit()
            batch = []
    if batch:
        connection.executemany("INSERT INTO selected VALUES (?, ?, ?)", batch)
    connection.commit()
    if len(cutoffs) != len(targets):
        raise ValueError("Could not derive all requested fragment-rank cutoffs.")
    return cutoffs


def write_nested_bams(
    connection: sqlite3.Connection,
    parent_bam: Path,
    targets: list[tuple[str, int]],
    outputs: dict[str, Path],
) -> dict[str, int]:
    """Write all nested BAMs in one second pass while preserving coordinate order."""
    handles = {}
    read_counts = {label: 0 for label, _ in targets}
    lookup = connection.cursor()
    with pysam.AlignmentFile(parent_bam, "rb") as source:
        for label, _ in targets:
            outputs[label].parent.mkdir(parents=True, exist_ok=True)
            handles[label] = pysam.AlignmentFile(outputs[label], "wb", template=source)
        try:
            for record in source.fetch(until_eof=True):
                if not eligible_record(record):
                    continue
                selected = lookup.execute(
                    "SELECT first_depth FROM selected WHERE fragment_id = ?", (fragment_id(record),)
                ).fetchone()
                if selected is None:
                    continue
                for index in range(int(selected[0]), len(targets)):
                    label = targets[index][0]
                    handles[label].write(record)
                    read_counts[label] += 1
        finally:
            for handle in handles.values():
                handle.close()
    for label, _ in targets:
        pysam.index(str(outputs[label]))
    return read_counts


def verify_selected_counts(connection: sqlite3.Connection, targets: list[tuple[str, int]]) -> None:
    """Require exact cumulative fragment counts at every nested depth."""
    previous = 0
    for index, (label, expected) in enumerate(targets):
        observed = connection.execute(
            "SELECT COUNT(*) FROM selected WHERE first_depth <= ?", (index,)
        ).fetchone()[0]
        if int(observed) != int(expected) or int(observed) <= previous:
            raise ValueError(
                f"Nested selection verification failed at {label}: expected {expected}, observed {observed}."
            )
        previous = int(observed)


def main() -> None:
    """Build exact nested subsets and their provenance manifest."""
    args = parse_args()
    if bool(args.source_checksum) == bool(args.source_checksum_file):
        raise ValueError("Provide exactly one of --source-checksum or --source-checksum-file.")
    if args.source_checksum_file:
        args.source_checksum = args.source_checksum_file.read_text().strip().split()[0]
    targets = sorted(
        [parse_key_value(item, int) for item in args.target], key=lambda item: int(item[1])
    )
    outputs = dict(parse_key_value(item, Path) for item in args.output)
    if set(outputs) != {label for label, _ in targets}:
        raise ValueError("Every target label must have exactly one output BAM.")
    if file_checksum(args.parent_bam) != args.source_checksum.lower():
        raise ValueError("Parent BAM SHA-256 does not match --source-checksum.")
    connection = initialize_ledger(args.ledger)
    try:
        eligible, duplicate_marked = populate_ledger(connection, args)
        if eligible < args.minimum_eligible_fragments or eligible < targets[-1][1]:
            raise ValueError(
                f"Parent library has {eligible} eligible fragments; locked design requires "
                f"at least {max(args.minimum_eligible_fragments, targets[-1][1])}."
            )
        cutoffs = mark_selected(connection, targets)
        verify_selected_counts(connection, targets)
        read_counts = write_nested_bams(connection, args.parent_bam, targets, outputs)
    finally:
        connection.close()
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    with args.manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "control_subsample_id", "control_parent_id", "study_id", "seed", "control_coverage_x", "coverage_label", "requested_fragments", "realized_fragments",
                "realized_reads", "realized_coverage_x", "eligible_parent_fragments",
                "duplicate_marked_parent_fragments", "source_sha256", "output_bam",
                "output_sha256", "rank_cutoff", "algorithm_version", "nested_verified",
            ],
        )
        writer.writeheader()
        for label, count in targets:
            writer.writerow(
                {
                    "control_subsample_id": f"{args.study_id}__control_{label}__seed_{args.seed}",
                    "control_parent_id": args.control_parent_id,
                    "study_id": args.study_id,
                    "seed": args.seed,
                    "control_coverage_x": label.removesuffix("x").replace("p", "."),
                    "coverage_label": label,
                    "requested_fragments": count,
                    "realized_fragments": count,
                    "realized_reads": read_counts[label],
                    "realized_coverage_x": (count * args.normalization_fragment_length_bp) / args.genome_size_bp,
                    "eligible_parent_fragments": eligible,
                    "duplicate_marked_parent_fragments": duplicate_marked,
                    "source_sha256": args.source_checksum.lower(),
                    "output_bam": outputs[label],
                    "output_sha256": file_checksum(outputs[label]),
                    "rank_cutoff": cutoffs[label],
                    "algorithm_version": args.algorithm_version,
                    "nested_verified": True,
                }
            )


main()
