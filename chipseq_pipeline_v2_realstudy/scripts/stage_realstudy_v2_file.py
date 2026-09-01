"""Stage or verify one registry-defined Realstudy v2 FASTQ."""

"""Imports"""

import argparse
import hashlib
from pathlib import Path

import pandas as pd
import requests


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse file-staging arguments."""
    parser = argparse.ArgumentParser(description="Stage one Realstudy v2 FASTQ with MD5 verification.")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--file-id", required=True)
    parser.add_argument("--marker", type=Path, required=True)
    parser.add_argument("--require-local", action="store_true")
    parser.add_argument("--chunk-size", type=int, default=1024 * 1024)
    return parser.parse_args()


def md5(path: Path) -> str:
    """Return a streaming MD5 for repository checksum verification."""
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    """Require a local FASTQ or download it explicitly, then verify its exact checksum."""
    args = parse_args()
    manifest = pd.read_csv(args.manifest, dtype=str).fillna("")
    selected = manifest[manifest["file_id"] == args.file_id]
    if len(selected) != 1:
        raise ValueError(f"Expected one registry row for {args.file_id}; observed {len(selected)}.")
    row = selected.iloc[0]
    output = Path(row["expected_local_path"])
    if not output.exists():
        if args.require_local:
            raise SystemExit(f"Missing required pre-staged FASTQ: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        temporary = output.with_suffix(output.suffix + ".partial")
        with requests.get(row["fastq_url"], stream=True, timeout=(60, 300)) as response:
            response.raise_for_status()
            with temporary.open("wb") as handle:
                for chunk in response.iter_content(args.chunk_size):
                    if chunk:
                        handle.write(chunk)
        temporary.replace(output)
    observed = md5(output)
    if observed != row["md5"]:
        raise ValueError(f"MD5 mismatch for {args.file_id}: expected {row['md5']}, observed {observed}")
    args.marker.parent.mkdir(parents=True, exist_ok=True)
    args.marker.write_text(f"file_id\tpath\tmd5\n{args.file_id}\t{output}\t{observed}\n")


main()
