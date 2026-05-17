"""Download one manifest-specified real-study file."""

"""Imports"""

import argparse
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd
import requests


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Download one file from the real-study manifest.")
    parser.add_argument("--manifest-csv", type=Path, default=Path("metadata/data_manifest.csv"))
    parser.add_argument("--study-id", required=True)
    parser.add_argument("--role", required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("data/raw"))
    parser.add_argument("--output-path", type=Path, default=None)
    parser.add_argument("--chunk-size", type=int, default=1024 * 1024)
    parser.add_argument("--connect-timeout", type=int, default=60)
    parser.add_argument("--read-timeout", type=int, default=300)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace any existing local file instead of resuming/skipping it.",
    )
    return parser.parse_args()


def filename_from_url(url: str) -> str:
    """Return a filename from the remote URL."""
    parsed = urlparse(url)
    name = Path(parsed.path).name
    return name or "download.bin"


def choose_transfer_mode(
    existing_size: int, remote_size: int | None, accept_ranges: str, overwrite: bool
) -> str:
    """Choose whether to skip, restart, or resume a download."""
    if overwrite:
        return "restart"
    if existing_size <= 0:
        return "start"
    if remote_size is not None and existing_size == remote_size:
        return "skip"
    if remote_size is not None and existing_size > remote_size:
        return "restart"
    if "bytes" in accept_ranges.lower():
        return "resume"
    return "restart"


def remote_metadata(session: requests.Session, remote_url: str, timeout: tuple[int, int]) -> tuple[int | None, str]:
    """Return remote size and range support when HEAD is accepted."""
    try:
        response = session.head(remote_url, allow_redirects=True, timeout=timeout)
        response.raise_for_status()
    except requests.RequestException:
        return None, ""
    size_text = response.headers.get("Content-Length", "").strip()
    remote_size = int(size_text) if size_text.isdigit() else None
    accept_ranges = response.headers.get("Accept-Ranges", "")
    return remote_size, accept_ranges


def main() -> None:
    """Download the requested manifest row into the raw-data tree."""
    args = parse_args()
    manifest = pd.read_csv(args.manifest_csv).fillna("")
    row = manifest[
        (manifest["study_id"] == args.study_id)
        & (manifest["role"] == args.role)
    ]
    if row.empty:
        raise SystemExit(f"Could not find manifest row for {args.study_id} / {args.role}")
    record = row.iloc[0]
    remote_url = str(record["remote_url"]).strip()
    if not remote_url:
        raise SystemExit("Manifest row has no remote_url")
    if args.output_path is not None:
        out_path = args.output_path
        out_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        study_dir = args.output_dir / args.study_id
        study_dir.mkdir(parents=True, exist_ok=True)
        out_path = study_dir / filename_from_url(remote_url)
    session = requests.Session()
    timeout = (args.connect_timeout, args.read_timeout)
    remote_size, accept_ranges = remote_metadata(session, remote_url, timeout)
    existing_size = out_path.stat().st_size if out_path.exists() else 0
    mode = choose_transfer_mode(existing_size, remote_size, accept_ranges, args.overwrite)
    if mode == "skip":
        print(str(out_path))
        return
    headers = {}
    file_mode = "wb"
    if mode == "resume":
        headers["Range"] = f"bytes={existing_size}-"
        file_mode = "ab"
    with session.get(remote_url, stream=True, timeout=timeout, headers=headers) as response:
        response.raise_for_status()
        with out_path.open(file_mode) as handle:
            for chunk in response.iter_content(chunk_size=args.chunk_size):
                if chunk:
                    handle.write(chunk)
    if remote_size is not None and out_path.stat().st_size != remote_size:
        raise SystemExit(
            f"Downloaded size mismatch for {out_path.name}: "
            f"expected {remote_size}, observed {out_path.stat().st_size}"
        )
    print(str(out_path))


if __name__ == "__main__":
    main()
