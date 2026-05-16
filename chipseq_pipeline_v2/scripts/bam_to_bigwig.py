#!/usr/bin/env python3
"""Convert a coordinate-sorted BAM file to BigWig using deepTools bamCoverage."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a BigWig from a sorted BAM + BAI using deepTools bamCoverage.",
        epilog=(
            "Example:\n"
            "  python scripts/bam_to_bigwig.py "
            "--bam aligned.sorted.bam "
            "--bai aligned.sorted.bam.bai "
            "--out aligned.sorted.bw"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--bam", required=True, help="Path to coordinate-sorted BAM file.")
    parser.add_argument("--bai", required=True, help="Path to BAM index (.bai).")
    parser.add_argument(
        "--out",
        required=True,
        help="Output BigWig path (.bw or .bigwig).",
    )
    parser.add_argument(
        "--bin-size",
        type=int,
        default=10,
        help="Bin size for bamCoverage (default: 10).",
    )
    parser.add_argument(
        "--normalize-using",
        default="RPKM",
        help="Normalization method passed to bamCoverage (default: RPKM).",
    )
    parser.add_argument(
        "--effective-genome-size",
        type=int,
        default=None,
        help="Effective genome size; only passed if provided.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads for bamCoverage (default: 1).",
    )
    return parser.parse_args()


def fail(message: str, code: int = 1) -> int:
    print(f"ERROR: {message}", file=sys.stderr)
    return code


def validate_inputs(bam: Path, bai: Path, out_bw: Path, bin_size: int, threads: int) -> int:
    if not bam.exists():
        return fail(f"BAM file not found: {bam}")
    if not bai.exists():
        return fail(f"BAI file not found: {bai}")
    if bam.suffix.lower() != ".bam":
        return fail(f"BAM file must end with .bam: {bam}")

    out_suffix = out_bw.suffix.lower()
    if out_suffix not in {".bw", ".bigwig"}:
        return fail(f"Output file must end with .bw or .bigwig: {out_bw}")
    if bin_size <= 0:
        return fail(f"--bin-size must be > 0 (got {bin_size}).")
    if threads <= 0:
        return fail(f"--threads must be > 0 (got {threads}).")

    return 0


def require_tool(name: str, hint: str) -> int:
    if shutil.which(name) is None:
        return fail(f"Required executable '{name}' not found in PATH. {hint}")
    return 0


def ensure_coordinate_sorted(bam: Path) -> int:
    cmd = ["samtools", "view", "-H", str(bam)]
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        return fail(
            f"Could not read BAM header with samtools (exit {exc.returncode}). "
            f"stderr: {exc.stderr.strip()}"
        )

    so_value = None
    for line in result.stdout.splitlines():
        if not line.startswith("@HD"):
            continue
        for field in line.split("\t"):
            if field.startswith("SO:"):
                so_value = field.split(":", 1)[1]
                break
        break

    if so_value != "coordinate":
        if so_value is None:
            return fail(
                "BAM header @HD does not contain SO:coordinate. "
                "Please sort BAM by coordinate and re-index."
            )
        return fail(
            f"BAM sort order is '{so_value}', not 'coordinate'. "
            "Please sort BAM by coordinate and re-index."
        )

    return 0


def run_bamcoverage(
    bam: Path,
    out_bw: Path,
    bin_size: int,
    normalize_using: str,
    effective_genome_size: int | None,
    threads: int,
) -> int:
    out_bw.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "bamCoverage",
        "-b",
        str(bam),
        "-o",
        str(out_bw),
        "--binSize",
        str(bin_size),
        "--normalizeUsing",
        normalize_using,
        "-p",
        str(threads),
    ]
    if effective_genome_size is not None:
        cmd.extend(["--effectiveGenomeSize", str(effective_genome_size)])

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        return fail(f"bamCoverage failed with exit code {exc.returncode}.")

    if not out_bw.exists() or out_bw.stat().st_size == 0:
        return fail(f"BigWig was not created or is empty: {out_bw}")

    print(f"BigWig created: {out_bw}")
    return 0


def main() -> int:
    args = parse_args()
    bam = Path(args.bam)
    bai = Path(args.bai)
    out_bw = Path(args.out)

    code = validate_inputs(bam, bai, out_bw, args.bin_size, args.threads)
    if code != 0:
        return code

    code = require_tool(
        "bamCoverage",
        "Install with: conda install -c bioconda deeptools",
    )
    if code != 0:
        return code

    code = require_tool(
        "samtools",
        "Install with: conda install -c bioconda samtools",
    )
    if code != 0:
        return code

    code = ensure_coordinate_sorted(bam)
    if code != 0:
        return code

    return run_bamcoverage(
        bam=bam,
        out_bw=out_bw,
        bin_size=args.bin_size,
        normalize_using=args.normalize_using,
        effective_genome_size=args.effective_genome_size,
        threads=args.threads,
    )


if __name__ == "__main__":
    sys.exit(main())
