"""Build simple binned coverage tracks from BAM files."""

"""Imports"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pysam


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Make binned treatment/control tracks from BAM files.")
    parser.add_argument("--treatment-bam", type=Path, required=True)
    parser.add_argument("--control-bam", type=Path, default=None)
    parser.add_argument("--bin-size", type=int, default=25)
    parser.add_argument("--output-prefix", type=Path, required=True)
    return parser.parse_args()


def depth_vector(bam_path: Path, chrom: str, start: int, end: int) -> np.ndarray:
    """Return per-base coverage for a genomic span."""
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        a_arr, c_arr, g_arr, t_arr = bam.count_coverage(chrom, start, end, quality_threshold=0)
    finally:
        bam.close()
    return np.asarray(a_arr) + np.asarray(c_arr) + np.asarray(g_arr) + np.asarray(t_arr)


def main() -> None:
    """Write binned treatment, control, and log2-ratio tables."""
    args = parse_args()
    bam = pysam.AlignmentFile(str(args.treatment_bam), "rb")
    chrom = bam.references[0]
    length = bam.lengths[0]
    bam.close()
    treat = depth_vector(args.treatment_bam, chrom, 0, length)
    control = np.zeros_like(treat)
    if args.control_bam is not None:
        control = depth_vector(args.control_bam, chrom, 0, length)

    starts = np.arange(0, len(treat), args.bin_size)
    rows = []
    for start in starts:
        end = min(start + args.bin_size, len(treat))
        treat_mean = float(np.mean(treat[start:end]))
        control_mean = float(np.mean(control[start:end]))
        rows.append(
            {
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "treatment": treat_mean,
                "control": control_mean,
                "log2_treat_control": float(np.log2((treat_mean + 1.0) / (control_mean + 1.0))),
            }
        )
    df = pd.DataFrame(rows)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_prefix.with_suffix(".csv"), index=False)
    args.output_prefix.with_suffix(".json").write_text(
        json.dumps({"bin_size": args.bin_size, "language": "reference / estimated, not true"}, indent=2),
        encoding="utf-8",
    )


main()
