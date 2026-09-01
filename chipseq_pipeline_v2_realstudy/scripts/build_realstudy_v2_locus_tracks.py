"""Build representative-locus signal tables from actual Realstudy v2 BAMs."""

"""Imports"""

import argparse
from pathlib import Path

import pandas as pd
import pysam


"""Functions"""


def parse_labeled_path(text: str) -> tuple[str, Path]:
    """Parse label=path arguments."""
    label, separator, path = text.partition("=")
    if not separator:
        raise argparse.ArgumentTypeError("Expected label=path")
    return label, Path(path)


def parse_args() -> argparse.Namespace:
    """Parse representative-locus arguments."""
    parser = argparse.ArgumentParser(description="Build representative locus signal source data.")
    parser.add_argument("--study-id", required=True)
    parser.add_argument("--reference-regions", type=Path, required=True)
    parser.add_argument("--track", action="append", required=True, help="track_label=bam_path")
    parser.add_argument("--peak-track", action="append", default=[], help="track_label=peak_bed_path")
    parser.add_argument("--flank-bp", type=int, default=2500)
    parser.add_argument("--bin-bp", type=int, default=50)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Choose the strongest anchor region and emit binned read-start coverage for each track."""
    args = parse_args()
    regions = pd.read_csv(args.reference_regions)
    regions = regions[regions["study_id"] == args.study_id].copy()
    if regions.empty:
        raise ValueError(f"No anchor regions available for {args.study_id}.")
    regions["score"] = pd.to_numeric(regions["score"], errors="coerce").fillna(0)
    anchor = regions.sort_values(["score", "end"], ascending=[False, True]).iloc[0]
    start = max(0, int(anchor["start"]) - args.flank_bp)
    end = int(anchor["end"]) + args.flank_bp
    tracks = [parse_labeled_path(item) for item in args.track]
    peak_tracks = [parse_labeled_path(item) for item in args.peak_track]
    rows = []
    for order, (label, path) in enumerate(tracks):
        with pysam.AlignmentFile(path, "rb") as bam:
            mapped_reads = max(1, bam.mapped)
            for bin_start in range(start, end, args.bin_bp):
                bin_end = min(end, bin_start + args.bin_bp)
                count = bam.count(str(anchor["chromosome"]), bin_start, bin_end, read_callback="all")
                rows.append(
                    {"study_id": args.study_id, "locus_id": f"{args.study_id}:{anchor['chromosome']}:{start}-{end}",
                     "track_label": label, "track_order": order, "position": (bin_start + bin_end) // 2,
                     "value": (count * 1_000_000) / (mapped_reads * max(1, bin_end - bin_start))}
                )
    for offset, (label, path) in enumerate(peak_tracks, start=len(tracks)):
        peak_frame = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2]) if path.stat().st_size else pd.DataFrame(columns=[0, 1, 2])
        peak_frame = peak_frame[peak_frame[0].astype(str) == str(anchor["chromosome"])]
        for bin_start in range(start, end, args.bin_bp):
            bin_end = min(end, bin_start + args.bin_bp)
            called = int(((peak_frame[1] < bin_end) & (peak_frame[2] > bin_start)).any()) if not peak_frame.empty else 0
            rows.append(
                {"study_id": args.study_id, "locus_id": f"{args.study_id}:{anchor['chromosome']}:{start}-{end}",
                 "track_label": label, "track_order": offset, "position": (bin_start + bin_end) // 2,
                 "value": called}
            )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False)


main()
