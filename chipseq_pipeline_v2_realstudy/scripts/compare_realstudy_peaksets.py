"""Compare one sampled-control MACS2 peak set with its full-control anchor."""

"""Imports"""

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse peak-comparison arguments."""
    parser = argparse.ArgumentParser(description="Compute anchor-relative peak stability metrics.")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--study-id", required=True)
    parser.add_argument("--mode", choices=["narrow", "broad"], required=True)
    parser.add_argument("--query-peaks", type=Path, required=True)
    parser.add_argument("--anchor-peaks", type=Path, required=True)
    parser.add_argument("--threshold", action="append", type=float, default=[])
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_peaks(path: Path, id_prefix: str, mode: str) -> pd.DataFrame:
    """Read native narrowPeak/broadPeak or normalized BED, including valid empty files."""
    columns = ["chromosome", "start", "end", "name", "bed_score", "strand", "score", "pvalue", "qvalue", "summit_offset"]
    if path.stat().st_size == 0:
        return pd.DataFrame(columns=["peak_id", "chromosome", "start", "end", "score", "summit", "width"])
    frame = pd.read_csv(path, sep="\t", header=None, comment="#")
    frame.columns = columns[: len(frame.columns)]
    for column in ["start", "end"]:
        frame[column] = pd.to_numeric(frame[column], errors="raise").astype(int)
    if (frame["start"] < 0).any() or (frame["end"] <= frame["start"]).any():
        raise ValueError(f"Illegal BED interval in {path}")
    scores = pd.to_numeric(frame["score"], errors="coerce") if "score" in frame else pd.Series(np.nan, index=frame.index)
    if mode == "narrow" and "summit_offset" in frame:
        offsets = pd.to_numeric(frame["summit_offset"], errors="coerce")
        summits = frame["start"] + offsets
    else:
        summits = pd.Series(np.nan, index=frame.index)
    return pd.DataFrame(
        {
            "peak_id": [f"{id_prefix}__peak_{index + 1}" for index in range(len(frame))],
            "chromosome": frame["chromosome"].astype(str),
            "start": frame["start"], "end": frame["end"], "score": scores,
            "summit": summits, "width": frame["end"] - frame["start"],
        }
    )


def interval_union_bp(frame: pd.DataFrame) -> int:
    """Return the genomic width of a peak-set union."""
    total = 0
    for _, group in frame.sort_values(["chromosome", "start", "end"]).groupby("chromosome"):
        current_start = None
        current_end = None
        for row in group.itertuples():
            if current_end is None or row.start > current_end:
                if current_end is not None:
                    total += current_end - current_start
                current_start, current_end = row.start, row.end
            else:
                current_end = max(current_end, row.end)
        if current_end is not None:
            total += current_end - current_start
    return int(total)


def overlap_rows(query: pd.DataFrame, reference: pd.DataFrame, run_id: str) -> list[dict]:
    """Enumerate interval overlaps with stable identifiers and fractions."""
    rows = []
    for chromosome, query_group in query.groupby("chromosome"):
        reference_group = reference[reference["chromosome"] == chromosome].sort_values("start")
        reference_records = list(reference_group.itertuples())
        left = 0
        for query_row in query_group.sort_values("start").itertuples():
            while left < len(reference_records) and reference_records[left].end <= query_row.start:
                left += 1
            index = left
            while index < len(reference_records) and reference_records[index].start < query_row.end:
                ref_row = reference_records[index]
                overlap = min(query_row.end, ref_row.end) - max(query_row.start, ref_row.start)
                if overlap > 0:
                    union = (query_row.end - query_row.start) + (ref_row.end - ref_row.start) - overlap
                    raw_id = f"{run_id}|{query_row.peak_id}|{ref_row.peak_id}"
                    rows.append(
                        {
                            "overlap_id": hashlib.sha256(raw_id.encode()).hexdigest()[:24],
                            "run_id": run_id, "peak_id": query_row.peak_id,
                            "reference_region_id": ref_row.peak_id, "overlap_bp": overlap,
                            "query_fraction": overlap / (query_row.end - query_row.start),
                            "reference_fraction": overlap / (ref_row.end - ref_row.start),
                            "jaccard": overlap / union,
                            "query_score": query_row.score, "reference_score": ref_row.score,
                            "query_summit": query_row.summit, "reference_summit": ref_row.summit,
                            "boundary_start_error_bp": query_row.start - ref_row.start,
                            "boundary_end_error_bp": query_row.end - ref_row.end,
                            "width_error_bp": (query_row.end - query_row.start) - (ref_row.end - ref_row.start),
                        }
                    )
                index += 1
    return rows


def qualifying(frame: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Apply the primary any-overlap or reciprocal-overlap sensitivity rule."""
    if threshold <= 0:
        return frame
    return frame[(frame["query_fraction"] >= threshold) & (frame["reference_fraction"] >= threshold)]


def metric_rows(
    run_id: str, mode: str, query: pd.DataFrame, reference: pd.DataFrame,
    overlaps: pd.DataFrame, thresholds: list[float],
) -> list[dict]:
    """Build normalized long-form anchor-relative metrics."""
    output = []
    query_bp = interval_union_bp(query)
    reference_bp = interval_union_bp(reference)
    for threshold in thresholds:
        valid = qualifying(overlaps, threshold)
        matched_query = valid["peak_id"].nunique() if not valid.empty else 0
        matched_reference = valid["reference_region_id"].nunique() if not valid.empty else 0
        intersection_bp = 0
        if not valid.empty:
            intersection_segments = []
            query_lookup = query.set_index("peak_id")
            reference_lookup = reference.set_index("peak_id")
            for row in valid.itertuples():
                qrow = query_lookup.loc[row.peak_id]
                rrow = reference_lookup.loc[row.reference_region_id]
                intersection_segments.append(
                    {"chromosome": qrow.chromosome, "start": max(qrow.start, rrow.start), "end": min(qrow.end, rrow.end)}
                )
            intersection_bp = interval_union_bp(pd.DataFrame(intersection_segments))
        union_bp = query_bp + reference_bp - intersection_bp
        best = valid.sort_values("jaccard").drop_duplicates("peak_id", keep="last") if not valid.empty else valid
        correlation = np.nan
        if len(best) >= 2 and best["query_score"].notna().all() and best["reference_score"].notna().all():
            correlation = float(spearmanr(best["query_score"], best["reference_score"]).statistic)
        metrics = {
            "anchor_peak_retention": matched_reference / len(reference) if len(reference) else np.nan,
            "query_peak_concordance": matched_query / len(query) if len(query) else np.nan,
            "genomic_base_jaccard": intersection_bp / union_bp if union_bp else np.nan,
            "peak_count_ratio": len(query) / len(reference) if len(reference) else np.nan,
            "called_base_ratio": query_bp / reference_bp if reference_bp else np.nan,
            "score_rank_correlation": correlation,
        }
        if mode == "narrow":
            distances = (best["query_summit"] - best["reference_summit"]).abs() if not best.empty else pd.Series(dtype=float)
            metrics["median_summit_displacement_bp"] = float(distances.median()) if distances.notna().any() else np.nan
        else:
            metrics["median_boundary_displacement_bp"] = (
                float(pd.concat([best["boundary_start_error_bp"].abs(), best["boundary_end_error_bp"].abs()]).median())
                if not best.empty else np.nan
            )
            metrics["median_width_displacement_bp"] = float(best["width_error_bp"].abs().median()) if not best.empty else np.nan
        status = "valid_zero_peak" if len(query) == 0 else "complete"
        for metric, value in metrics.items():
            output.append({"run_id": run_id, "metric": metric, "threshold": threshold, "value": value, "status": status})
    return output


def main() -> None:
    """Write normalized peaks, overlaps, recovery rows, and run metrics."""
    args = parse_args()
    thresholds = sorted(set(args.threshold or [0.0, 0.1, 0.5]))
    query = read_peaks(args.query_peaks, args.run_id, args.mode)
    reference = read_peaks(args.anchor_peaks, f"{args.study_id}__reference", args.mode)
    overlaps = pd.DataFrame(overlap_rows(query, reference, args.run_id))
    if overlaps.empty:
        overlaps = pd.DataFrame(columns=[
            "overlap_id", "run_id", "peak_id", "reference_region_id", "overlap_bp",
            "query_fraction", "reference_fraction", "jaccard", "query_score", "reference_score",
            "query_summit", "reference_summit", "boundary_start_error_bp", "boundary_end_error_bp", "width_error_bp",
        ])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    query.assign(run_id=args.run_id)[["peak_id", "run_id", "chromosome", "start", "end", "score", "summit", "width"]].to_csv(
        args.output_dir / "peaks.csv", index=False
    )
    reference.assign(study_id=args.study_id, reference_run_id=f"{args.study_id}__full_control_anchor", ontology_class="")\
        .rename(columns={"peak_id": "reference_region_id"})[[
            "reference_region_id", "study_id", "reference_run_id", "chromosome", "start", "end", "ontology_class", "score", "summit"
        ]].to_csv(args.output_dir / "reference_regions.csv", index=False)
    overlap_exports = []
    recovery_exports = []
    for threshold in thresholds:
        valid = qualifying(overlaps, threshold).copy()
        valid["threshold"] = threshold
        valid["overlap_id"] = valid["overlap_id"].map(
            lambda value: hashlib.sha256(f"{value}|{threshold}".encode()).hexdigest()[:24]
        )
        overlap_exports.append(valid[["overlap_id", "run_id", "peak_id", "reference_region_id", "overlap_bp", "query_fraction", "reference_fraction", "jaccard", "threshold"]])
        recovered = set(valid["reference_region_id"])
        recovery_exports.extend(
            {"run_id": args.run_id, "reference_region_id": region_id, "threshold": threshold, "recovered": int(region_id in recovered)}
            for region_id in reference["peak_id"]
        )
    pd.concat(overlap_exports, ignore_index=True).to_csv(args.output_dir / "peak_overlaps.csv", index=False)
    pd.DataFrame(recovery_exports, columns=["run_id", "reference_region_id", "threshold", "recovered"]).to_csv(
        args.output_dir / "reference_region_recovery.csv", index=False
    )
    pd.DataFrame(metric_rows(args.run_id, args.mode, query, reference, overlaps, thresholds)).to_csv(
        args.output_dir / "run_metrics.csv", index=False
    )


main()
