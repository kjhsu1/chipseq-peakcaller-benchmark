"""Annotate anchor regions by empirical ce11 genomic context and summarize recovery."""

"""Imports"""

import argparse
from pathlib import Path

import pandas as pd


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse genomic-context annotation arguments."""
    parser = argparse.ArgumentParser(description="Annotate Realstudy v2 anchor regions with ce11 genomic classes.")
    parser.add_argument("--reference-regions", type=Path, required=True)
    parser.add_argument("--recovery", type=Path, required=True)
    parser.add_argument("--run-table", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--annotated-regions", type=Path, required=True)
    parser.add_argument("--stratified-metrics", type=Path, required=True)
    parser.add_argument("--promoter-bp", type=int, default=1000)
    return parser.parse_args()


def read_gff(path: Path, promoter_bp: int) -> dict[str, dict[str, list[tuple[int, int]]]]:
    """Index exon, gene-body, and strand-aware promoter intervals by chromosome."""
    indexed = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 9 or fields[2] not in {"gene", "exon"}:
                continue
            chromosome, feature = fields[0], fields[2]
            start, end, strand = int(fields[3]) - 1, int(fields[4]), fields[6]
            bucket = indexed.setdefault(chromosome, {"promoter": [], "exon": [], "gene_body": []})
            if feature == "exon":
                bucket["exon"].append((start, end))
            else:
                bucket["gene_body"].append((start, end))
                tss = start if strand != "-" else end
                bucket["promoter"].append((max(0, tss - promoter_bp), tss + promoter_bp))
    for chromosome in indexed:
        for feature in indexed[chromosome]:
            indexed[chromosome][feature].sort()
    return indexed


def overlaps(intervals: list[tuple[int, int]], start: int, end: int) -> bool:
    """Return whether one interval overlaps any sorted annotation interval."""
    for feature_start, feature_end in intervals:
        if feature_start >= end:
            return False
        if feature_end > start:
            return True
    return False


def genomic_class(index: dict, chromosome: str, start: int, end: int) -> str:
    """Assign promoter, exon, gene-body, or intergenic by fixed priority."""
    bucket = index.get(chromosome, {})
    for label in ["promoter", "exon", "gene_body"]:
        if overlaps(bucket.get(label, []), start, end):
            return label
    return "intergenic"


def main() -> None:
    """Write annotated anchor regions and long-form stratified recovery metrics."""
    args = parse_args()
    regions = pd.read_csv(args.reference_regions).drop_duplicates("reference_region_id")
    index = read_gff(args.gff, args.promoter_bp)
    regions["ontology_class"] = [
        genomic_class(index, row.chromosome, int(row.start), int(row.end)) for row in regions.itertuples()
    ]
    recovery = pd.read_csv(args.recovery)
    runs = pd.read_csv(args.run_table, dtype=str).fillna("")
    joined = recovery.merge(regions[["reference_region_id", "ontology_class"]], on="reference_region_id", validate="many_to_one")
    joined = joined.merge(runs[["run_id", "study_id", "control_coverage_x"]], on="run_id", validate="many_to_one")
    summary = joined.groupby(
        ["run_id", "study_id", "control_coverage_x", "ontology_class", "threshold"], as_index=False
    )["recovered"].mean()
    summary = summary.rename(columns={"ontology_class": "stratum_value", "recovered": "value"})
    summary["stratum_type"] = "ce11_genomic_context"
    summary["metric"] = "anchor_peak_retention"
    args.annotated_regions.parent.mkdir(parents=True, exist_ok=True)
    regions.to_csv(args.annotated_regions, index=False)
    summary[["run_id", "study_id", "control_coverage_x", "stratum_type", "stratum_value", "metric", "threshold", "value"]].to_csv(
        args.stratified_metrics, index=False
    )


main()
