"""
Curate archived FP/FN example runs into one visualization pack.
"""

"""Imports"""

import argparse
import csv
from datetime import datetime
from pathlib import Path
import shutil
from typing import Dict, List, Tuple


"""Constants"""

DEFAULT_INVESTIGATION_DIR = Path("analysis_outputs/peak_recovery_investigation_20260422")
DEFAULT_SELECTIONS = [
    (
        "balanced_flatearth_peaks_broad_integrated_288",
        "0119",
        "recall_seed37_perfect_precision_flatearth_broad_run0119",
        "Perfect precision, poor recall, seed 37; inspect repeated missed planted centers.",
        ["I:24022", "I:106027", "II:10531", "III:110778", "IV:11262"],
    ),
    (
        "balanced_realistic_peaks_wavy_narrow_integrated_288",
        "0119",
        "recall_seed37_perfect_precision_wavy_narrow_run0119",
        "Perfect precision, poor recall, seed 37; inspect repeated missed planted centers.",
        ["I:24022", "I:106027", "II:10531", "III:110778", "IV:11262"],
    ),
    (
        "balanced_flatearth_peaks_broad_integrated_288",
        "0071",
        "fp_hotspot_iv_9p6kb_seed23_flatearth_broad_run0071",
        "False-positive hotspot example, seed 23; inspect IV:9620-9773.",
        ["IV:9620-9773"],
    ),
    (
        "balanced_realistic_peaks_wavy_narrow_integrated_288",
        "0071",
        "fp_hotspot_iv_9p6kb_seed23_wavy_narrow_run0071",
        "False-positive hotspot example, seed 23; inspect IV:9620-9773.",
        ["IV:9620-9773"],
    ),
    (
        "balanced_flatearth_peaks_broad_integrated_288",
        "0167",
        "fp_hotspot_iii_100p5kb_seed53_flatearth_broad_run0167",
        "False-positive hotspot example, seed 53; inspect III:100397-100726.",
        ["III:100397-100726"],
    ),
    (
        "balanced_realistic_peaks_wavy_narrow_integrated_288",
        "0167",
        "fp_hotspot_iii_100p5kb_seed53_wavy_narrow_run0167",
        "False-positive hotspot example, seed 53; inspect III:100397-100726.",
        ["III:100397-100726"],
    ),
]


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Copy selected archived FP/FN runs into one curated visualization pack."
    )
    parser.add_argument(
        "--investigation-dir",
        type=Path,
        default=DEFAULT_INVESTIGATION_DIR,
        help="Directory containing peak recovery investigation CSV outputs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Defaults to analysis_outputs/curated_fp_fn_visual_pack_TIMESTAMP.",
    )
    return parser.parse_args()


def resolve_output_dir(output_dir: Path) -> Path:
    """Return explicit output directory or a timestamped default."""
    if output_dir is not None:
        return output_dir
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return Path("analysis_outputs") / f"curated_fp_fn_visual_pack_{timestamp}"


def read_csv(path: Path) -> List[Dict[str, str]]:
    """Read CSV rows as dictionaries."""
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def copy_required(src: Path, dst: Path) -> None:
    """Copy one required file."""
    if not src.exists():
        raise FileNotFoundError(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def maybe_copy(src: Path, dst: Path) -> bool:
    """Copy one optional file if present."""
    if not src.exists():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return True


def run_lookup(rows: List[Dict[str, str]]) -> Dict[Tuple[str, str], Dict[str, str]]:
    """Build category/run lookup."""
    return {(row["category"], row["run_id"]): row for row in rows}


def filter_rows(rows: List[Dict[str, str]], category: str, run_id: str) -> List[Dict[str, str]]:
    """Filter detail rows by category and run."""
    return [row for row in rows if row["category"] == category and row["run_id"] == run_id]


def peak_path_for(row: Dict[str, str]) -> Path:
    """Return copied peak path from per-run summary row."""
    return Path(row["called_peak_path"])


def source_run_dir(row: Dict[str, str]) -> Path:
    """Return source archived run directory."""
    return Path(row["source_dir"]) / row["run_id"]


def format_param_lines(row: Dict[str, str]) -> List[str]:
    """Return parameter summary lines."""
    keys = [
        "category",
        "run_id",
        "source_dir",
        "aligner",
        "peakcaller",
        "macs2_mode",
        "seed",
        "tf_seed",
        "map_seed",
        "coverage_treat",
        "coverage_ctrl",
        "tf_enrich",
        "tf_sigma",
        "tf_peak_count_treat",
        "tf_exp",
        "gc_exp",
        "acc_exp",
        "map_coverage_pct",
        "map_sigma",
        "map_enrich",
        "map_exp",
        "precision",
        "recall",
        "f1",
        "false_positive_called",
        "missed_planted",
    ]
    return [f"{key}: {row.get(key, '')}" for key in keys]


def write_run_info(
    path: Path,
    row: Dict[str, str],
    selection_reason: str,
    loci: List[str],
    missed_rows: List[Dict[str, str]],
    fp_rows: List[Dict[str, str]],
    copied_pmf: bool,
) -> None:
    """Write one per-run text summary."""
    lines = [
        "Curated peak recovery visualization run",
        "",
        f"selection_reason: {selection_reason}",
        f"suggested_loci: {', '.join(loci)}",
        f"pmf_csv_copied: {copied_pmf}",
        "",
        "Parameters and metrics:",
        *format_param_lines(row),
        "",
        "Missed planted centers for this run:",
    ]
    if missed_rows:
        lines.extend(
            f"- {detail['planted_chrom']}:{detail['planted_center']}"
            for detail in missed_rows
        )
    else:
        lines.append("- none")
    lines.extend(["", "False-positive called intervals for this run:"])
    if fp_rows:
        lines.extend(
            f"- {detail['called_chrom']}:{detail['called_start']}-{detail['called_end']}"
            for detail in fp_rows
        )
    else:
        lines.append("- none")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_root_readme(path: Path, manifest_rows: List[Dict[str, str]]) -> None:
    """Write top-level README for the curated pack."""
    lines = [
        "# Curated FP/FN Visualization Pack",
        "",
        "This folder reuses archived balanced-288 sweep outputs. No runs were rerun.",
        "",
        "BigWig files were not generated in this pack because the current environment does not provide a BigWig conversion tool (`bamCoverage`, `bedGraphToBigWig`, or `pyBigWig`).",
        "Each run folder instead includes coordinate-sorted BAM/BAI pairs for treatment and control, which can be opened directly in IGV.",
        "",
        "## Included Runs",
        "",
        "| label | category | run_id | seed | precision | recall | reason | loci |",
        "| --- | --- | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    for row in manifest_rows:
        lines.append(
            f"| {row['label']} | {row['category']} | {row['run_id']} | {row['seed']} | "
            f"{row['precision']} | {row['recall']} | {row['selection_reason']} | {row['suggested_loci']} |"
        )
    lines.extend(
        [
            "",
            "## Contents Per Run",
            "",
            "- `treat/aligned.sorted.bam` and `.bai`",
            "- `control/aligned.sorted.bam` and `.bai`",
            "- `planted_peaks.bed`",
            "- `called_peaks.bed`",
            "- `run_info.txt`",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def curate_one_run(
    output_dir: Path,
    row: Dict[str, str],
    selection_label: str,
    selection_reason: str,
    loci: List[str],
    missed_rows: List[Dict[str, str]],
    fp_rows: List[Dict[str, str]],
) -> Dict[str, str]:
    """Copy one selected run into the curated output."""
    src_run_dir = source_run_dir(row)
    run_dir = output_dir / selection_label

    copy_required(
        src_run_dir / row["aligner"] / "treat" / "aligned.sorted.bam",
        run_dir / "treat" / "aligned.sorted.bam",
    )
    copy_required(
        src_run_dir / row["aligner"] / "treat" / "aligned.sorted.bam.bai",
        run_dir / "treat" / "aligned.sorted.bam.bai",
    )
    copy_required(
        src_run_dir / row["aligner"] / "con" / "aligned.sorted.bam",
        run_dir / "control" / "aligned.sorted.bam",
    )
    copy_required(
        src_run_dir / row["aligner"] / "con" / "aligned.sorted.bam.bai",
        run_dir / "control" / "aligned.sorted.bam.bai",
    )
    copy_required(src_run_dir / "treat" / "planted_peaks.bed", run_dir / "planted_peaks.bed")
    copy_required(peak_path_for(row), run_dir / "called_peaks.bed")

    copied_pmf = maybe_copy(src_run_dir / "treat" / "pmf.csv", run_dir / "treat" / "pmf.csv")
    maybe_copy(src_run_dir / "con" / "pmf.csv", run_dir / "control" / "pmf.csv")

    write_run_info(
        run_dir / "run_info.txt",
        row,
        selection_reason,
        loci,
        missed_rows,
        fp_rows,
        copied_pmf,
    )

    return {
        "label": selection_label,
        "category": row["category"],
        "run_id": row["run_id"],
        "seed": row["seed"],
        "precision": row["precision"],
        "recall": row["recall"],
        "selection_reason": selection_reason,
        "suggested_loci": ", ".join(loci),
        "source_dir": row["source_dir"],
    }


def main() -> None:
    """Curate the selected archived runs."""
    args = parse_args()
    investigation_dir = args.investigation_dir
    output_dir = resolve_output_dir(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    per_run_rows = read_csv(investigation_dir / "per_run_recovery_summary.csv")
    missed_rows = read_csv(investigation_dir / "missed_planted_peaks.csv")
    fp_rows = read_csv(investigation_dir / "false_positive_called_peaks.csv")
    row_by_key = run_lookup(per_run_rows)

    manifest_rows = []
    for category, run_id, label, reason, loci in DEFAULT_SELECTIONS:
        row = row_by_key[(category, run_id)]
        manifest_rows.append(
            curate_one_run(
                output_dir,
                row,
                label,
                reason,
                loci,
                filter_rows(missed_rows, category, run_id),
                filter_rows(fp_rows, category, run_id),
            )
        )

    manifest_path = output_dir / "pilot_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "label",
                "category",
                "run_id",
                "seed",
                "precision",
                "recall",
                "selection_reason",
                "suggested_loci",
                "source_dir",
            ],
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    write_root_readme(output_dir / "README.md", manifest_rows)


main()
