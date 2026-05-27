"""Validate a production ChIPs ontology output root."""

"""Imports"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


"""Functions"""


REQUIRED_SUMMARY_CSVS = [
    "per_run_overall_metrics.csv",
    "per_ontology_metrics.csv",
    "control_response_by_ontology.csv",
    "failure_mode_metrics.csv",
]
REQUIRED_SUMMARY_PNGS = [
    "ontology_f1_heatmap.png",
    "control_response_by_ontology.png",
]
REPRO_REQUIRED_FILES = [
    "config.yaml",
    "configs/chips_cluster_full.yaml",
    "metadata/prototype_run_table.csv",
    "metadata/prototype_run_table.summary.json",
    "slurm/submit_chips_realsim_ontology_128cpu_2tb.sh",
    "slurm/chips_realsim_ontology_128cpu_2tb.sbatch",
    "commit_hash.txt",
    "output_root.txt",
    "log_file.txt",
    "validation_report.md",
    "validation_summary.json",
    "runtime_resource_report.txt",
    "runtime_resource_report.csv",
    "top_conclusions.md",
]
COVERAGE_CLASS_MIN_RUNS = 10
COVERAGE_MIN_CLASSES = 3


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Validate a production ChIPs ontology output root.")
    parser.add_argument("--run-table", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--summary-dir", type=Path)
    parser.add_argument("--write-report", type=Path, required=True)
    return parser.parse_args()


def load_expected_runs(run_table: Path) -> list[str]:
    """Return sorted unique expected run ids from the run table."""
    run_df = pd.read_csv(run_table)
    return sorted(str(value) for value in run_df["run_id"].dropna().astype(str).unique())


def readable_csv(path: Path) -> tuple[bool, int, list[str], pd.DataFrame | None]:
    """Return readability, row count, columns, and the parsed dataframe."""
    try:
        df = pd.read_csv(path)
    except Exception:
        return False, 0, [], None
    return True, len(df.index), list(df.columns), df


def readable_png(path: Path) -> bool:
    """Return whether a PNG file can be loaded."""
    try:
        plt.imread(path)
    except Exception:
        return False
    return True


def add_result(results: list[dict], criterion: int, title: str, passed: bool, evidence: list[str], details: list[str]) -> None:
    """Append one criterion result row."""
    results.append(
        {
            "criterion": criterion,
            "title": title,
            "status": "PASS" if passed else "FAIL",
            "evidence": evidence,
            "details": details,
        }
    )


def file_list_for_run(output_root: Path, run_id: str) -> dict[str, Path]:
    """Return the per-run required file paths for one run id."""
    return {
        "treat_bam": output_root / f"results_chips/{run_id}/bowtie2/treat/aligned.sorted.bam",
        "treat_bai": output_root / f"results_chips/{run_id}/bowtie2/treat/aligned.sorted.bam.bai",
        "con_bam": output_root / f"results_chips/{run_id}/bowtie2/con/aligned.sorted.bam",
        "con_bai": output_root / f"results_chips/{run_id}/bowtie2/con/aligned.sorted.bam.bai",
        "peak_bed": output_root / f"results_chips/{run_id}/peaks/macs2/{run_id}_peaks.bed",
        "region_metrics": output_root / f"analysis_outputs/chips_ontology/{run_id}/region_metrics.csv",
        "classified": output_root / f"analysis_outputs/chips_ontology/{run_id}/classified.csv",
    }


def build_report_text(
    output_root: Path,
    expected_runs: list[str],
    completed_runs: list[str],
    results: list[dict],
    summary_payload: dict,
) -> str:
    """Render the Markdown validation report."""
    lines = [
        "# ChIPs Ontology Production Validation Report",
        "",
        f"- output root: `{output_root}`",
        f"- expected run count: `{len(expected_runs)}`",
        f"- completed run count: `{len(completed_runs)}`",
        "",
        "## Criterion Status",
    ]
    for result in results:
        lines.append(
            f"- {result['criterion']}. {result['title']}: `{result['status']}`"
        )
        for detail in result["details"]:
            lines.append(f"  - {detail}")
        for evidence in result["evidence"]:
            lines.append(f"  - evidence: `{evidence}`")
    lines.extend(
        [
            "",
            "## Summary",
            f"- passed: `{summary_payload['passed_count']}`",
            f"- failed: `{summary_payload['failed_count']}`",
            f"- overall status: `{summary_payload['overall_status']}`",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    """Validate the production output root and emit Markdown plus JSON reports."""
    args = parse_args()
    output_root = args.output_root.resolve()
    summary_dir = (args.summary_dir or output_root / "analysis_outputs/chips_ontology/summary").resolve()
    repro_root = output_root / "reproducibility"
    json_path = args.write_report.with_name("validation_summary.json")
    expected_runs = load_expected_runs(args.run_table)

    results: list[dict] = []
    completed_runs: list[str] = []
    missing_alignment_files: list[str] = []
    missing_peak_beds: list[str] = []
    empty_peak_beds: list[str] = []
    missing_region_metrics: list[str] = []
    missing_classified: list[str] = []
    unreadable_region_metrics: list[str] = []
    unreadable_classified: list[str] = []

    run_complete = output_root / "RUN_COMPLETE"
    run_failed = output_root / "RUN_FAILED"
    add_result(
        results,
        1,
        "Successful end-to-end realstudy production run",
        run_complete.exists() and not run_failed.exists(),
        [str(run_complete), str(run_failed)],
        [
            f"RUN_COMPLETE exists: {run_complete.exists()}",
            f"RUN_FAILED exists: {run_failed.exists()}",
        ],
    )

    for run_id in expected_runs:
        paths = file_list_for_run(output_root, run_id)
        run_missing_alignment = [str(path) for name, path in paths.items() if name.endswith(("bam", "bai")) and not path.exists()]
        if run_missing_alignment:
            missing_alignment_files.extend(run_missing_alignment)
        peak_bed = paths["peak_bed"]
        if not peak_bed.exists():
            missing_peak_beds.append(str(peak_bed))
        elif peak_bed.stat().st_size == 0:
            empty_peak_beds.append(str(peak_bed))

        region_metrics = paths["region_metrics"]
        classified = paths["classified"]
        if not region_metrics.exists():
            missing_region_metrics.append(str(region_metrics))
        else:
            readable, _, _, _ = readable_csv(region_metrics)
            if not readable:
                unreadable_region_metrics.append(str(region_metrics))
        if not classified.exists():
            missing_classified.append(str(classified))
        else:
            readable, _, _, _ = readable_csv(classified)
            if not readable:
                unreadable_classified.append(str(classified))

        if not any(
            [
                run_missing_alignment,
                not peak_bed.exists(),
                peak_bed.exists() and peak_bed.stat().st_size == 0,
                not region_metrics.exists(),
                not classified.exists(),
                str(region_metrics) in unreadable_region_metrics,
                str(classified) in unreadable_classified,
            ]
        ):
            completed_runs.append(run_id)

    missing_run_ids = sorted(set(expected_runs) - set(completed_runs))
    extra_run_ids = sorted(set(completed_runs) - set(expected_runs))
    add_result(
        results,
        2,
        "Complete expected run set produced",
        len(expected_runs) == len(completed_runs) and not missing_run_ids and not extra_run_ids,
        [str(args.run_table)],
        [
            f"expected run count: {len(expected_runs)}",
            f"completed run count: {len(completed_runs)}",
            f"missing run ids: {', '.join(missing_run_ids[:10]) if missing_run_ids else 'none'}",
            f"extra run ids: {', '.join(extra_run_ids[:10]) if extra_run_ids else 'none'}",
        ],
    )
    add_result(
        results,
        3,
        "All simulated alignment outputs present",
        len(missing_alignment_files) == 0,
        [str(output_root / "results_chips")],
        [
            f"missing alignment file count: {len(missing_alignment_files)}",
        ],
    )
    add_result(
        results,
        4,
        "All peak-call outputs present",
        len(missing_peak_beds) == 0 and len(empty_peak_beds) == 0,
        [str(output_root / "results_chips")],
        [
            f"missing peak bed count: {len(missing_peak_beds)}",
            f"empty peak bed count: {len(empty_peak_beds)}",
        ],
    )
    add_result(
        results,
        5,
        "Ontology scoring complete for all runs",
        not any([missing_region_metrics, missing_classified, unreadable_region_metrics, unreadable_classified]),
        [str(output_root / "analysis_outputs/chips_ontology")],
        [
            f"missing region_metrics count: {len(missing_region_metrics)}",
            f"missing classified count: {len(missing_classified)}",
            f"unreadable region_metrics count: {len(unreadable_region_metrics)}",
            f"unreadable classified count: {len(unreadable_classified)}",
        ],
    )

    combined_path = output_root / "analysis_outputs/chips_ontology/combined_region_metrics.csv"
    combined_readable, combined_rows, _, combined_df = readable_csv(combined_path) if combined_path.exists() else (False, 0, [], None)
    combined_run_ids = []
    if combined_df is not None and "run_id" in combined_df.columns:
        combined_run_ids = sorted(str(value) for value in combined_df["run_id"].dropna().astype(str).unique())
    add_result(
        results,
        6,
        "Combined ontology table built",
        combined_readable and combined_rows > 0 and sorted(combined_run_ids) == expected_runs,
        [str(combined_path)],
        [
            f"combined table exists: {combined_path.exists()}",
            f"combined table readable: {combined_readable}",
            f"combined table rows: {combined_rows}",
            f"combined run count: {len(combined_run_ids)}",
        ],
    )

    summary_csv_failures: list[str] = []
    for name in REQUIRED_SUMMARY_CSVS:
        path = summary_dir / name
        readable, rows, _, _ = readable_csv(path) if path.exists() else (False, 0, [], None)
        if not path.exists() or not readable or rows == 0:
            summary_csv_failures.append(f"{path} (exists={path.exists()} readable={readable} rows={rows})")
    add_result(
        results,
        7,
        "Final performance summary tables built",
        len(summary_csv_failures) == 0,
        [str(summary_dir)],
        summary_csv_failures or ["all required summary CSVs exist, are readable, and have data rows"],
    )

    summary_png_failures: list[str] = []
    for name in REQUIRED_SUMMARY_PNGS:
        path = summary_dir / name
        readable = path.exists() and path.stat().st_size > 0 and readable_png(path)
        if not readable:
            size = path.stat().st_size if path.exists() else 0
            summary_png_failures.append(f"{path} (exists={path.exists()} size={size})")
    add_result(
        results,
        8,
        "Final performance plots built",
        len(summary_png_failures) == 0,
        [str(summary_dir)],
        summary_png_failures or ["all required PNGs exist, are non-empty, and are readable"],
    )

    coverage_path = summary_dir / "ontology_class_coverage.csv"
    coverage_readable, coverage_rows, _, coverage_df = readable_csv(coverage_path) if coverage_path.exists() else (False, 0, [], None)
    covered_class_count = 0
    coverage_pass = False
    if coverage_df is not None and "run_count" in coverage_df.columns:
        covered_class_count = int((coverage_df["run_count"] >= COVERAGE_CLASS_MIN_RUNS).sum())
        coverage_pass = covered_class_count >= COVERAGE_MIN_CLASSES
    add_result(
        results,
        9,
        "Coverage of ontology classes is acceptable",
        coverage_readable and coverage_rows > 0 and coverage_pass,
        [str(coverage_path)],
        [
            f"coverage table exists: {coverage_path.exists()}",
            f"coverage table readable: {coverage_readable}",
            f"covered classes >= {COVERAGE_CLASS_MIN_RUNS} runs: {covered_class_count}",
            f"threshold requires at least {COVERAGE_MIN_CLASSES} classes",
        ],
    )

    repro_failures = []
    for rel_path in REPRO_REQUIRED_FILES:
        path = repro_root / rel_path
        if path == args.write_report or path == json_path:
            continue
        if not path.exists():
            repro_failures.append(str(path))
    add_result(
        results,
        11,
        "Reproducibility package complete",
        len(repro_failures) == 0,
        [str(repro_root)],
        repro_failures or ["all required reproducibility files are present"],
    )

    failed_count = sum(1 for result in results if result["status"] == "FAIL")
    passed_count = sum(1 for result in results if result["status"] == "PASS")
    summary_payload = {
        "output_root": str(output_root),
        "expected_run_count": len(expected_runs),
        "completed_run_count": len(completed_runs),
        "completed_run_ids": completed_runs,
        "criteria": results,
        "passed_count": passed_count,
        "failed_count": failed_count,
        "overall_status": "PASS" if failed_count == 0 else "FAIL",
    }

    args.write_report.parent.mkdir(parents=True, exist_ok=True)
    report_text = build_report_text(output_root, expected_runs, completed_runs, results, summary_payload)
    args.write_report.write_text(report_text, encoding="utf-8")
    json_path.write_text(json.dumps(summary_payload, indent=2), encoding="utf-8")

    raise SystemExit(0 if failed_count == 0 else 1)


main()
