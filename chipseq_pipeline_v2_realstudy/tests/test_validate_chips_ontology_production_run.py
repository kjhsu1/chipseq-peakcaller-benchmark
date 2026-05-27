"""Tests for production validation of ChIPs ontology outputs."""

"""Imports"""

from pathlib import Path
import json
import subprocess
import sys

import pandas as pd


"""Functions"""


def write_csv(path: Path, rows: list[dict]) -> None:
    """Write a list of dictionaries to a CSV path."""
    pd.DataFrame(rows).to_csv(path, index=False)


def test_validator_passes_on_complete_minimal_output(tmp_path: Path):
    """The production validator should pass on a complete minimal output root."""
    repo_root = Path(__file__).resolve().parents[1]
    run_table = tmp_path / "prototype_run_table.csv"
    output_root = tmp_path / "output_root"
    summary_dir = output_root / "analysis_outputs/chips_ontology/summary"
    repro_root = output_root / "reproducibility"
    run_id = "run_a"

    write_csv(
        run_table,
        [
            {
                "run_id": run_id,
                "study_id": "study_a",
                "coverage_treat": 10,
                "coverage_ctrl": 4,
                "seed": 11,
                "fragment_length": 150,
                "read_length": 38,
                "aligner": "bowtie2",
                "peakcaller": "macs2",
                "macs2_mode": "broad",
                "assembly": "ce11",
            }
        ],
    )

    (output_root / "RUN_COMPLETE").parent.mkdir(parents=True, exist_ok=True)
    (output_root / "RUN_COMPLETE").write_text("", encoding="utf-8")

    for rel_path in [
        f"results_chips/{run_id}/bowtie2/treat/aligned.sorted.bam",
        f"results_chips/{run_id}/bowtie2/treat/aligned.sorted.bam.bai",
        f"results_chips/{run_id}/bowtie2/con/aligned.sorted.bam",
        f"results_chips/{run_id}/bowtie2/con/aligned.sorted.bam.bai",
    ]:
        path = output_root / rel_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("bam", encoding="utf-8")

    peak_bed = output_root / f"results_chips/{run_id}/peaks/macs2/{run_id}_peaks.bed"
    peak_bed.parent.mkdir(parents=True, exist_ok=True)
    peak_bed.write_text("chrI\t0\t10\n", encoding="utf-8")

    region_rows = [
        {
            "run_id": run_id,
            "study_id": "study_a",
            "coverage_ctrl": 4,
            "precision": 1.0,
            "recall": 1.0,
            "f1": 1.0,
            "failure_mode": "none",
            "ontology_class": "flat_background__strong_peak",
        }
    ]
    for name in ["region_metrics.csv", "classified.csv"]:
        path = output_root / f"analysis_outputs/chips_ontology/{run_id}/{name}"
        path.parent.mkdir(parents=True, exist_ok=True)
        write_csv(path, region_rows)

    combined_path = output_root / "analysis_outputs/chips_ontology/combined_region_metrics.csv"
    combined_path.parent.mkdir(parents=True, exist_ok=True)
    combined_rows = []
    for class_name in ["class_a", "class_b", "class_c"]:
        for index in range(10):
            combined_rows.append(
                {
                    "run_id": run_id,
                    "study_id": "study_a",
                    "coverage_ctrl": index + 1,
                    "precision": 0.8,
                    "recall": 0.8,
                    "f1": 0.8,
                    "failure_mode": "none",
                    "ontology_class": class_name,
                }
            )
    write_csv(combined_path, combined_rows)

    summary_dir.mkdir(parents=True, exist_ok=True)
    for name in [
        "per_run_overall_metrics.csv",
        "per_ontology_metrics.csv",
        "control_response_by_ontology.csv",
        "failure_mode_metrics.csv",
    ]:
        write_csv(summary_dir / name, combined_rows[:3])
    write_csv(
        summary_dir / "ontology_class_coverage.csv",
        [
            {"ontology_class": "class_a", "run_count": 10, "region_count": 10, "coverage_threshold_runs": 10, "meets_min_runs": True},
            {"ontology_class": "class_b", "run_count": 10, "region_count": 10, "coverage_threshold_runs": 10, "meets_min_runs": True},
            {"ontology_class": "class_c", "run_count": 10, "region_count": 10, "coverage_threshold_runs": 10, "meets_min_runs": True},
        ],
    )
    for name in ["ontology_f1_heatmap.png", "control_response_by_ontology.png"]:
        (summary_dir / name).write_bytes(
            b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01\x08\x02\x00\x00\x00\x90wS\xde\x00\x00\x00\x0cIDAT\x08\xd7c\xf8\xff\xff?\x00\x05\xfe\x02\xfeA\xa4\x9d\xb1\x00\x00\x00\x00IEND\xaeB`\x82"
        )

    repro_root.mkdir(parents=True, exist_ok=True)
    for rel_path in [
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
    ]:
        path = repro_root / rel_path
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix == ".json":
            path.write_text(json.dumps({"ok": True}), encoding="utf-8")
        else:
            path.write_text("ok\n", encoding="utf-8")

    report_path = repro_root / "validation_report.md"
    subprocess.run(
        [
            sys.executable,
            "scripts/validate_chips_ontology_production_run.py",
            "--run-table",
            str(run_table),
            "--output-root",
            str(output_root),
            "--summary-dir",
            str(summary_dir),
            "--write-report",
            str(report_path),
        ],
        check=True,
        cwd=repo_root,
    )

    assert report_path.exists()
    assert (repro_root / "validation_summary.json").exists()
