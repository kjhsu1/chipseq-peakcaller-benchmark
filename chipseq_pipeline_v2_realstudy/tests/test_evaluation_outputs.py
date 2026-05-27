"""Tests for ontology evaluation aggregation behavior."""

"""Imports"""

from pathlib import Path
import subprocess
import sys

import pandas as pd


def test_control_response_shape_is_groupable(tmp_path: Path):
    """Prototype evaluation inputs should group cleanly by ontology and control depth."""
    df = pd.DataFrame(
        [
            {
                "run_id": "r1",
                "study_id": "s1",
                "ontology_class": "flat_background__weak_peak",
                "failure_mode": "none",
                "coverage_ctrl": 1,
                "precision": 0.5,
                "recall": 0.4,
                "f1": 0.44,
            },
            {
                "run_id": "r2",
                "study_id": "s1",
                "ontology_class": "flat_background__weak_peak",
                "failure_mode": "none",
                "coverage_ctrl": 4,
                "precision": 0.6,
                "recall": 0.5,
                "f1": 0.54,
            },
        ]
    )
    grouped = df.groupby(["ontology_class", "coverage_ctrl"], as_index=False)[["precision", "recall", "f1"]].mean()
    assert len(grouped) == 2


def test_evaluation_script_writes_required_outputs(tmp_path: Path):
    """The ontology evaluator should emit the required tables and PNGs."""
    input_csv = tmp_path / "combined_region_metrics.csv"
    output_dir = tmp_path / "summary"
    pd.DataFrame(
        [
            {
                "run_id": "r1",
                "study_id": "s1",
                "ontology_class": "flat_background__weak_peak",
                "failure_mode": "none",
                "coverage_ctrl": 1,
                "precision": 0.5,
                "recall": 0.4,
                "f1": 0.44,
            },
            {
                "run_id": "r2",
                "study_id": "s1",
                "ontology_class": "flat_background__weak_peak",
                "failure_mode": "none",
                "coverage_ctrl": 4,
                "precision": 0.6,
                "recall": 0.5,
                "f1": 0.54,
            },
            {
                "run_id": "r3",
                "study_id": "s2",
                "ontology_class": "bumpy_background__strong_peak",
                "failure_mode": "background_dominant",
                "coverage_ctrl": 4,
                "precision": 0.7,
                "recall": 0.8,
                "f1": 0.75,
            },
        ]
    ).to_csv(input_csv, index=False)

    subprocess.run(
        [
            sys.executable,
            "scripts/evaluate_by_region_ontology.py",
            "--input-csv",
            str(input_csv),
            "--output-dir",
            str(output_dir),
        ],
        check=True,
        cwd=Path(__file__).resolve().parents[1],
    )

    assert (output_dir / "control_response_by_ontology.csv").exists()
    assert (output_dir / "control_response_by_ontology.png").exists()
    assert (output_dir / "ontology_class_coverage.csv").exists()
    assert (output_dir / "ontology_f1_heatmap.png").exists()
