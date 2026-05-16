"""Tests for ontology evaluation aggregation behavior."""

"""Imports"""

from pathlib import Path

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
