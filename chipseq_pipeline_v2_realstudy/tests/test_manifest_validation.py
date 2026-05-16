"""Tests for real-study manifest helpers."""

"""Imports"""

from pathlib import Path

import pandas as pd

from scripts.realstudy_manifest_lib import (
    load_study_selection,
    parse_simple_yaml_list,
    validate_selection_table,
)


"""Functions"""


def test_parse_simple_yaml_list(tmp_path):
    """The lightweight YAML parser should recover study records."""
    path = tmp_path / "study.yaml"
    path.write_text("studies:\n  - study_id: test\n    has_treatment: true\n", encoding="utf-8")
    studies = parse_simple_yaml_list(path)
    assert studies[0]["study_id"] == "test"


def test_validate_selection_table_adds_pairing_flags():
    """Validation should expose explicit eligibility flags."""
    df = pd.DataFrame(
        [
            {"has_treatment": "true", "has_control": "true", "processed_bam_available": "true", "usable_for_realism_scorecard": "true", "usable_for_realstudy_simulation": "true"}
        ]
    )
    validated = validate_selection_table(df)
    assert bool(validated.iloc[0]["eligible_pairing"])
    assert bool(validated.iloc[0]["eligible_files"])


def test_load_study_selection_fills_expected_columns(tmp_path):
    """Selection loader should emit the canonical column set."""
    path = tmp_path / "study.yaml"
    path.write_text("studies:\n  - study_id: test\n    source: ENCODE\n", encoding="utf-8")
    loaded = load_study_selection(path)
    assert "selection_status" in loaded.columns
    assert loaded.iloc[0]["study_id"] == "test"
