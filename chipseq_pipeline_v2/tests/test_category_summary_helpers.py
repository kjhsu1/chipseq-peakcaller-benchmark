"""Tests for controlled-category summary helpers."""

"""Imports"""

from pathlib import Path

import pandas as pd

from scripts.category_summary_lib import (
    best_control_by_mean_f1,
    canonical_label,
    parse_simple_yaml_map,
    precision_dropoff_points,
    recall_saturation_points,
    zero_metric_row_count,
)


"""Functions"""


def test_parse_simple_yaml_map(tmp_path):
    """Simple YAML key-value parsing should preserve mappings."""
    path = tmp_path / "map.yaml"
    path.write_text("a: alpha\nb: beta\n", encoding="utf-8")
    assert parse_simple_yaml_map(path) == {"a": "alpha", "b": "beta"}


def test_canonical_label_uses_mapping():
    """Configured labels should override raw config names."""
    assert canonical_label("raw_name", {"raw_name": "clean_name"}) == "clean_name"


def test_zero_metric_row_count_flags_any_zero_metric():
    """Rows with zero-valued core metrics should be counted."""
    df = pd.DataFrame(
        [
            {"total_called": 10, "tp_called": 9, "tp_planted": 9, "precision": 0.9, "recall": 0.9, "f1": 0.9},
            {"total_called": 0, "tp_called": 0, "tp_planted": 0, "precision": 0.0, "recall": 0.0, "f1": 0.0},
        ]
    )
    assert zero_metric_row_count(df) == 1


def test_best_control_by_mean_f1_prefers_highest_mean():
    """Best control depth should come from grouped mean F1."""
    df = pd.DataFrame(
        [
            {"coverage_ctrl": 1.0, "f1": 0.4},
            {"coverage_ctrl": 1.0, "f1": 0.6},
            {"coverage_ctrl": 2.0, "f1": 0.9},
        ]
    )
    assert best_control_by_mean_f1(df) == 2.0


def test_recall_saturation_and_precision_dropoff_points():
    """Threshold helpers should return first control-depth crossings."""
    plot_df = pd.DataFrame(
        [
            {"coverage_treat": 5, "coverage_ctrl": 0.5, "precision": 1.0, "recall": 0.5},
            {"coverage_treat": 5, "coverage_ctrl": 1.0, "precision": 0.9, "recall": 0.98},
            {"coverage_treat": 5, "coverage_ctrl": 2.0, "precision": 0.7, "recall": 1.0},
        ]
    )
    saturation = recall_saturation_points(plot_df, fraction_of_max=0.98)
    dropoff = precision_dropoff_points(plot_df, retained_fraction=0.8)
    assert saturation.iloc[0]["recall_saturation_ctrl"] == 1.0
    assert dropoff.iloc[0]["precision_dropoff_ctrl"] == 2.0
