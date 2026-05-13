"""Tests for BAM-only realism metric helpers."""

"""Imports"""

import numpy as np
import pandas as pd

from scripts.realism_metrics_lib import (
    gini_coefficient,
    local_bumpiness,
    metaprofile_metrics,
    score_distance,
)


"""Functions"""


def test_gini_coefficient_zero_for_uniform_values():
    """Uniform vectors should have zero Gini."""
    assert gini_coefficient(np.array([3.0, 3.0, 3.0])) == 0.0


def test_local_bumpiness_detects_adjacent_changes():
    """Bumpier vectors should score above zero."""
    assert local_bumpiness(np.array([1.0, 3.0, 1.0])) > 0.0


def test_metaprofile_metrics_capture_shape_summary():
    """Metaprofile helper should expose core summary fields."""
    stats = metaprofile_metrics(np.array([1.0, 2.0, 4.0, 2.0, 1.0]))
    assert stats["summit_height"] == 4.0
    assert stats["half_max_width"] >= 1.0
    assert stats["auc"] == 10.0


def test_score_distance_zero_for_identical_rows():
    """Distance to reference should be zero for identical vectors."""
    row = pd.Series({"a": 4.0, "b": 8.0})
    assert score_distance(row, row, ["a", "b"]) == 0.0
