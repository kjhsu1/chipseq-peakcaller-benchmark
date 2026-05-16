"""Tests for prototype sampling-table helpers."""

"""Imports"""

import pandas as pd

from scripts.realstudy_sampling_lib import build_run_id, control_sampling_label, sampling_label


def test_build_run_id_is_deterministic():
    """Run identifiers should preserve study and sweep coordinates."""
    assert build_run_id("study-a", 10, 4, 11) == "study_a_t10_c4_s11"


def test_control_sampling_label_matches_depth_regimes():
    """Control labels should mirror the same below/at/above semantics."""
    assert control_sampling_label(2, 4) == "sampling_below_observed_depth"
    assert control_sampling_label(4, 4) == "sampling_at_observed_depth"
    assert control_sampling_label(8, 4) == "model_resampled_above_observed_depth"
