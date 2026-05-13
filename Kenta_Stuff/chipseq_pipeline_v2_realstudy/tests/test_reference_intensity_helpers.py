"""Tests for reference-intensity sampling helpers."""

"""Imports"""

from scripts.realstudy_sampling_lib import sampling_label


"""Functions"""


def test_sampling_label_covers_below_at_and_above():
    """Sampling labels should distinguish depth regimes."""
    assert sampling_label(5, 10) == "sampling_below_observed_depth"
    assert sampling_label(10, 10) == "sampling_at_observed_depth"
    assert sampling_label(20, 10) == "model_resampled_above_observed_depth"
