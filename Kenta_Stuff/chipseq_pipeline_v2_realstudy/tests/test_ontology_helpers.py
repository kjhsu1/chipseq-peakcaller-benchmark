"""Tests for ontology classification helpers."""

"""Imports"""

import pandas as pd

from scripts.ontology_lib import classify_row


"""Functions"""


def test_classify_row_uses_signal_and_bumpiness_thresholds():
    """Ontology classes should distinguish hotspot, peak, and background regions."""
    assert classify_row(pd.Series({"signal_mean": 5.0, "bumpiness": 0.9})) == "hilly_hotspot"
    assert classify_row(pd.Series({"signal_mean": 5.0, "bumpiness": 0.2})) == "strong_peak"
    assert classify_row(pd.Series({"signal_mean": 1.0, "bumpiness": 0.1})) == "background"
