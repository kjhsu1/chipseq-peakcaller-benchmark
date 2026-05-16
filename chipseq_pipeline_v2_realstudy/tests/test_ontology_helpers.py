"""Tests for ontology classification helpers."""

"""Imports"""

import pandas as pd

from scripts.ontology_lib import classify_row


"""Functions"""


def test_classify_row_uses_signal_and_bumpiness_thresholds():
    """Ontology classes should distinguish hotspot, peak, and background regions."""
    hotspot = classify_row(
        pd.Series({"ctrl_cv": 0.5, "ctrl_bump_z": 4.0, "treat_peak_z": 3.5, "log2_enrichment": 2.0})
    )
    peak = classify_row(
        pd.Series({"ctrl_cv": 0.1, "ctrl_bump_z": 0.2, "treat_peak_z": 7.0, "log2_enrichment": 2.0})
    )
    background = classify_row(
        pd.Series({"ctrl_cv": 0.1, "ctrl_bump_z": 0.1, "treat_peak_z": 1.0, "log2_enrichment": 0.2})
    )
    assert hotspot["background_class"] == "hilly_background"
    assert peak["foreground_class"] == "strong_peak"
    assert background["foreground_class"] == "background_only"
