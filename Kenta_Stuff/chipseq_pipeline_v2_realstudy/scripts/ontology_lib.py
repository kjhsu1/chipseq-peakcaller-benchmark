"""Helpers for simple ontology classification."""

"""Imports"""

import pandas as pd


"""Functions"""


def classify_row(row: pd.Series) -> str:
    """Assign a coarse ontology class."""
    if row["signal_mean"] >= 4 and row["bumpiness"] >= 0.8:
        return "hilly_hotspot"
    if row["signal_mean"] >= 4:
        return "strong_peak"
    if row["signal_mean"] >= 2:
        return "weak_peak"
    return "background"
