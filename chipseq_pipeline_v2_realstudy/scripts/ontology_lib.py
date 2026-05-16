"""Helpers for ontology classification and summaries."""

"""Imports"""

import pandas as pd


"""Functions"""


DEFAULT_THRESHOLDS = {
    "ctrl_flat_cv_max": 0.25,
    "ctrl_wavy_cv_min": 0.25,
    "ctrl_hilly_bump_z_min": 3.0,
    "treat_peak_z_min": 3.0,
    "treat_strong_peak_z_min": 6.0,
    "treat_plateau_width_min_bp": 1000.0,
    "log2_enrichment_min": 1.0,
}


def _value(row: pd.Series, name: str, default: float = 0.0) -> float:
    """Return a numeric row value with a safe default."""
    if name not in row:
        return default
    value = row[name]
    if pd.isna(value):
        return default
    return float(value)


def classify_row(row: pd.Series, thresholds: dict | None = None) -> dict:
    """Assign ontology classes, confounders, and combined labels."""
    limits = dict(DEFAULT_THRESHOLDS)
    if thresholds:
        limits.update(thresholds)

    ctrl_cv = _value(row, "ctrl_cv")
    ctrl_bump_z = _value(row, "ctrl_bump_z")
    treat_peak_z = _value(row, "treat_peak_z", _value(row, "signal_mean"))
    log2_enrichment = _value(row, "log2_enrichment")
    plateau_width_bp = _value(row, "plateau_width_bp")

    if ctrl_bump_z >= limits["ctrl_hilly_bump_z_min"]:
        background_class = "hilly_background"
    elif ctrl_cv >= limits["ctrl_wavy_cv_min"]:
        background_class = "wavy_background"
    elif ctrl_cv <= limits["ctrl_flat_cv_max"]:
        background_class = "flat_background"
    else:
        background_class = "mixed_background"

    if plateau_width_bp >= limits["treat_plateau_width_min_bp"] and log2_enrichment >= limits["log2_enrichment_min"]:
        foreground_class = "plateau_like_signal"
    elif treat_peak_z >= limits["treat_strong_peak_z_min"]:
        foreground_class = "strong_peak"
    elif treat_peak_z >= limits["treat_peak_z_min"]:
        foreground_class = "weak_peak"
    else:
        foreground_class = "background_only"

    confounders = []
    if background_class == "hilly_background" and foreground_class in {"weak_peak", "background_only"}:
        confounders.append("hotspot_like_background")
    if plateau_width_bp >= limits["treat_plateau_width_min_bp"] and foreground_class != "plateau_like_signal":
        confounders.append("plateau_width_conflict")
    if log2_enrichment < limits["log2_enrichment_min"] and foreground_class != "background_only":
        confounders.append("weak_enrichment")

    failure_mode = "none"
    if foreground_class == "background_only" and background_class != "flat_background":
        failure_mode = "background_dominant"
    elif confounders:
        failure_mode = "confounded_signal"

    ontology_class = f"{background_class}__{foreground_class}"
    return {
        "background_class": background_class,
        "foreground_class": foreground_class,
        "confounders": ";".join(confounders) if confounders else "none",
        "failure_mode": failure_mode,
        "ontology_class": ontology_class,
    }
