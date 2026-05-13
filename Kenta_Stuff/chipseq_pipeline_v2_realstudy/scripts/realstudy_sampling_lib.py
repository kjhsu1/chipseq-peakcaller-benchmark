"""Helpers for labeling realstudy sampling regimes."""

"""Functions"""


def sampling_label(target_depth: float, observed_depth: float) -> str:
    """Label whether requested sampling is below, at, or above observed depth."""
    if target_depth < observed_depth:
        return "sampling_below_observed_depth"
    if target_depth == observed_depth:
        return "sampling_at_observed_depth"
    return "model_resampled_above_observed_depth"
