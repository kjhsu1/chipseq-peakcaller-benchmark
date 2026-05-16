"""Helpers for labeling realstudy sampling regimes."""

"""Functions"""


def sampling_label(target_depth: float, observed_depth: float) -> str:
    """Label whether requested sampling is below, at, or above observed depth."""
    if target_depth < observed_depth:
        return "sampling_below_observed_depth"
    if target_depth == observed_depth:
        return "sampling_at_observed_depth"
    return "model_resampled_above_observed_depth"


def control_sampling_label(target_depth: float, observed_depth: float) -> str:
    """Label the requested control-depth regime."""
    return sampling_label(target_depth, observed_depth)


def build_run_id(study_id: str, cov_treat: float, cov_ctrl: float, seed: int) -> str:
    """Build a deterministic run identifier."""
    normalized_study = study_id.replace("-", "_")
    return f"{normalized_study}_t{str(cov_treat).replace('.', 'p')}_c{str(cov_ctrl).replace('.', 'p')}_s{seed}"
