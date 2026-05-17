"""Helpers for labeling realstudy sampling regimes."""

"""Imports"""

from pathlib import Path

import pandas as pd

"""Functions"""


def format_run_id_value(value: float | int) -> str:
    """Format numeric run-id axes consistently across YAML and argparse inputs."""
    numeric = float(value)
    if numeric.is_integer():
        return str(int(numeric))
    return str(numeric).replace(".", "p")


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


def build_run_id(
    study_id: str,
    cov_treat: float,
    cov_ctrl: float,
    seed: int,
    fragment_length: int = 150,
    read_length: int = 38,
    aligner: str = "bowtie2",
    peakcaller: str = "macs2",
    macs2_mode: str = "narrow",
) -> str:
    """Build a deterministic run identifier."""
    normalized_study = study_id.replace("-", "_")
    return (
        f"{normalized_study}"
        f"_t{format_run_id_value(cov_treat)}"
        f"_c{format_run_id_value(cov_ctrl)}"
        f"_s{seed}"
        f"_fl{fragment_length}"
        f"_rl{read_length}"
        f"_al{aligner}"
        f"_pc{peakcaller}"
        f"_m{macs2_mode}"
    )


def load_observed_depths(path: Path | None) -> dict[str, dict]:
    """Load observed-depth metadata keyed by study identifier."""
    if path is None or not Path(path).exists():
        return {}
    observed_df = pd.read_csv(path)
    return {row["study_id"]: row for row in observed_df.to_dict(orient="records")}


def build_run_table_rows(
    *,
    study_ids: list[str],
    observed_by_study: dict[str, dict],
    coverage_treat: list[float],
    coverage_ctrl: list[float],
    seeds: list[int],
    fragment_lengths: list[int],
    read_lengths: list[int],
    aligners: list[str],
    peakcallers: list[str],
    macs2_modes: list[str],
    observed_treatment_depth: float | None = None,
    observed_control_depth: float = 1.0,
) -> list[dict]:
    """Build prototype run rows directly from source config and study depths."""
    rows = []
    for study_id in study_ids:
        observed = observed_by_study.get(
            study_id,
            {
                "observed_treatment_depth": observed_treatment_depth,
                "observed_control_depth": observed_control_depth,
                "assembly": "",
                "macs2_mode": "",
            },
        )
        mode_axis = [observed["macs2_mode"]] if str(observed.get("macs2_mode", "")).strip() else macs2_modes
        for cov_treat in coverage_treat:
            for cov_ctrl in coverage_ctrl:
                for seed in seeds:
                    for fragment_length in fragment_lengths:
                        for read_length in read_lengths:
                            for aligner in aligners:
                                for peakcaller in peakcallers:
                                    for macs2_mode in mode_axis:
                                        rows.append(
                                            {
                                                "run_id": build_run_id(
                                                    study_id,
                                                    cov_treat,
                                                    cov_ctrl,
                                                    seed,
                                                    fragment_length=fragment_length,
                                                    read_length=read_length,
                                                    aligner=aligner,
                                                    peakcaller=peakcaller,
                                                    macs2_mode=macs2_mode,
                                                ),
                                                "study_id": study_id,
                                                "coverage_treat": cov_treat,
                                                "coverage_ctrl": cov_ctrl,
                                                "seed": seed,
                                                "fragment_length": fragment_length,
                                                "read_length": read_length,
                                                "aligner": aligner,
                                                "peakcaller": peakcaller,
                                                "macs2_mode": macs2_mode,
                                                "assembly": observed.get("assembly", ""),
                                                "sampling_label": sampling_label(
                                                    cov_treat, float(observed["observed_treatment_depth"])
                                                ),
                                                "control_sampling_label": control_sampling_label(
                                                    cov_ctrl, float(observed.get("observed_control_depth", 1.0))
                                                ),
                                            }
                                        )
    return rows
