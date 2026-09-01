"""Helpers for labeling realstudy sampling regimes."""

"""Imports"""

import hashlib
from pathlib import Path

import pandas as pd

"""Functions"""


REALSTUDY_V2_ALGORITHM_VERSION = "sha256-rg-qname-v1"


def eligible_alignment(record) -> bool:
    """Return whether a BAM record is mapped, primary, and QC-passing."""
    return not (
        record.is_unmapped
        or record.is_secondary
        or record.is_supplementary
        or record.is_qcfail
    )


def fragment_identifier(record) -> str:
    """Return stable read-group plus query-name identity for single or paired reads."""
    read_group = record.get_tag("RG") if record.has_tag("RG") else "NO_RG"
    return f"{read_group}\x1f{record.query_name}"


def fragment_rank_digest(version: str, checksum: str, seed: int, identifier: str) -> str:
    """Return deterministic SHA-256 rank material for one fragment."""
    payload = "\x1e".join([version, checksum.lower(), str(seed), identifier]).encode()
    return hashlib.sha256(payload).hexdigest()


def coverage_label(value: float | int) -> str:
    """Return a filesystem-safe control-coverage label."""
    return f"{format_run_id_value(value)}x"


def target_fragments(genome_size_bp: int, fragment_length_bp: int, coverage: float) -> int:
    """Return the locked nearest-integer fragment target for one coverage."""
    return int(round((float(genome_size_bp) * float(coverage)) / float(fragment_length_bp)))


def build_realstudy_v2_run_rows(design: dict) -> list[dict]:
    """Expand the empirical two-study design into 42 samples and two anchors."""
    genome_size_bp = int(design["genome_size_bp"])
    fragment_length_bp = int(design["normalization_fragment_length_bp"])
    coverages = [float(value) for value in design["control_coverages_x"]]
    seeds = [int(value) for value in design["seeds"]]
    studies = design["studies"]
    rows = []
    for study_id, study in studies.items():
        anchor_id = f"{study_id}__full_control_anchor"
        common = {
            "study_id": study_id,
            "signal_class": study["signal_class"],
            "macs2_mode": study["macs2_mode"],
            "treatment_parent_id": f"{study_id}__treatment_parent",
            "control_parent_id": f"{study_id}__control_parent",
            "reference_run_id": anchor_id,
            "assembly": design["assembly"],
            "genome_size_bp": genome_size_bp,
            "normalization_fragment_length_bp": fragment_length_bp,
            "peakcaller": "macs2",
        }
        rows.append(
            {
                **common,
                "run_id": anchor_id,
                "matched_block_id": f"{study_id}__anchor_block",
                "run_type": "full_control_anchor",
                "seed": "",
                "control_coverage_x": "full",
                "requested_control_fragments": "",
                "control_subsample_id": "",
            }
        )
        for seed in seeds:
            block_id = f"{study_id}__seed_{seed}"
            for coverage in coverages:
                label = coverage_label(coverage)
                rows.append(
                    {
                        **common,
                        "run_id": f"{study_id}__control_{label}__seed_{seed}",
                        "matched_block_id": block_id,
                        "run_type": "control_subsample",
                        "seed": seed,
                        "control_coverage_x": coverage,
                        "requested_control_fragments": target_fragments(
                            genome_size_bp, fragment_length_bp, coverage
                        ),
                        "control_subsample_id": f"{study_id}__control_{label}__seed_{seed}",
                    }
                )
    return rows


def audit_realstudy_v2_design(rows: list[dict], design: dict) -> dict:
    """Validate the locked run expansion and return a compact design audit."""
    anchors = [row for row in rows if row["run_type"] == "full_control_anchor"]
    sampled = [row for row in rows if row["run_type"] == "control_subsample"]
    expected_studies = len(design["studies"])
    expected_depths = len(design["control_coverages_x"])
    expected_seeds = len(design["seeds"])
    expected_sampled = expected_studies * expected_depths * expected_seeds
    run_ids = [row["run_id"] for row in rows]
    if len(rows) != expected_sampled + expected_studies:
        raise ValueError("Realstudy v2 must expand to exactly 42 sampled calls and two anchors.")
    if len(sampled) != expected_sampled or len(anchors) != expected_studies:
        raise ValueError("Realstudy v2 run types do not match the locked design.")
    if len(set(run_ids)) != len(run_ids):
        raise ValueError("Realstudy v2 run identifiers are not unique.")
    maximum_target = max(int(row["requested_control_fragments"]) for row in sampled)
    if maximum_target != int(design["minimum_eligible_control_fragments"]):
        raise ValueError("Maximum target and parent-library preflight threshold disagree.")
    return {
        "studies": expected_studies,
        "control_depths": expected_depths,
        "seeds": expected_seeds,
        "sampled_peak_calls": len(sampled),
        "anchor_peak_calls": len(anchors),
        "total_peak_calls": len(rows),
        "maximum_target_fragments": maximum_target,
        "status": "valid",
    }


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
