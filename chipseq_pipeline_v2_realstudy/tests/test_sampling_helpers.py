"""Tests for prototype sampling-table helpers."""

"""Imports"""

from pathlib import Path

import pandas as pd

from scripts.realstudy_sampling_lib import (
    build_run_id,
    build_run_table_rows,
    control_sampling_label,
    load_observed_depths,
    sampling_label,
)


def test_build_run_id_is_deterministic():
    """Run identifiers should preserve study and sweep coordinates."""
    assert (
        build_run_id("study-a", 10, 4, 11)
        == "study_a_t10_c4_s11_fl150_rl38_albowtie2_pcmacs2_mnarrow"
    )


def test_build_run_id_normalizes_integer_valued_floats():
    """Run identifiers should match whether depths come from YAML or argparse."""
    assert (
        build_run_id("study-a", 5.0, 0.5, 11)
        == "study_a_t5_c0p5_s11_fl150_rl38_albowtie2_pcmacs2_mnarrow"
    )


def test_control_sampling_label_matches_depth_regimes():
    """Control labels should mirror the same below/at/above semantics."""
    assert control_sampling_label(2, 4) == "sampling_below_observed_depth"
    assert control_sampling_label(4, 4) == "sampling_at_observed_depth"
    assert control_sampling_label(8, 4) == "model_resampled_above_observed_depth"


def test_load_observed_depths_reads_csv(tmp_path: Path):
    """Observed-depth helper should preserve study metadata by study id."""
    path = tmp_path / "depths.csv"
    pd.DataFrame(
        [
            {
                "study_id": "study-a",
                "observed_treatment_depth": 10,
                "observed_control_depth": 4,
                "assembly": "ce11",
                "macs2_mode": "broad",
            }
        ]
    ).to_csv(path, index=False)
    observed = load_observed_depths(path)
    assert observed["study-a"]["assembly"] == "ce11"


def test_build_run_table_rows_uses_all_axes_in_run_id():
    """Run-table rows should stay unique when non-default axes are present."""
    rows = build_run_table_rows(
        study_ids=["study-a"],
        observed_by_study={
            "study-a": {
                "study_id": "study-a",
                "observed_treatment_depth": 10,
                "observed_control_depth": 4,
                "assembly": "ce11",
                "macs2_mode": "broad",
            }
        },
        coverage_treat=[10],
        coverage_ctrl=[4],
        seeds=[11],
        fragment_lengths=[150, 200],
        read_lengths=[38],
        aligners=["bowtie2"],
        peakcallers=["macs2"],
        macs2_modes=["narrow", "broad"],
    )
    assert len(rows) == 2
    assert len({row["run_id"] for row in rows}) == 2
    assert {row["macs2_mode"] for row in rows} == {"broad"}
