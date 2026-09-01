"""Tests for preregistered Realstudy v2 sufficiency outcomes."""

"""Imports"""

from pathlib import Path
import subprocess
import sys

import pandas as pd
import yaml


"""Functions"""


def write_fixture(tmp_path: Path, pass_depths: set[float], parent_fragments: int = 30000000) -> list[str]:
    """Write a three-seed seven-depth anchor-relative metric fixture."""
    root = Path(__file__).resolve().parents[1]
    config = yaml.safe_load((root / "configs/realstudy_v2.yaml").read_text())
    config["realstudy_v2"]["studies"] = {"study_a": {"signal_class": "narrow", "macs2_mode": "narrow"}}
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config))
    runs = []
    metrics = []
    for depth in [0.5, 1, 2, 4, 8, 16, 32]:
        for seed in [11, 23, 37]:
            run_id = f"study_a__{depth}__{seed}"
            runs.append({"run_id": run_id, "study_id": "study_a", "run_type": "control_subsample", "control_coverage_x": depth, "seed": seed})
            passing = depth in pass_depths
            values = {
                "anchor_peak_retention": 0.95 if passing else 0.70,
                "query_peak_concordance": 0.95 if passing else 0.70,
                "genomic_base_jaccard": 0.90 if passing else 0.60,
                "peak_count_ratio": 1.0 if passing else 0.70,
            }
            for metric, value in values.items():
                metrics.append({"run_id": run_id, "metric": metric, "threshold": 0.0, "value": value, "status": "complete"})
    paths = {"runs": tmp_path / "runs.csv", "metrics": tmp_path / "metrics.csv", "parents": tmp_path / "parents.csv"}
    pd.DataFrame(runs).to_csv(paths["runs"], index=False)
    pd.DataFrame(metrics).to_csv(paths["metrics"], index=False)
    pd.DataFrame([{"study_id": "study_a", "role": "control", "eligible_fragments": parent_fragments}]).to_csv(paths["parents"], index=False)
    output = tmp_path / "output"
    command = [sys.executable, "scripts/realstudy_v2_analysis.py", "--config", str(config_path), "--run-table", str(paths["runs"]), "--run-metrics", str(paths["metrics"]), "--parent-libraries", str(paths["parents"]), "--output-dir", str(output)]
    subprocess.run(command, cwd=root, check=True)
    return [str(output), str(config_path)]


def test_sufficiency_selects_lowest_passing_plateau_depth(tmp_path: Path):
    """A passing depth with a passing next depth and small gain is selected."""
    output, _ = write_fixture(tmp_path, {16, 32})
    decision = pd.read_csv(Path(output) / "enough_control_decisions.csv").iloc[0]
    assert decision["outcome"] == "sufficient_depth_identified"
    assert decision["selected_coverage_x"] == 16


def test_highest_only_pass_is_agreement_without_plateau(tmp_path: Path):
    """Agreement first reached at 32× cannot establish the next-depth plateau."""
    output, _ = write_fixture(tmp_path, {32})
    decision = pd.read_csv(Path(output) / "enough_control_decisions.csv").iloc[0]
    assert decision["outcome"] == "agreement_without_plateau"


def test_insufficient_parent_is_not_silently_downscaled(tmp_path: Path):
    """A parent below the 32× target receives the explicit insufficient outcome."""
    output, _ = write_fixture(tmp_path, {16, 32}, parent_fragments=1000000)
    decision = pd.read_csv(Path(output) / "enough_control_decisions.csv").iloc[0]
    assert decision["outcome"] == "parent_library_insufficient"
