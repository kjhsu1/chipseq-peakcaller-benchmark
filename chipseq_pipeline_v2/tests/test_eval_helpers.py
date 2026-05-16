"""Tests for shared evaluation helpers."""

from pathlib import Path

import pandas as pd

from scripts.eval_helpers import (
    aggregate_counts_summary,
    aggregate_mean_of_runs_summary,
    resolve_peak_path,
)


def test_aggregate_counts_summary_matches_ratio_of_summed_counts():
    """Counts-based aggregation should recompute metrics from summed counts."""
    data = pd.DataFrame(
        [
            {"run_id": "0001", "coverage_ctrl": 1.0, "tp_called": 4, "total_called": 8, "tp_planted": 6, "total_planted": 10, "precision": 0.5, "recall": 0.6, "f1": 0.5455},
            {"run_id": "0002", "coverage_ctrl": 1.0, "tp_called": 1, "total_called": 2, "tp_planted": 2, "total_planted": 10, "precision": 0.5, "recall": 0.2, "f1": 0.2857},
        ]
    )
    summary = aggregate_counts_summary(data, ["coverage_ctrl"])
    row = summary.iloc[0]
    assert row["tp_called"] == 5
    assert row["total_called"] == 10
    assert row["tp_planted"] == 8
    assert row["total_planted"] == 20
    assert row["precision"] == 0.5
    assert row["recall"] == 0.4


def test_aggregate_mean_of_runs_summary_keeps_secondary_view():
    """Mean-of-runs aggregation should stay available as a separate view."""
    data = pd.DataFrame(
        [
            {"run_id": "0001", "coverage_ctrl": 1.0, "precision": 1.0, "recall": 0.5, "f1": 2 / 3},
            {"run_id": "0002", "coverage_ctrl": 1.0, "precision": 0.0, "recall": 1.0, "f1": 0.0},
        ]
    )
    summary = aggregate_mean_of_runs_summary(data, ["coverage_ctrl"])
    row = summary.iloc[0]
    assert row["precision"] == 0.5
    assert row["recall"] == 0.75
    assert row["n_runs"] == 2


def test_resolve_peak_path_supports_epic2_and_macs2(tmp_path):
    """Peak resolution should follow caller-specific layouts."""
    results_dir = tmp_path / "results"
    macs2_peak = results_dir / "0001" / "peaks" / "macs2" / "0001_peaks.bed"
    macs2_peak.parent.mkdir(parents=True)
    macs2_peak.write_text("chr1\t1\t2\n", encoding="utf-8")
    epic2_peak = results_dir / "0002" / "peaks" / "epic2" / "0002_domains.bed"
    epic2_peak.parent.mkdir(parents=True)
    epic2_peak.write_text("chr1\t1\t5\n", encoding="utf-8")

    assert resolve_peak_path(results_dir, "0001", "macs2") == macs2_peak
    assert resolve_peak_path(results_dir, "0002", "epic2") == epic2_peak
