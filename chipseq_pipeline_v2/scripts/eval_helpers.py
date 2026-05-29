"""Shared helpers for peak-evaluation summaries and path resolution."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd


PEAKLIKE_MAX_SIGMA = 5
PLATEAULIKE_MIN_SIGMA = 15


def resolve_peak_path(results_dir: Path, run_id: str, peakcaller: str) -> Path:
    """Return the expected called-peak path for one run."""
    if peakcaller == "epic2":
        return results_dir / run_id / "peaks" / "epic2" / f"{run_id}_domains.bed"
    if peakcaller == "homer":
        return results_dir / run_id / "peaks" / "homer" / f"{run_id}_peaks.bed"

    macs2_dir = results_dir / run_id / "peaks" / "macs2"
    candidate_paths = [
        macs2_dir / f"{run_id}_peaks.bed",
        macs2_dir / f"{run_id}_peaks.narrowPeak",
        macs2_dir / f"{run_id}_peaks.broadPeak",
    ]
    return next((path for path in candidate_paths if path.exists()), candidate_paths[0])


def expected_decode_mode(tf_sigma: float) -> str:
    """Return the expected broad/narrow decoding mode for one planted shape."""
    if tf_sigma <= PEAKLIKE_MAX_SIGMA:
        return "narrow"
    if tf_sigma >= PLATEAULIKE_MIN_SIGMA:
        return "broad"
    raise ValueError(
        "Cannot infer expected decode mode for tf_sigma values between "
        f"{PEAKLIKE_MAX_SIGMA} and {PLATEAULIKE_MIN_SIGMA}."
    )


def validate_decode_modes(rows: Iterable[Dict[str, object]]) -> None:
    """Raise if any run row mixes planted shape with the wrong decode mode."""
    mismatches = []
    for row in rows:
        tf_sigma = float(row["tf_sigma"])
        observed = str(row.get("macs2_mode", "narrow"))
        expected = expected_decode_mode(tf_sigma)
        if observed != expected:
            mismatches.append(
                {
                    "run_id": row.get("run_id", "unknown"),
                    "tf_sigma": tf_sigma,
                    "observed_mode": observed,
                    "expected_mode": expected,
                }
            )
    if not mismatches:
        return
    first = mismatches[0]
    raise ValueError(
        "Decode-mode mismatch detected between planted shape and caller mode. "
        f"First mismatch: run_id={first['run_id']}, tf_sigma={first['tf_sigma']}, "
        f"observed={first['observed_mode']}, expected={first['expected_mode']}. "
        "Peak-like shapes (tf_sigma <= 5) must use narrow mode; "
        "plateau-like shapes (tf_sigma >= 15) must use broad mode."
    )


def truth_mode_for_row(row: Dict[str, object]) -> str:
    """Return the planted-truth representation used for one run row."""
    peakcaller = str(row.get("peakcaller", "macs2"))
    if peakcaller == "epic2":
        return "interval"
    return expected_decode_mode(float(row["tf_sigma"]))


def counts_to_metrics(
    tp_called: float,
    total_called: float,
    tp_planted: float,
    total_planted: float,
) -> Tuple[float, float, float]:
    """Convert overlap counts into precision, recall, and F1."""
    precision = tp_called / total_called if total_called else 0.0
    recall = tp_planted / total_planted if total_planted else 0.0
    if precision + recall == 0:
        return precision, recall, 0.0
    f1 = 2 * (precision * recall) / (precision + recall)
    return precision, recall, f1


def aggregate_counts_summary(
    data: pd.DataFrame,
    group_cols: List[str],
) -> pd.DataFrame:
    """Aggregate overlap counts using ratio-of-summed-counts."""
    summary = (
        data.groupby(group_cols, as_index=False)
        .agg(
            tp_called=("tp_called", "sum"),
            total_called=("total_called", "sum"),
            tp_planted=("tp_planted", "sum"),
            total_planted=("total_planted", "sum"),
            n_runs=("run_id", "count"),
        )
        .sort_values(group_cols)
    )
    metrics = summary.apply(
        lambda row: counts_to_metrics(
            row["tp_called"],
            row["total_called"],
            row["tp_planted"],
            row["total_planted"],
        ),
        axis=1,
        result_type="expand",
    )
    metrics.columns = ["precision", "recall", "f1"]
    return pd.concat([summary, metrics], axis=1)


def aggregate_mean_of_runs_summary(
    data: pd.DataFrame,
    group_cols: List[str],
) -> pd.DataFrame:
    """Aggregate mean per-run metrics for transparency only."""
    summary = (
        data.groupby(group_cols)[["precision", "recall", "f1"]]
        .mean()
        .reset_index()
        .sort_values(group_cols)
    )
    summary["n_runs"] = data.groupby(group_cols)["run_id"].count().values
    return summary


def metric_definition_lines(
    *,
    truth_mode: str,
    aggregation_rule: str,
    caller_settings: str,
    evaluation_scope: str,
) -> List[str]:
    """Return standard metadata lines describing metric semantics."""
    return [
        "# Metric Definition",
        "",
        f"- truth mode: `{truth_mode}`",
        f"- aggregation rule: `{aggregation_rule}`",
        f"- caller settings: `{caller_settings}`",
        f"- evaluation scope: `{evaluation_scope}`",
    ]
