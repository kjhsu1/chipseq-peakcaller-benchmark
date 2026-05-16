"""Helpers for controlled-category cleanup summaries."""

"""Imports"""

from pathlib import Path
from typing import Dict, Iterable, List

import pandas as pd


"""Constants"""

IGNORE_SWEEP_COLUMNS = {
    "run_id",
    "id_ctrl",
    "id_treat",
    "total_called",
    "total_planted",
    "tp_called",
    "tp_planted",
    "precision",
    "recall",
    "f1",
}


"""Functions"""


def parse_simple_yaml_map(path: Path) -> Dict[str, str]:
    """Parse a simple key-value YAML mapping file."""
    mapping: Dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        key, value = line.split(":", 1)
        mapping[key.strip()] = value.strip()
    return mapping


def canonical_label(config_name: str, mapping: Dict[str, str]) -> str:
    """Return the configured canonical label or the original name."""
    return mapping.get(config_name, config_name)


def expected_run_count(per_run_df: pd.DataFrame) -> int:
    """Estimate expected run count from observed sweep columns."""
    expected = 1
    for column in per_run_df.columns:
        if column in IGNORE_SWEEP_COLUMNS:
            continue
        unique_count = per_run_df[column].nunique(dropna=True)
        if unique_count > 1:
            expected *= unique_count
    return expected


def zero_metric_row_count(per_run_df: pd.DataFrame) -> int:
    """Count rows with any core metric at zero."""
    metric_cols = ["total_called", "tp_called", "tp_planted", "precision", "recall", "f1"]
    return int((per_run_df[metric_cols] == 0).any(axis=1).sum())


def best_control_by_mean_f1(per_run_df: pd.DataFrame) -> float:
    """Return the control coverage with the best mean F1."""
    grouped = per_run_df.groupby("coverage_ctrl", as_index=False)["f1"].mean()
    grouped = grouped.sort_values(["f1", "coverage_ctrl"], ascending=[False, True])
    return float(grouped.iloc[0]["coverage_ctrl"])


def summarize_counts_by_control(per_run_df: pd.DataFrame) -> pd.DataFrame:
    """Summarize mean called, TP, and FP counts by control depth."""
    summary = (
        per_run_df.assign(fp_called=per_run_df["total_called"] - per_run_df["tp_called"])
        .groupby("coverage_ctrl", as_index=False)
        .agg(
            mean_called=("total_called", "mean"),
            mean_tp_called=("tp_called", "mean"),
            mean_fp_called=("fp_called", "mean"),
            mean_precision=("precision", "mean"),
            mean_recall=("recall", "mean"),
            mean_f1=("f1", "mean"),
        )
        .sort_values("coverage_ctrl")
    )
    return summary


def recall_saturation_points(plot_df: pd.DataFrame, fraction_of_max: float = 0.99) -> pd.DataFrame:
    """Return first control depth where recall reaches a near-maximum plateau."""
    rows: List[dict] = []
    for coverage_treat, group in plot_df.groupby("coverage_treat"):
        ordered = group.sort_values("coverage_ctrl")
        max_recall = float(ordered["recall"].max())
        threshold = max_recall * fraction_of_max
        hit = ordered[ordered["recall"] >= threshold].iloc[0]
        rows.append(
            {
                "coverage_treat": float(coverage_treat),
                "recall_saturation_ctrl": float(hit["coverage_ctrl"]),
                "max_recall": max_recall,
            }
        )
    return pd.DataFrame(rows).sort_values("coverage_treat")


def precision_dropoff_points(plot_df: pd.DataFrame, retained_fraction: float = 0.8) -> pd.DataFrame:
    """Return first control depth where precision drops materially from the best point."""
    rows: List[dict] = []
    for coverage_treat, group in plot_df.groupby("coverage_treat"):
        ordered = group.sort_values("coverage_ctrl")
        max_precision = float(ordered["precision"].max())
        threshold = max_precision * retained_fraction
        below = ordered[ordered["precision"] <= threshold]
        hit = below.iloc[0] if not below.empty else ordered.iloc[-1]
        rows.append(
            {
                "coverage_treat": float(coverage_treat),
                "precision_dropoff_ctrl": float(hit["coverage_ctrl"]),
                "max_precision": max_precision,
            }
        )
    return pd.DataFrame(rows).sort_values("coverage_treat")


def observed_parameter_values(per_run_df: pd.DataFrame) -> Dict[str, List[str]]:
    """Collect observed parameter values for varying columns."""
    observed: Dict[str, List[str]] = {}
    for column in per_run_df.columns:
        if column in IGNORE_SWEEP_COLUMNS:
            continue
        values = sorted(per_run_df[column].dropna().astype(str).unique().tolist())
        if len(values) > 1:
            observed[column] = values
    return observed


def warning_lines() -> Iterable[str]:
    """Return standard scientific interpretation warnings."""
    return [
        "plateau vs peak naming is only partly trustworthy",
        "MACS2 broad/narrow mode is confounded with TF-shape settings in some categories",
        "current hilly labeling corresponds to a dense hotspot process, not '3% of peaks'",
    ]
