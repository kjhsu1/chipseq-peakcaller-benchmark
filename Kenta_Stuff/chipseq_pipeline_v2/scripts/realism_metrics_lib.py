"""Helpers for BAM-only realism scoring."""

"""Imports"""

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
import pysam


"""Data Structures"""


@dataclass(frozen=True)
class DepthStats:
    """Container for simple depth distribution summaries."""

    mean_depth: float
    median_depth: float
    cv_depth: float
    gini_depth: float
    bumpiness: float


"""Functions"""


def load_bed_intervals(path: Path) -> List[Tuple[str, int, int]]:
    """Load basic BED intervals."""
    intervals: List[Tuple[str, int, int]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, start, end, *_ = line.rstrip().split("\t")
            intervals.append((chrom, int(start), int(end)))
    return intervals


def interval_lengths(intervals: Sequence[Tuple[str, int, int]]) -> np.ndarray:
    """Return interval lengths as a numpy array."""
    return np.asarray([end - start for _, start, end in intervals], dtype=float)


def gini_coefficient(values: np.ndarray) -> float:
    """Return the Gini coefficient for nonnegative values."""
    if values.size == 0:
        return 0.0
    sorted_values = np.sort(np.asarray(values, dtype=float))
    if np.allclose(sorted_values, 0):
        return 0.0
    n = sorted_values.size
    index = np.arange(1, n + 1)
    return float(np.sum((2 * index - n - 1) * sorted_values) / (n * np.sum(sorted_values)))


def local_bumpiness(values: np.ndarray) -> float:
    """Return mean adjacent absolute change scaled by the mean depth."""
    if values.size < 2:
        return 0.0
    mean_value = float(np.mean(values))
    if mean_value == 0:
        return 0.0
    return float(np.mean(np.abs(np.diff(values))) / mean_value)


def metaprofile_metrics(profile: np.ndarray) -> dict:
    """Summarize a metaprofile with simple shape descriptors."""
    if profile.size == 0:
        return {
            "summit_height": 0.0,
            "half_max_width": 0.0,
            "auc": 0.0,
            "symmetry": 0.0,
            "shoulder_score": 0.0,
        }
    summit_height = float(np.max(profile))
    half_max = summit_height / 2.0
    above = np.where(profile >= half_max)[0]
    half_max_width = float(above[-1] - above[0] + 1) if above.size else 0.0
    mid = profile.size // 2
    left = profile[:mid]
    right = profile[-mid:][::-1]
    symmetry = float(1.0 - (np.mean(np.abs(left - right)) / summit_height)) if summit_height else 0.0
    shoulder_window = max(1, profile.size // 10)
    shoulder_values = np.concatenate([profile[:shoulder_window], profile[-shoulder_window:]])
    shoulder_score = float(np.mean(shoulder_values) / summit_height) if summit_height else 0.0
    return {
        "summit_height": summit_height,
        "half_max_width": half_max_width,
        "auc": float(np.sum(profile)),
        "symmetry": symmetry,
        "shoulder_score": shoulder_score,
    }


def estimate_depths(bam_path: Path, intervals: Sequence[Tuple[str, int, int]]) -> np.ndarray:
    """Estimate per-base depths inside the provided intervals."""
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    depth_chunks: List[np.ndarray] = []
    try:
        for chrom, start, end in intervals:
            a_arr, c_arr, g_arr, t_arr = bam.count_coverage(chrom, start, end, quality_threshold=0)
            depth_chunks.append(
                np.asarray(a_arr, dtype=float)
                + np.asarray(c_arr, dtype=float)
                + np.asarray(g_arr, dtype=float)
                + np.asarray(t_arr, dtype=float)
            )
    finally:
        bam.close()
    if not depth_chunks:
        return np.asarray([], dtype=float)
    return np.concatenate(depth_chunks)


def summarize_depths(depths: np.ndarray) -> DepthStats:
    """Summarize a depth vector."""
    if depths.size == 0:
        return DepthStats(0.0, 0.0, 0.0, 0.0, 0.0)
    mean_depth = float(np.mean(depths))
    std_depth = float(np.std(depths))
    return DepthStats(
        mean_depth=mean_depth,
        median_depth=float(np.median(depths)),
        cv_depth=(std_depth / mean_depth) if mean_depth else 0.0,
        gini_depth=gini_coefficient(depths),
        bumpiness=local_bumpiness(depths),
    )


def build_centered_windows(
    intervals: Sequence[Tuple[str, int, int]], half_window: int
) -> List[Tuple[str, int, int]]:
    """Build summit-centered windows from peak intervals."""
    windows: List[Tuple[str, int, int]] = []
    for chrom, start, end in intervals:
        center = (start + end) // 2
        windows.append((chrom, max(0, center - half_window), center + half_window))
    return windows


def mean_metaprofile(bam_path: Path, intervals: Sequence[Tuple[str, int, int]], half_window: int) -> np.ndarray:
    """Return an average treatment-depth profile around interval centers."""
    windows = build_centered_windows(intervals, half_window)
    profiles: List[np.ndarray] = []
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        for chrom, start, end in windows:
            a_arr, c_arr, g_arr, t_arr = bam.count_coverage(chrom, start, end, quality_threshold=0)
            profile = (
                np.asarray(a_arr, dtype=float)
                + np.asarray(c_arr, dtype=float)
                + np.asarray(g_arr, dtype=float)
                + np.asarray(t_arr, dtype=float)
            )
            if profile.size == (2 * half_window):
                profiles.append(profile)
    finally:
        bam.close()
    if not profiles:
        return np.asarray([], dtype=float)
    return np.mean(np.vstack(profiles), axis=0)


def score_distance(sim_metrics: pd.Series, ref_metrics: pd.Series, metric_names: Iterable[str]) -> float:
    """Return a simple normalized distance to a reference metric vector."""
    parts = []
    for metric_name in metric_names:
        sim_value = float(sim_metrics[metric_name])
        ref_value = float(ref_metrics[metric_name])
        scale = max(abs(ref_value), 1.0)
        parts.append(abs(sim_value - ref_value) / scale)
    return float(np.mean(parts)) if parts else 0.0
