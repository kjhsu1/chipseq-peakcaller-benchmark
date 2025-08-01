"""Simplified CHIP-seq helpers used for tests."""

import os
from typing import List, Dict

import numpy as np
import pandas as pd
from scipy.stats import norm

k = 1


def build_tf_bias_pmf(length: int, peaks: List[int], sigma: float, enrichment: float) -> np.ndarray:
    """Return normalized TF-binding bias PMF."""
    bias = np.ones(length, dtype=float)
    positions = np.arange(length)
    for p in peaks:
        kernel = norm.pdf(positions, loc=p, scale=sigma)
        bias += enrichment * kernel
    return bias / bias.sum()


def build_gc_bias_pmf(sequence: str, loess_params: Dict, fragment_length: int = k) -> np.ndarray:
    """Return normalized GC-content bias PMF."""
    if not loess_params or 'csv' not in loess_params or loess_params['csv'] is None:
        bias = np.ones(len(sequence) - fragment_length + 1, dtype=float)
        return bias / bias.sum()
    csv_path = loess_params['csv']
    if not os.path.exists(csv_path):
        raise FileNotFoundError(csv_path)
    table = pd.read_csv(csv_path)
    gc_vals = table.iloc[:, 0].to_numpy()
    weights = table.iloc[:, 1].to_numpy()
    seq_arr = np.frombuffer(sequence.upper().encode('ascii'), dtype='S1')
    gc_mask = np.isin(seq_arr, [b'G', b'C'])
    cumsum = np.cumsum(gc_mask)
    counts = cumsum[fragment_length - 1:] - np.concatenate(([0], cumsum[:-fragment_length]))
    gc_percent = counts / fragment_length
    bias = np.interp(gc_percent, gc_vals, weights)
    return bias / bias.sum()


def build_accessibility_bias_pmf(length: int, accessibility_bed: str, acc_weight: float) -> np.ndarray:
    """Return normalized accessibility bias PMF."""
    bias = np.ones(length, dtype=float)
    if accessibility_bed:
        if not os.path.exists(accessibility_bed):
            raise FileNotFoundError(accessibility_bed)
        with open(accessibility_bed) as fh:
            for line in fh:
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.strip().split()[:3]
                if len(fields) < 3:
                    continue
                start = max(int(fields[1]), 0)
                end = min(int(fields[2]), length)
                bias[start:end] *= acc_weight
    return bias / bias.sum()
