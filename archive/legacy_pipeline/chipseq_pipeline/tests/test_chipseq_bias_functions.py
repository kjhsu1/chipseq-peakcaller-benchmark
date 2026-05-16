"""Simplified CHIP-seq helpers used for tests."""

"""Imports"""

import os
from typing import List, Dict

import numpy as np
import pandas as pd
from scipy.stats import norm

"""Constants"""

k = 1


def build_tf_bias_pmf(length: int, peaks: List[int], sigma: float,
                      enrichment: float, exp: float = 1.0) -> np.ndarray:
    """Return TF-binding PMF reshaped by exponent."""
    bias = np.ones(length, dtype=float)
    positions = np.arange(length)
    for p in peaks:
        kernel = norm.pdf(positions, loc=p, scale=sigma)
        bias += enrichment * kernel
    bias /= bias.sum()
    if exp != 1.0:
        bias = np.power(bias, exp)
        bias /= bias.sum()
    return bias


def build_gc_bias_pmf(sequence: str, loess_params: Dict,
                      fragment_length: int = k, exp: float = 1.0) -> np.ndarray:
    """Return GC-content PMF reshaped by exponent."""
    if not loess_params or 'csv' not in loess_params or loess_params['csv'] is None:
        bias = np.ones(len(sequence) - fragment_length + 1, dtype=float)
        bias /= bias.sum()
        if exp != 1.0:
            bias = np.power(bias, exp)
            bias /= bias.sum()
        return bias
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
    bias /= bias.sum()
    if exp != 1.0:
        bias = np.power(bias, exp)
        bias /= bias.sum()
    return bias


def build_accessibility_bias_pmf(length: int, accessibility_bed: str,
                                 acc_weight: float, chrom_id: str,
                                 exp: float = 1.0) -> np.ndarray:
    """Return accessibility PMF reshaped by exponent."""
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
                if fields[0] != chrom_id:
                    continue
                start = max(int(fields[1]), 0)
                end = min(int(fields[2]), length)
                bias[start:end] *= acc_weight
    bias /= bias.sum()
    if exp != 1.0:
        bias = np.power(bias, exp)
        bias /= bias.sum()
    return bias


def test_tf_bias_exponent_effect():
    """Ensure TF exponent reshapes distribution."""
    flat = build_tf_bias_pmf(10, [5], 1.0, 1.0, exp=0.5)
    sharp = build_tf_bias_pmf(10, [5], 1.0, 1.0, exp=2.0)
    assert np.isclose(flat.sum(), 1.0)
    assert np.isclose(sharp.sum(), 1.0)
    assert sharp.max() > flat.max()


def test_gc_bias_exponent_effect(tmp_path):
    """Ensure GC exponent reshapes distribution."""
    table = pd.DataFrame({'gc': [0.0, 1.0], 'w': [1, 2]})
    csv = tmp_path / 'gc.csv'
    table.to_csv(csv, index=False)
    seq = 'GGGGGAAAAA'
    flat = build_gc_bias_pmf(seq, {'csv': str(csv)}, 5, exp=0.5)
    sharp = build_gc_bias_pmf(seq, {'csv': str(csv)}, 5, exp=2.0)
    assert np.isclose(flat.sum(), 1.0)
    assert np.isclose(sharp.sum(), 1.0)
    assert sharp.max() > flat.max()


def test_accessibility_bias_exponent_effect(tmp_path):
    """Ensure accessibility exponent reshapes distribution."""
    bed = tmp_path / 'acc.bed'
    bed.write_text('chr1\t2\t5\n')
    flat = build_accessibility_bias_pmf(10, str(bed), 2.0, 'chr1', exp=0.5)
    sharp = build_accessibility_bias_pmf(10, str(bed), 2.0, 'chr1', exp=2.0)
    assert np.isclose(flat.sum(), 1.0)
    assert np.isclose(sharp.sum(), 1.0)
    assert sharp.max() > flat.max()
