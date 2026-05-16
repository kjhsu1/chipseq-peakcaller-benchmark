"""
This program models a Chip-Seq experiment as a sampler drawing from a combined probability mass function
of the background and foreground

This program will print FASTA of reads generated based on user defined arguments (coverage, number of peaks in fg, bg, etc)
"""

"""Imports"""

import random
import os
from typing import List, Dict, Tuple

import numpy as np
import pandas as pd
from scipy.stats import norm

import scripts.lib as lib
import argparse


"""
Assumptions for our Model
_________________________

- We will just make the assumption that all fragments are uniform in this model
- CHIP-seq reads are anywhere from 100-500bp approx.
"""


"""
Command line Arguments
________________
"""

parser = argparse.ArgumentParser(
    description="CHIP-seq Experiment Simulator with bias modeling"
)
parser.add_argument('--fasta', type=str, required=False, default='',
    help='Path to the FASTA file')
parser.add_argument('--coverage', type=float, default=1.0,
    help='Coverage depth to simulate')
parser.add_argument('--tf_peak_count', type=int, default=1,
    help='Number of TF peaks to simulate')
parser.add_argument('--fragment_length', type=int, default=150,
    help='Length (bp) of the dsDNA fragment; assumed constant.')
parser.add_argument('--read_length', type=int, default=38,
    help='Length (bp) of each mate-pair read (typical kits: 38 bp or 100 bp).')
parser.add_argument('--tf_sigma', type=float, default=5.0,
    help='Standard deviation for TF-binding bias')
parser.add_argument('--tf_enrich', '--tf_enrichment', dest='tf_enrichment', type=float, default=1.0,
    help='Enrichment factor for TF-binding bias')
parser.add_argument('--accessibility_bed', type=str, default=None,
    help='BED file describing open chromatin intervals')
parser.add_argument('--acc_weight', type=float, default=1.0,
    help='Weight multiplier for accessible regions')
parser.add_argument('--gc_bias_params', type=str, default=None,
    help='CSV file for GC bias lookup table')
parser.add_argument('--tf_exp', type=float, default=1.0,
    help='Exponent to reshape TF PMF (<1 flattens, >1 sharpens).')
parser.add_argument('--gc_exp', type=float, default=1.0,
    help='Exponent to reshape GC PMF (<1 flattens, >1 sharpens).')
parser.add_argument('--acc_exp', type=float, default=1.0,
    help='Exponent to reshape accessibility PMF (<1 flattens, >1 sharpens).')
parser.add_argument('--map_coverage_pct', type=float, default=0.0,
    help='Percent of fragment-start bins used as mappability Gaussian centers.')
parser.add_argument('--map_sigma', type=float, default=5.0,
    help='Standard deviation for mappability Gaussian bias.')
parser.add_argument('--map_enrich', type=float, default=1.0,
    help='Enrichment factor for mappability Gaussian bias.')
parser.add_argument('--map_exp', type=float, default=1.0,
    help='Exponent to reshape mappability PMF (<1 flattens, >1 sharpens).')
parser.add_argument('--seed', type=int, default=42,
    help='Fallback random seed for reproducible simulation components.')
parser.add_argument('--tf_seed', type=int, default=None,
    help='Random seed for reproducible TF peak placement. Defaults to --seed.')
parser.add_argument('--map_seed', type=int, default=None,
    help='Random seed for reproducible mappability center placement. Defaults to --seed.')
parser.add_argument('--read_seed', type=int, default=None,
    help='Random seed for read-count sampling. Defaults to --seed.')
parser.add_argument('--nb_k', type=float, default=10.0,
    help='Inverse-dispersion (size) parameter k for negative-binomial noise; smaller k ⇒ more variance.')
parser.add_argument('--output_fasta1', required=True,
    help='Path to write mate 1 (R1) reads in FASTA format.')
parser.add_argument('--output_fasta2', required=True,
    help='Path to write mate 2 (R2) reads in FASTA format.')
parser.add_argument('--pmf_csv', type=str, default=None,
    help='Path to CSV file storing PMF and variance per bin')
parser.add_argument('--skip_pmf_csv', action='store_true',
    help='Do not write a PMF CSV, even when FASTA input is provided.')
parser.add_argument('--planted_peaks_bed', type=str, default=None,
    help='Path to BED file storing planted peak centers (1-bp intervals).')
parser.add_argument('--planted_peak_intervals_bed', type=str, default=None,
    help='Path to BED file storing planted peak support intervals.')

args, _ = parser.parse_known_args()

# this particular genome has 3 chroms, chr1, chr2, chr3, with lengths 100, 200, 300 respectively
fasta = args.fasta
coverage = args.coverage
tf_peak_count = args.tf_peak_count
k = args.fragment_length
tf_sigma = args.tf_sigma
tf_enrichment = args.tf_enrichment
accessibility_bed = args.accessibility_bed
acc_weight = args.acc_weight
gc_bias_params = args.gc_bias_params
tf_exp = args.tf_exp
gc_exp = args.gc_exp
acc_exp = args.acc_exp
map_coverage_pct = args.map_coverage_pct
map_sigma = args.map_sigma
map_enrich = args.map_enrich
map_exp = args.map_exp
seed = args.seed
tf_seed = args.tf_seed
map_seed = args.map_seed
read_seed = args.read_seed
read_length = args.read_length
nb_k = args.nb_k
output_fasta1 = args.output_fasta1
output_fasta2 = args.output_fasta2
pmf_csv = args.pmf_csv
skip_pmf_csv = args.skip_pmf_csv
planted_peaks_bed = args.planted_peaks_bed
planted_peak_intervals_bed = args.planted_peak_intervals_bed

if not skip_pmf_csv and not pmf_csv and fasta:
    base = os.path.splitext(os.path.basename(fasta))[0]
    pmf_csv = f'{base}_pmf.csv'

"""PMF CSV Structure
___________________

- bin_idx: zero-based index of fragment start bin
- pmf: probability assigned to the bin
- variance: pmf * (1 - pmf)
"""

if seed is not None:
    random.seed(seed)
    np.random.seed(seed)



"""
Functions
_________
"""


def resolve_component_seed(component_seed: int, fallback_seed: int) -> int:
    """Return explicit component seed or the shared fallback seed."""
    return fallback_seed if component_seed is None else component_seed


def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    table = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(table)[::-1]

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

def build_gc_bias_pmf(sequence: str, loess_params: Dict, fragment_length: int,
                      exp: float = 1.0) -> np.ndarray:
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
    s = bias.sum()
    if s == 0:
        bias = np.ones_like(bias, dtype=float)
        s = bias.sum()
    bias /= s
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


def build_mappability_bias_pmf(length: int, map_coverage_pct: float,
                               map_sigma: float, map_enrich: float,
                               exp: float, rng: np.random.Generator) -> np.ndarray:
    """Return mappability PMF with deterministic center count from coverage percent."""
    bias = np.ones(length, dtype=float)
    if length <= 0:
        return bias

    pct = min(max(float(map_coverage_pct), 0.0), 100.0)
    num_map_peaks = int(round(length * pct / 100.0))
    num_map_peaks = min(max(num_map_peaks, 0), length)

    if num_map_peaks == 0:
        bias /= bias.sum()
        if exp != 1.0:
            bias = np.power(bias, exp)
            bias /= bias.sum()
        return bias
    if map_sigma <= 0:
        raise ValueError('map_sigma must be > 0 when map_coverage_pct > 0')

    centers = rng.choice(length, size=num_map_peaks, replace=False)
    positions = np.arange(length)
    for center in centers:
        kernel = norm.pdf(positions, loc=center, scale=map_sigma)
        bias += map_enrich * kernel

    bias /= bias.sum()
    if exp != 1.0:
        bias = np.power(bias, exp)
        bias /= bias.sum()
    return bias

def create_pmf(chrom_len: int, k: int) -> List[float]:
    """Initialize uniform PMF array for one chromosome."""

    num_bins = chrom_len - k + 1

    pmf = [1] * num_bins

    return pmf

def create_pmf_all_chroms(
    fasta: str,
    fragment_length: int,
    tf_peak_count: int,
    tf_sigma: float,
    tf_enrichment: float,
    accessibility_bed: str,
    acc_weight: float,
    seed: int,
    tf_seed: int,
    gc_bias_params: str,
    tf_exp: float,
    gc_exp: float,
    acc_exp: float,
    map_coverage_pct: float,
    map_sigma: float,
    map_enrich: float,
    map_exp: float,
    map_seed: int,
) -> Tuple[Dict[str, List[float]], Dict[str, List[int]]]:
    """Build PMF dictionary for all chromosomes with bias modeling."""

    genome_pmfs = {}
    planted_peaks = {}
    tf_rng = np.random.default_rng(resolve_component_seed(tf_seed, seed))
    map_rng = np.random.default_rng(resolve_component_seed(map_seed, seed))
    gc_params = {'csv': gc_bias_params}
    for chrom_id, seq in lib.read_fasta(fasta):
        if len(seq) < fragment_length:
            continue
        base = np.array(
            create_pmf(
                len(seq),
                fragment_length,
            ),
            dtype=float,
        )
        length = base.shape[0]
        tf_centers = (
            tf_rng.integers(0, length, size=tf_peak_count)
            if tf_peak_count > 0
            else np.array([], dtype=int)
        )
        planted_peaks[chrom_id] = tf_centers.astype(int).tolist()
        tf_bias = build_tf_bias_pmf(
            length, tf_centers.tolist(), tf_sigma, tf_enrichment, exp=tf_exp
        )
        gc_bias = build_gc_bias_pmf(seq, gc_params, fragment_length, exp=gc_exp)
        acc_bias = build_accessibility_bias_pmf(length, accessibility_bed, acc_weight,
                                               chrom_id.split()[0], exp=acc_exp)
        map_bias = build_mappability_bias_pmf(
            length, map_coverage_pct, map_sigma, map_enrich, map_exp, map_rng
        )
        combined = base * tf_bias * gc_bias * acc_bias * map_bias
        pmf = combined / combined.sum()
        genome_pmfs[chrom_id] = pmf.tolist()
    return genome_pmfs, planted_peaks


def write_pmf_csv(genome_pmfs: Dict[str, List[float]], path: str) -> None:
    """Write PMF and variance arrays to CSV, including chromosome."""
    rows = []
    for chrom_id, pmf in genome_pmfs.items():
        arr = np.asarray(pmf, dtype=float)
        var = arr * (1 - arr)
        for idx, (p, v) in enumerate(zip(arr, var)):
            rows.append((chrom_id, idx, p, v))
    df = pd.DataFrame(rows, columns=['chrom', 'bin_idx', 'pmf', 'variance'])
    df.to_csv(path, index=False)


def write_planted_peaks_bed(planted_peaks: Dict[str, List[int]], path: str) -> None:
    """Write planted peak centers as a 1-bp BED file."""
    peak_idx = 1
    with open(path, 'w') as fh:
        for chrom_id, centers in planted_peaks.items():
            chrom = chrom_id.split()[0]
            for center in centers:
                start = int(center)
                end = start + 1
                fh.write(f"{chrom}\t{start}\t{end}\tpeak_{peak_idx}\n")
                peak_idx += 1


def write_planted_peak_intervals_bed(
    planted_peaks: Dict[str, List[int]],
    path: str,
    sigma: float,
    fragment_length: int,
    fasta: str,
) -> None:
    """Write planted support intervals for interval-aware broad evaluation."""
    chrom_lengths = {
        chrom_id.split()[0]: len(seq)
        for chrom_id, seq in lib.read_fasta(fasta)
    }
    half_width = max(int(round(3 * sigma)), fragment_length // 2)
    peak_idx = 1
    with open(path, 'w') as fh:
        for chrom_id, centers in planted_peaks.items():
            chrom = chrom_id.split()[0]
            chrom_length = chrom_lengths.get(chrom, 0)
            for center in centers:
                start = max(int(center) - half_width, 0)
                end = min(int(center) + half_width + 1, chrom_length)
                fh.write(f"{chrom}\t{start}\t{end}\tpeak_interval_{peak_idx}\n")
                peak_idx += 1

def sample_genome(
    fasta: str,
    genome_pmfs: Dict[str, List[float]],
    coverage: float,
    fragment_length: int,
    read_length: int,
    nb_k: float,
    seed: int,
) -> (List[tuple], Dict[str, np.ndarray]):
    """Return paired-end reads and negative-binomial counts."""

    chrom_bias = {}
    seqs = {}
    read_rng = np.random.default_rng(resolve_component_seed(read_seed, seed))
    total_bp = 0
    for chrom_id, seq in lib.read_fasta(fasta):
        seqs[chrom_id] = seq
        if len(seq) < fragment_length:
            continue
        total_bp += len(seq)
        chrom_bias[chrom_id] = len(seq)
    for cid in chrom_bias:
        chrom_bias[cid] /= total_bp

    total_reads = int((total_bp * coverage) / fragment_length)

    paired_reads = []
    nb_counts_dict = {}

    for chrom_id, pmf_list in genome_pmfs.items():
        if chrom_id not in chrom_bias:
            continue
        pmf = np.array(pmf_list, dtype=float)
        expected_counts = pmf * (chrom_bias[chrom_id] * total_reads)
        nb_counts = read_rng.negative_binomial(
            n=nb_k,
            p=nb_k / (nb_k + expected_counts)
        )
        nb_counts_dict[chrom_id] = nb_counts
        seq = seqs[chrom_id]
        for start_idx, count in enumerate(nb_counts):
            for _ in range(int(count)):
                frag_start = start_idx
                frag_end = frag_start + fragment_length - 1
                r1_seq = seq[frag_start: frag_start + read_length]
                r2_seq = reverse_complement(
                    seq[frag_end - read_length + 1: frag_end + 1]
                )
                paired_reads.append((r1_seq, r2_seq))

    return paired_reads, nb_counts_dict

def write_paired_fasta(paired_reads, output_path):
    """Write paired-end reads in interleaved FASTA format."""

    with open(output_path, 'w') as fh:
        for i, (r1, r2) in enumerate(paired_reads, start=1):
            fh.write(f">read_{i:06d}/1\n{r1}\n")
            fh.write(f">read_{i:06d}/2\n{r2}\n")


def write_r1_r2_fastas(paired_reads, r1_path, r2_path):
    """Write paired-end reads to two separate FASTA files: R1 and R2."""
    with open(r1_path, 'w') as f1, open(r2_path, 'w') as f2:
        for i, (r1, r2) in enumerate(paired_reads, start=1):
            f1.write(f">read_{i:06d}/1\n{r1}\n")
            f2.write(f">read_{i:06d}/2\n{r2}\n")


def ensure_parent_dir(path: str) -> None:
    """Create parent directory for path when needed."""
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)

"""Below code will print the FASTA for the reads generated from experiment"""

if fasta:
    if read_length > k:
        raise ValueError('read_length must not exceed fragment_length')
    genome_pmf, planted_peaks = create_pmf_all_chroms(
        fasta,
        k,
        tf_peak_count,
        tf_sigma,
        tf_enrichment,
        accessibility_bed,
        acc_weight,
        seed,
        tf_seed,
        gc_bias_params,
        tf_exp,
        gc_exp,
        acc_exp,
        map_coverage_pct,
        map_sigma,
        map_enrich,
        map_exp,
        map_seed,
    )

    paired_reads, nb_counts = sample_genome(
        fasta,
        genome_pmf,
        coverage,
        k,
        read_length,
        nb_k,
        seed,
    )
    ensure_parent_dir(output_fasta1)
    ensure_parent_dir(output_fasta2)
    write_r1_r2_fastas(paired_reads, output_fasta1, output_fasta2)
    
    if pmf_csv and not skip_pmf_csv:
        ensure_parent_dir(pmf_csv)
        write_pmf_csv(genome_pmf, pmf_csv)
    if planted_peaks_bed:
        ensure_parent_dir(planted_peaks_bed)
        write_planted_peaks_bed(planted_peaks, planted_peaks_bed)
    if planted_peak_intervals_bed:
        ensure_parent_dir(planted_peak_intervals_bed)
        write_planted_peak_intervals_bed(
            planted_peaks,
            planted_peak_intervals_bed,
            tf_sigma,
            k,
            fasta,
        )

