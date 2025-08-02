"""
This program models a Chip-Seq experiment as a sampler drawing from a combined probability mass function
of the background and foreground

This program will print FASTA of reads generated based on user defined arguments (coverage, number of peaks in fg, bg, etc)
"""

"""Imports"""

import random
import os
from typing import List, Dict

import numpy as np
import pandas as pd
from scipy.stats import norm

from . import lib
import argparse


"""
Assumptions for our Model
_________________________

- We will just make the assumption that all fragments are uniform in this model
- CHIP-seq reads are anywhere from 100-500bp approx.
"""


"""
Global Variables
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
parser.add_argument('--tf_enrichment', type=float, default=1.0,
    help='Enrichment factor for TF-binding bias')
parser.add_argument('--accessibility_bed', type=str, default=None,
    help='BED file describing open chromatin intervals')
parser.add_argument('--acc_weight', type=float, default=1.0,
    help='Weight multiplier for accessible regions')
parser.add_argument('--gc_bias_params', type=str, default=None,
    help='CSV file for GC bias lookup table')
parser.add_argument('--seed', type=int, default=42,
    help='Random seed for reproducible TF peak placement')
parser.add_argument('--nb_k', type=float, default=10.0,
    help='Inverse-dispersion (size) parameter k for negative-binomial noise; smaller k ⇒ more variance.')
parser.add_argument('--output_fasta', required=True,
    help='Path to write paired-end reads in FASTA format.')
parser.add_argument('--pmf_csv', type=str, default=None,
    help='Path to CSV file storing PMF and variance per bin')

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
seed = args.seed
read_length = args.read_length
nb_k = args.nb_k
output_fasta = args.output_fasta
pmf_csv = args.pmf_csv

if not pmf_csv and fasta:
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


'''
Notes for random_genome_1.fa.gz:
______________________

- total bp is 600
- with coverage =1, and k=20, we should see 30 samples getting pulled (CHECKED)
- also the chrom frequencies should be 1/6, 2/6, 3/6
    - meaning when we look at percentage of reads generated with respect to which chromosomes...
    - we should observe that it converges to above probabilities (CHECKED)

- to check if peaks are added properly we can just graph them (CHECKED)
- we should also check if experiment samples have similar distributions to experiment pmf (CHECKED)
'''


"""
Functions
_________
"""


def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    table = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(table)[::-1]




def build_tf_bias_pmf(length: int, peaks: List[int], sigma: float,
                      enrichment: float) -> np.ndarray:
    """Return normalized TF-binding bias PMF."""
    bias = np.ones(length, dtype=float)
    positions = np.arange(length)
    for p in peaks:
        kernel = norm.pdf(positions, loc=p, scale=sigma)
        bias += enrichment * kernel
    bias /= bias.sum()
    return bias

def build_gc_bias_pmf(sequence: str, loess_params: Dict, fragment_length: int) -> np.ndarray:
    """Return normalized GC-content bias PMF."""
    if not loess_params or 'csv' not in loess_params or loess_params['csv'] is None:
        bias = np.ones(len(sequence) - fragment_length + 1, dtype=float)
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
    return bias

def build_accessibility_bias_pmf(length: int, accessibility_bed: str,
                                 acc_weight: float) -> np.ndarray:
    """Return normalized chromatin accessibility bias PMF."""
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
                start = int(fields[1])
                end = int(fields[2])
                start = max(start, 0)
                end = min(end, length)
                bias[start:end] *= acc_weight
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
    gc_bias_params: str,
) -> Dict[str, List[float]]:
    """Build PMF dictionary for all chromosomes with bias modeling."""

    genome_pmfs = {}
    rng = np.random.default_rng(seed)
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
        tf_centers = rng.integers(0, length, size=max(1, tf_peak_count))
        tf_bias = build_tf_bias_pmf(length, tf_centers.tolist(), tf_sigma, tf_enrichment)
        gc_bias = build_gc_bias_pmf(seq, gc_params, fragment_length)
        acc_bias = build_accessibility_bias_pmf(length, accessibility_bed, acc_weight)
        combined = base * tf_bias * gc_bias * acc_bias
        pmf = combined / combined.sum()
        genome_pmfs[chrom_id] = pmf.tolist()
    return genome_pmfs


def write_pmf_csv(genome_pmfs: Dict[str, List[float]], path: str) -> None:
    """Write PMF and variance arrays to CSV."""
    rows = []
    for pmf in genome_pmfs.values():
        arr = np.asarray(pmf, dtype=float)
        var = arr * (1 - arr)
        for idx, (p, v) in enumerate(zip(arr, var)):
            rows.append((idx, p, v))
    df = pd.DataFrame(rows, columns=['bin_idx', 'pmf', 'variance'])
    df.to_csv(path, index=False)

def sample_genome(
    fasta: str,
    genome_pmfs: Dict[str, List[float]],
    coverage: float,
    fragment_length: int,
    read_length: int,
    nb_k: float,
) -> (List[tuple], Dict[str, np.ndarray]):
    """Return paired-end reads and negative-binomial counts."""

    chrom_bias = {}
    seqs = {}
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
        nb_counts = np.random.negative_binomial(
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

"""Below code will print the FASTA for the reads generated from experiment"""

if fasta:
    if read_length > k:
        raise ValueError('read_length must not exceed fragment_length')
    genome_pmf = create_pmf_all_chroms(
        fasta,
        k,
        tf_peak_count,
        tf_sigma,
        tf_enrichment,
        accessibility_bed,
        acc_weight,
        seed,
        gc_bias_params,
    )

    paired_reads, nb_counts = sample_genome(
        fasta,
        genome_pmf,
        coverage,
        k,
        read_length,
        nb_k,
    )
    write_paired_fasta(paired_reads, output_fasta)
    
    if pmf_csv:
        write_pmf_csv(genome_pmf, pmf_csv)

    '''
    uncomment both to compare pmf graph with actual experiment graph
    NOTE: x-axis for pmf is not base coordinate for all bases in genome, rather base coordinates for first base of kmer/fragment
    '''

    # graph_all_genome_pmf(genome_pmf)
    # graph_experiment(exp, genome_pmf)









