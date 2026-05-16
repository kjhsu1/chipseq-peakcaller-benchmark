"""
This program models a Chip-Seq experiment as a sampler drawing from a combined probability mass function
of the background and foreground

This program will print FASTA of reads generated based on user defined arguments (coverage, number of peaks in fg, bg, etc)
"""

import sys
import random
import matplotlib.pyplot as plt
import LIB
import json
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

parser = argparse.ArgumentParser(description="CHIP-seq Experiment Simulator")
parser.add_argument('--fasta', type=str, required=True,
    help='Path to the FASTA file')
parser.add_argument('--coverage', type=float, default=1.0,
    help='Coverage depth to simulate')
parser.add_argument('--num_bg_peaks', type=int, default=1,
    help='Number of random background peaks to add')
parser.add_argument('--num_fg_peaks', type=int, default=1,
    help='Number of random foreground peaks to add')
parser.add_argument('--fragment_length', type=int, default=20,
    help='The fragment length of each read (k-mer size); usually between 100=500bp')
parser.add_argument('--peak_broadness', type=int, default=9,
    help='Number of bins that make up each peak')
parser.add_argument('--tallness', type=int, default=10,
    help='Height multiplier for peak center relative to baseline')

args = parser.parse_args()

# this particular genome has 3 chroms, chr1, chr2, chr3, with lengths 100, 200, 300 respectively
fasta = args.fasta
coverage = args.coverage
num_bg_peaks = args.num_bg_peaks
num_fg_peaks = args.num_fg_peaks
k = args.fragment_length
peak_broadness = args.peak_broadness
tallness = args.tallness


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


def add_peaks(pmf, num_peaks, peak_broadness, tallness):
    """Add num_peaks peaks using broadness and tallness"""
    center = max(peak_broadness // 2, 1)
    peaks = []
    for i in range(peak_broadness):
        distance = abs(i - center)
        height = 1 + (tallness - 1) * (1 - distance / center)
        peaks.append(height)
    used_peaks = []
    p_index = random.randint(0, len(pmf) - len(peaks))

    for peak in range(num_peaks):
        while p_index in used_peaks:
            p_index = random.randint(0, len(pmf) - len(peaks)) # can't create peaks on the last 6 bases
        for i in range(p_index, p_index + len(peaks)): # create the peaks in the bin
            pmf[i] = pmf[i] + peaks[i-p_index] # should we add or multiply?
        # don't want peaks to overlap, and want a space of at least 1 between peaks
        start = p_index - 1
        end = p_index + len(peaks) + 1
        used_peaks.extend(range(start, end))
    return pmf 

def normalize_bins(pmf):
    """
    Takes in  a bins pmf, normalize so all bins sum to 1
    """
    total = sum(pmf)
    for i, bin in enumerate(pmf):
        pmf[i] = (pmf[i]/total)
    return pmf

def sample_from_bins(pmf, num_sample):
    """
    Takes in a bins pmf, and samples num_sample times from the distribution.
    """
    bins_indices = list(range(len(pmf))) # ex. range(3) = [0, 1, 2]
    samples = random.choices(bins_indices, weights=pmf, k=num_sample)
    return samples

def create_pmf(chrom_len, num_bg_peaks, num_fg_peaks, k, peak_broadness, tallness):
    """Create combined pmf for one chromosome"""
    
    num_bins = chrom_len - k + 1

    # initialize bins 
    pmf = [1] * num_bins

    # add background peaks
    pmf = add_peaks(pmf, num_bg_peaks, peak_broadness, tallness)

    # add foreground peaks
    pmf = add_peaks(pmf, num_fg_peaks, peak_broadness, tallness)

    pmf = normalize_bins(pmf)

    return pmf

def create_pmf_all_chroms(fasta, peak_broadness, tallness):
    """Build pmf dictionary for all chromosomes"""
    genome_pmfs = {}
    for id, seq in LIB.read_fasta(fasta):
        if len(seq) < k: continue
        genome_pmfs[id] = create_pmf(len(seq), num_bg_peaks, num_fg_peaks,
            k, peak_broadness, tallness)
    # print(json.dumps(genome_pmfs, indent=4))
    return genome_pmfs

def chrom_bias(fasta):
    """
    Find the sampling bias for each chrom based on length relative to the sum total of genomic bps
    """
    chrom_bias = {}
    total_bp = 0
    for id, seq in LIB.read_fasta(fasta):
        total_bp += len(seq)
        chrom_bias[id] = len(seq)
    for id in chrom_bias:
        chrom_bias[id] = chrom_bias[id] / total_bp
    return chrom_bias

def sample_genome(fasta, genome_pmfs):
    """
    Sample from a genome pmf
    
    Parameters: Path to fasta, genomic P.M.F

    Output: dictionary of reads
        key: chrom (ex. 'chr1', 'chr2')
        value: list of lists of reads, (ex. [ [1,2,3,4], [4,5,6,7] ])
            - in the example, it denotes two reads; both 4 bps long spanning those coords 
    """
    # find chrom bias
    chrom_bias = {}
    total_bp = 0
    for id, seq in LIB.read_fasta(fasta):
        if len(seq) < k: continue
        total_bp += len(seq)
        chrom_bias[id] = len(seq)
    for id in chrom_bias:
        chrom_bias[id] = chrom_bias[id] / total_bp
    
    # Do the experiment
    num_sample = int((total_bp * coverage) / k)
    chroms = list(chrom_bias.keys())
    biases = list(chrom_bias.values())
    reads_dict = {} # store samples for each chrom
    for chrom in chroms:
        reads_dict[chrom] = []

    for i in range(num_sample):
        picked_chrom = random.choices(chroms, weights=biases)[0]
        sample_index = sample_from_bins(genome_pmfs[picked_chrom], 1)[0]

        strand = '+'
        coords = list(range(sample_index, sample_index + k))
        reads_dict[picked_chrom].append((strand, coords))

    # print(json.dumps(reads_dict, indent=4))
    return reads_dict

def sample_to_fasta(exp, fasta):
    """Convert sampled indices into FASTA format"""

    frag_num = 1

    for id, seq in LIB.read_fasta(fasta):
        for strand, coords in exp[id]:
            dna = [seq[index] for index in coords]

            print(f'>read_{frag_num}:{id}:{strand}')
            print(''.join(dna))
            frag_num += 1

"""
Functions for Debugging/Checking
_______________________
"""

def total_num_reads(exp_sample):
    """
    Get total number of reads from a sampling experiment
    """
    read_sum = 0
    for chrom in exp_sample:
        read_sum += len(exp_sample[chrom])
    return read_sum

def chrom_distribution(exp_sample):
    """
    Get distribution of chromosomes in a sampling experiment

    WILL NOT RETURN VALUE, WILL JUST PRINT
    """
    total = total_num_reads(exp_sample)
    dict = {}
    for chrom in exp_sample:
        chrom_total = len(exp_sample[chrom])
        dict[chrom] = chrom_total / total
    print('Experimental Distribution', '\n',json.dumps(dict, indent=4))
    
    # can compare with expected distribution
    c_bias = chrom_bias(fasta)
    print('Expected Distribution', '\n', json.dumps(c_bias, indent=4))

def graph_pmf(pmf, title="Experimental P.M.F"):
    """
    Graphs the P.M.F
    """
    bin_coords = range(len(pmf))
    probs = pmf

    plt.figure()
    plt.plot(bin_coords, probs)
    plt.xlabel('Read/K-mer index')
    plt.ylabel('Probability')
    plt.title(title)
    plt.show()

def graph_all_genome_pmf(genome_pmf):
    for chrom in genome_pmf.keys():
        pmf = genome_pmf[chrom]
        graph_pmf(pmf, title=f'P.M.F for Experimental, {chrom}')

def graph_experiment(exp, genome_pmf):
    """
    Graph the sample base distribution
    """
    # initialize dict to track number of samples
    sample_distribution = {}
    for chrom in genome_pmf.keys():
        sample_distribution[chrom] = {}
        for i in range(len(genome_pmf[chrom]) + k-1):
            sample_distribution[chrom][i] = 0

    # print(json.dumps(sample_distribution, indent=4)) # debug

    # count each sample
    for chrom in exp.keys():
        samples = exp[chrom]
        for sample in samples:
            for base in sample:
                sample_distribution[chrom][base] += 1

    # print(json.dumps(sample_distribution, indent=4)) # debug
    
    # graph each chrom
    for chrom in sample_distribution.keys():
        x = list(sample_distribution[chrom].keys())
        y = list(sample_distribution[chrom].values())
        # print(x) # debug
        # print(y) #debug
        plt.figure()
        plt.plot(x, y)
        plt.xlabel('Base Coordinate')
        plt.ylabel('Counts')
        plt.title(f'Base Distribution in Sample: {chrom}')
        plt.show()

'''
Below code will print the FASTA for the reads generated from experiment
'''

genome_pmf = create_pmf_all_chroms(fasta, peak_broadness, tallness)
exp = sample_genome(fasta, genome_pmf)
# print(json.dumps(genome_pmf, indent=4)) # for debugging
# print(json.dumps(exp, indent=4)) # for debugging
sample_to_fasta(exp, fasta)

'''
uncomment both to compare pmf graph with actual experiment graph
NOTE: x-axis for pmf is not base coordinate for all bases in genome, rather base coordinates for first base of kmer/fragment
'''

# graph_all_genome_pmf(genome_pmf)
# graph_experiment(exp, genome_pmf)









