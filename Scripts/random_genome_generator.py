import random
import json

"""
create a random genome
"""

# number of elements=num chromosomes, value=length of chrom
genome_specs = [250000, 250000, 250000, 250000]

def generator(genome_specs):
    """
    Create dictionary with specified # of chromosomes and its bp_lengths

    Input: genome_specs; list with bp_length specifications for each chromosome (ex. [5000, 4000, 9000] which is a 3 chrom genome)

    Output: prints fasta with all chromosomes
    """

    # create dictionary to store generated chroms 
    genome = {}
            
    # randomly generate chroms
    for i, chrom in enumerate(genome_specs):
        chrom_num = i + 1
        genome[chrom_num] = ''.join(random.choices('ATGC', k=chrom))
    
    for chrom in genome:
        print(f'>chr{chrom}')
        print(genome[chrom])

# print the fasta
genome = generator(genome_specs)

