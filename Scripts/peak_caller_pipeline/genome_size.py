"""
Find Genome Size Sum
"""

import LIB
import sys
import gzip

genome = sys.argv[1]

'''
total = 0
with open gzip.open(genome, 'rt') as f:
    for line in f:
        if line.startswith('>'):  continue
        total += len(line)

print(total)
'''
total = 0
for id, seq in LIB.read_fasta(genome):
    total += len(seq)

print(total)
