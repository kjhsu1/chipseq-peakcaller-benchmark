"""
Program Will Run chip_seq.py On Given List of Coverages 

NOTE
- program takes command line arguments for genome path and output dir 
- however, list of coverages are hard to provide in command line. 
    - therefore have to directly change this script to alter coverages, # of peaks, fragment length
    - edit the cmd variable in the run_all_covs function below
"""

"""
Imports
"""

import argparse
import subprocess
import os

"""
User Arguments
"""

parser = argparse.ArgumentParser(description='Program Will Run a Chip-Seq Experiment On All Coverages Specified')

parser.add_argument('--genome_path', type=str, required=True, help=
                    'path to gzipped genome file')
parser.add_argument('--output', type=str, required=True, help=
                    'path to the output directory for Reads FASTA')

arg = parser.parse_args()

genome_path = arg.genome_path
output_path = arg.output
covs = [1, 5, 10, 25] # list of coverages; EDIT THIS IF YOU WANT TO RUN DIFFERENT COVERAGES

"""
Functions
"""

def run_all_covs(genome_path, covs, output_path):
    for cov in covs:
        file_path = os.path.join(output_path, f'cov_{cov}_chip_seq_reads.fa')
        
        cmd = [
            "python3", "chip_seq.py",
            "--fasta", genome_path,
            "--coverage", str(cov),
            "--num_bg_peaks", "0",
            "--num_fg_peaks", "1",
            "--fragment_length", "300"
        ]

        try:
            with open(file_path, 'w') as out:
                subprocess.run(cmd, check=True, stdout=out)
            print(f'Successfully Ran Coverage: {cov}')
        except subprocess.CalledProcessError:
            print(f'Failed At Coverage: {cov}')

# run 
run_all_covs(genome_path, covs, output_path)