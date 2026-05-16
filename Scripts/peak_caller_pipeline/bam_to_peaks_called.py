"""
Calls Peaks On Sorted BAM File Using MACS2
"""

"""
Imports
_______
"""

import os
import argparse
import subprocess

"""
User Input Arguments
__________
"""

parser = argparse.ArgumentParser(description='Calling Peaks On Sorted BAM File Using MACS2')
parser.add_argument('--exp_sorted_bam_path', type=str, required=True,
                    help='Path to Sorted BAM File For Experiment/Treatment')
parser.add_argument('--control_sorted_bam_path', type=str, required=True,
                    help='Path to Sorted BAM File For Control/Background')
parser.add_argument('--genome_size', type=str, required=True, 
                    help='Total Genome Size; Sum # of Bases of All Chromosomes')
parser.add_argument('--output_dir', type=str, required=True,
                    help='Path to Output Directory; Do Not Need File Name')
parser.add_argument('--name', type=str, required=True, help=
                    'file prefix you want to use for output peaks file')
arg = parser.parse_args()

exp_sorted_bam_path = arg.exp_sorted_bam_path
control_sorted_bam_path = arg.control_sorted_bam_path
genome_size = arg.genome_size
output_dir = arg.output_dir
file_prefix = arg.name

"""
Functions
_______
"""

def bam_to_peaks(exp_sorted_bam_path, control_sorted_bam_path, genome_size, output_dir, file_prefix):
    '''
    basename = os.path.basename(exp_sorted_bam_path) # ex. sorted.exp1.bam
    file_prefix, _ = os.path.splitext(basename) # sorted.exp1, bam
    split = file_prefix.split('.')
    file_prefix = split[1]
    '''

    cmd = [
        'macs2',
        'callpeak',
        '-t', exp_sorted_bam_path, # treatment/experiment 
        '-c', control_sorted_bam_path, # control/background
        '-f', 'BAM', # file type
        '-g', genome_size,
        '-n', file_prefix,
        '--outdir', output_dir, # output directory 
        '--nomodel', # COMMENT THIS WHEN RUNNING BIG ENOUGH GENOME + READ NUMBER
        '--extsize', '300' # COMMENT THIS WHEN RUNNING BIG ENOUGH GENOME + READ NUMBER
    ]

    try:
        subprocess.run(cmd, check=True)
        print('Sucessfully Called Peaks')
    except subprocess.CalledProcessError as e:
        print(f'Unsuccessful Peak Calling: {e}')


"""
Main
____
"""

def main():
    bam_to_peaks(exp_sorted_bam_path, control_sorted_bam_path, genome_size, output_dir, file_prefix)

main()