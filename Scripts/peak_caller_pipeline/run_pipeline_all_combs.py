"""
This Program Will Call Peaks On All Combinations of Control and Experiment Reads FASTAs

Main Input: 2 Directories: 
    - One which stores all control reads FASTAs of various coverages,
    - Second which stores all experiment reads FASTAs of various coverages

Output: All of the peak files for all combinations

ex. control directory contains: control_5X.fa, control_10X.fa, control_15X.fa
    experiment directory contains: 5X.fa, 10X.fa, 15X.fa

    Then program will return peak files of 3*3 = 9 total peak files

NOTE
******PLEASE READ TOP COMMENT IN PIPELINE.PY BEFORE RUNNING ******
"""


"""
Imports
"""

import subprocess
import argparse
import os

"""
User Arguments
"""

parser = argparse.ArgumentParser(description='Call Peaks On All Combinations of Control+Experiment')

parser.add_argument('--control_dir', type=str, required=True,
                    help='path to directory containing all control reads FASTA files')
parser.add_argument('--experiment_dir', type=str, required=True,
                    help='path to directory containing all experiment reads FASTA files')
parser.add_argument('--genome_index_path', type=str, required=True,
                    help='path to genome index file')
parser.add_argument('--exp_alignment_output_path', type=str, required=True,
                    help='path to directory where you want experiment SAM/BAM/etc stored')
parser.add_argument('--control_alignment_output_path', type=str, required=True,
                    help='path to directory where you want control SAM/BAM/etc stored')
parser.add_argument('--peaks_called_directory_path', type=str, required=True,
                    help='path to directory where you want peak files to be stored')
parser.add_argument('--genome_size', type=str, required=True,
                    help='total genome size')

arg = parser.parse_args()

control_dir = arg.control_dir
experiment_dir = arg.experiment_dir
genome_index_path            = arg.genome_index_path
exp_alignment_output_path    = arg.exp_alignment_output_path
control_alignment_output_path = arg.control_alignment_output_path
peaks_called_directory_path  = arg.peaks_called_directory_path
genome_size                  = arg.genome_size

"""
Functions
"""

def get_files_in_dir(directory):
    """
    Return List of Paths of All Files in Given Directory
    """
    #directory = '../../../Shared/Reads_FASTA/control_unbiased_reads'
    filenames = os.listdir(directory)

    paths = []
    for file in filenames:
        paths.append(os.path.join(directory, file))
    return paths

def call_pipeline(
        exp_fasta,
        control_fasta,
        genome_index,
        exp_output_path,
        control_output_path,
        peaks_output_path,
        genome_size):
    """
    Calls the end-to-end CHIP-seq pipeline script in a single subprocess.

    Parameters
    ----------
    - exp_fasta        : Path to experiment / treatment reads FASTA
    - control_fasta    : Path to control / background reads FASTA
    - genome_index     : Path (prefix only) to the Bowtie2 genome index
    - exp_output_path  : Directory that will hold SAM/BAMs for the experiment
    - control_output_path : Directory that will hold SAM/BAMs for the control
    - peaks_output_path   : Directory where MACS2 peak files should be written
    - genome_size      : Effective genome size string passed through to MACS2
    """

    cmd = [
        "python3", 'pipeline.py',
        "--exp_fasta", exp_fasta,
        "--control_fasta", control_fasta,
        "--genome_index", genome_index,
        "--exp_output_path", exp_output_path,
        "--control_output_path", control_output_path,
        "--peaks_output_path", peaks_output_path,
        "--genome_size", genome_size
    ]

    try:
        subprocess.run(cmd, check=True)
        print("CHIP-seq pipeline finished successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Pipeline run failed: {e}")
    
def peak_call_all_combs(control_dir, 
                        experiment_dir, 
                        genome_index_path, 
                        exp_alignment_output_path, 
                        control_alignment_output_path, 
                        peaks_called_directory_path, 
                        genome_size):
    
    """
    Takes in Control and Experiment Directories with CHIP-Seq Reads and Calls Peaks On All Possible Combinations
    """

    all_control = get_files_in_dir(control_dir)
    all_exp = get_files_in_dir(experiment_dir)

    for control in all_control:
        for exp in all_exp:
            call_pipeline(
                exp,
                control,
                genome_index_path,
                exp_alignment_output_path,
                control_alignment_output_path,
                peaks_called_directory_path,
                genome_size
            )

# run
peak_call_all_combs(control_dir, 
                    experiment_dir, 
                    genome_index_path, 
                    exp_alignment_output_path, 
                    control_alignment_output_path, 
                    peaks_called_directory_path, 
                    genome_size)





