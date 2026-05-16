"""
Convert a single SAM file to BAM, sorted BAM, sorted indexed BAM
"""

"""
Imports
_______
"""

import argparse
import subprocess
import os

"""
User Input Arguments
_________
"""

parser = argparse.ArgumentParser(description="Convert SAM to BAM, sorted BAM, and sorted indexed BAM")
parser.add_argument("--sam_path", type=str, required=True, 
    help="Path to SAM file")
parser.add_argument("--bam_output", type=str, required=True, 
    help="Path to BAM file; Make Sure Basename ends with BAM Filename; \
        NOTE: Sorted BAM and Sorted Indexed BAM Will Also Be Placed in the Same Directory")
arg = parser.parse_args()

sam_path = arg.sam_path
output_path = arg.bam_output # Newly created BAM file path (including the file name itself)

def sam_to_bam(sam_path, output_path):
    """
    Converts SAM to BAM
    """

    cmd = [
        "samtools",
        "view",
        "-Sb",
        sam_path,
        "-o",
        output_path
    ]

    try: 
        subprocess.run(cmd, check=True)
        print('Sucessful SAM to BAM conversion')
    except subprocess.CalledProcessError as e:
        print(f'Unsucessful SAM to BAM conversion: {e}')

def bam_to_sorted_bam(output_path):
    """
    Converts form BAM to Sorted BAM
    """
    # some path string manipulations to get new path for sorted BAM file
    bam_path = output_path
    dir_name = os.path.dirname(bam_path)
    base_name = os.path.basename(bam_path)
    sorted_bam_name = "sorted." + base_name 
    sorted_bam_path = os.path.join(dir_name, sorted_bam_name)
    cmd = [
        "samtools",
        "sort",
        output_path,
        '-o',
        sorted_bam_path
    ]

    try:
        subprocess.run(cmd, check=True)
        print('Successfull BAM to Sorted BAM')
    except subprocess.CalledProcessError as e:
        print(f'Unsuccessfull BAM to Sorted BAM: {e}')

def sorted_bam_to_sorted_bam_index(output_path):
    """"
    Converts Sorted BAM to Sorted Indexed BAM
    """

    # file path string manipulation
    dir_name = os.path.dirname(output_path)
    base_name = os.path.basename(output_path)
    sorted_bam_file_name = 'sorted.' + base_name
    sorted_bam_path = os.path.join(dir_name, sorted_bam_file_name)

    cmd = [
        'samtools',
        'index',
        sorted_bam_path
    ]

    try: 
        subprocess.run(cmd, check=True)
        print('Sucessfully Created BAM Index File')
    except subprocess.CalledProcessError as e:
        print(f'Unsuccessful Creation of BAM Index File: {e}')


"""
Main
____
"""

def main():
    """
    After Running, Should End Up With BAM, Sorted BAM, and BAM Index for the Specified SAM File
    """
    sam_to_bam(sam_path, output_path)
    bam_to_sorted_bam(output_path)
    sorted_bam_to_sorted_bam_index(output_path)

main()






