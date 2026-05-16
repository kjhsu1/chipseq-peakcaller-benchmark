#!/bin/zsh

# input: directory where you stored all SAM alignment files you want to covert
# output: in same directory, creates BAM, sorted BAM, and BAM index files (.bai)

# change this appropriately
directory="../Generated_Data/Alignments"


for sam in ${directory}/*.sam
do
	base=$(basename ${sam} .sam) # remove path and .sam extension

	# BAM output in the same directory 
	bam=${directory}/${base}.bam

	# SAM to BAM 
	samtools view -Sb ${sam} > ${bam}

	# sorting
	samtools sort ${bam} -o ${directory}/sorted.${base}.bam

	# indexing
	samtools index ${directory}/sorted.${base}.bam
done

