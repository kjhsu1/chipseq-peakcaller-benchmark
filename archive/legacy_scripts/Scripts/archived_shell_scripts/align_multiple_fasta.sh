#!/bin/zsh

# change path to directory containing all the chip-seq reads FASTA that you want to align
# all reads must be coming from the same genome
for filepath in ../Generated_Data/Reads_FASTA/Different_Coverages_Unbiased/*
do
	filename=$(basename "$filepath")
	echo "Aligning ${filename}"

	# change to appropriate genome index directory
	bowtie2 -x ../Bowtie2_Stuff/index/random_genome_1 \
	-f -U ${filepath} \
	-S ../Generated_Data/Alignments/${filename}.sam
done