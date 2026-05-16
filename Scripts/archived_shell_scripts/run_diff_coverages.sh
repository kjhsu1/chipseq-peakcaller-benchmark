#!/bin/zsh

for coverage in 1 5 10 25 50 100
do
	echo='Running chip_seq.py at coverage: $coverage'
	python3 chip_seq.py \
	--fasta ../Genomes/random_genome_1.fa.gz \
	--coverage $coverage \
	--num_bg_peaks 0 \
	--num_fg_peaks 0 \
	--fragment_length 20 \
	> ../Generated_Data/random_genome_1_cov_${coverage}_reads.fa
done


