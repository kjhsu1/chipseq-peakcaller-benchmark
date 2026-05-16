# imports
import os

"""
VARS FOR SWEEP (that all aligner x pc combo uses)
"""

genomes = ['600bp', '1pct']
GENOME_PATHS = {
	'600bp': 'genomes/600bp.fa',
	'1pct': 'genomes/1pct.fa'
}
covs = [1, 2, 3, 4, 5]
# pvals = [0.05, 0.01, 0.001]
peak_nums = [1, 2, 3]
peak_tallness = [10, 100, 1000]
peak_broadness = [10, 100, 1000]
#aligner = ['bowtie2', 'bwa-mem']
#peak_callers = ['macs2', 'epic2']

"""
VARS FOR SWEEP (specific to aligner x pc combo)
"""

# stores both genome name, genome.fa, and genome index file in a tuple; for bowtie+bwa mem
bowtie_idx = {
	'600bp': '/path/to/600bp_genome_index', 
	'1pct': '/path/to/1pct_c_elegans_index'
	}
bwa_idx = {
	'600bp': '/path/to/600bp_genome_index', 
	'1pct': '/path/to/1pct_c_elegans_index'
	}
nomodel = [True, False] # only for macs2


"""
CONSTANT VARS
"""

RESULTS = '/results/'

# chipseq reads
READS = os.path.join(RESULTS, 'reads/')

# reads aligned
ALIGN = os.path.join(RESULTS, 'align/')

# peaks
PEAKS = os.path.join(RESULTS, 'peaks/')


"""
Rules
"""

# note: have one rule each for bowtie2 and bwa-mem alignment
	# no matter how complex the DAG is, snakemake should be able to handle it
	# as long as the output file paths are unique
	
rule all:
	# final name includes all vars for sweep
	# concern: we would have to include parameters specific to specifc aligner x pc combos
		# solution: could just create separate expand() for each specific aligner x pc combo

	input:
		# for bowtie2 x macs2
			# {nomodel} is unique to macs2
		expand('bowtie2_macs2'
				'_g_{g}'
				'_cov_{cov}'
				'_pnum_{pnum}'
				'_ptall_{ptall}'
				'_pbroad_{pbroad}'
				'_nomodel_{nomodel}'
				'.peaks',
				g=genomes, 
				cov=covs, 
				pnum=peak_nums, 
				ptall=peak_tallness,
				pbroad=peak_broadness,
				nomodel=nomodel
			),

		# bowtie2 x epic2
			# no unique vars
		expand('bowtie2_epic2'
				'_g_{g}'
				'_cov_{cov}'
				'_pnum_{pnum}'
				'_ptall_{ptall}'
				'_pbroad_{pbroad}'
				'.peaks',
				g=genomes, 
				cov=covs, 
				pnum=peak_nums, 
				ptall=peak_tallness,
				pbroad=peak_broadness
			),

		# bwa x macs2
			# {nomodel} is unique to macs2
		expand('bwa_macs2'
				'_g_{g}'
				'_cov_{cov}'
				'_pnum_{pnum}'
				'_ptall_{ptall}'
				'_pbroad_{pbroad}'
				'_nomodel_{nomodel}'
				'.peaks',
				g=genomes, 
				cov=covs, 
				pnum=peak_nums, 
				ptall=peak_tallness,
				pbroad=peak_broadness,
				nomodel=nomodel
			),

		# bwa x epic2
			# no unique vars
		expand('bwa_epic2'
			'_g_{g}'
			'_cov_{cov}'
			'_pnum_{pnum}'
			'_ptall_{ptall}'
			'_pbroad_{pbroad}'
			'.peaks',
			g=genomes, 
			cov=covs, 
			pnum=peak_nums, 
			ptall=peak_tallness,
			pbroad=peak_broadness
		)
	
rule chipseq:
	params:
		'''
		-essentially I wanted a wildcard genome path var that is synced with {g}, 
		the genome name wildcard
		- below does that by creating a function that takes in w, the wildcard object
		and then extracting the genome name of the particular run, w.g (ex. w.g=600bp)
		- then it gets the genome path of the genome name from a dictionary.

		'''
		g_path = lambda w: GENOME_PATHS[w.g] 
	output: 
		exp_reads = (
			READS
			+ 	'/g_{g}'
			+	'_cov_{cov}'
			+	'_pnum_{pnum}'
			+	'_ptall_{ptall}'
			+	'_pbroad_{pbroad}'
			+	'exp_chipseq_reads.fa'
			),

		con_reads = (
			READS
			+ 	'/g_{g}'
			+	'_cov_{cov}'
			+	'_pnum_{pnum}'
			+	'_ptall_{ptall}'
			+	'_pbroad_{pbroad}'
			+	'control_chipseq_reads.fa'
			)

	shell:
		"""
		python3 /scripts/chip_seq.py \
			--genome {params.g_path} \
			--cov {cov} \
			--pnum {pnum} \
			--ptall {ptall} \
			--pbroad {pbroad} \
			> {output.exp_reads}

		python3 /scripts/chip_seq.py \
			--genome {params.g_path} \
			--cov {cov} \
			--pnum {pnum} \
			--ptall {ptall} \
			--pbroad {pbroad} \
			> {output.control_reads}
		"""

rule align_bowtie2:
	input:
		exp_reads = chipseq.output.exp_reads
		con_reads = chipseq.output.con_reads
	params:
		# params are evaluated right before running the shellc commands
		g_path = lambda w: GENOME_PATHS[w.g]
		g_index = lambda w: bowtie_idx[w.g]
	output:
		# NOTE: most peakcallers will not use default SAM or BAM, so exclude
		
		# sorted bam
		sorted_bams=(
			ALIGN
			+	'{exp_con}_' # Q. is there a way to store these universal parameter sweep variables
			+	'g_{g}_'				# in a var somewhere so I don't have to write it out everytime?
			+	'cov_{cov}_'
			+	'pnum_{pnum}_'
			+	'ptall_{ptall}_'
			+	'pbroad_{pbroad}'
			+	'.sorted.bam'
			)

		# sorted bai
		sorted_bais=(
			ALIGN
			+	'{exp_con}_' # Q. is there a way to store these universal parameter sweep variables
			+	'g_{g}_'				# in a var somewhere so I don't have to write it out everytime?
			+	'cov_{cov}_'
			+	'pnum_{pnum}_'
			+	'ptall_{ptall}_'
			+	'pbroad_{pbroad}'
			+	'.sorted.bai'
			)

	shell:
		# first need to create genome index for the two genomes in {g}
		'samtools createidx {params.g_path} > {params.g_index}' 

		# For EXP
		# then align using created genome index
		bowtie2 align -genome_index (ALIGN, '/{g}_index') -reads_fasta align_bowtie2.input[1] | \
		bowtie2 ... | \ # stdout (the SAM) pipe into BAM 
		bowtie2 ... > # convert to sorted.bam save the file in ALIGNMENT
		
		# Q. How does snakemake treat when the grid size differs in one command
			# ex. bowtie2 align -genome_index {os.path.join(ALIGN, '{gen}_index')} -reads_fasta align.input[1] | \
			# in this example, we are creating the aligning and creating SAM file
			# however, there are two sample space grid with different dimensions
			# for the genome_index argument the grid contains just two things; {gen} = 1pct_index or 600bp_index
			# however, the align_bowtie2.input[1] grid space is all combos of nonspecific parameter sweep vars
				# ie. peak_num x peak_tallness x genome x peak_broadness
			# how does snakemake run these differing grid spaces?

rule align_bwa:








	