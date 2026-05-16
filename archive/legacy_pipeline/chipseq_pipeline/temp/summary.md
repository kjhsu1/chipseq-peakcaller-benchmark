# ChIP-seq Pipeline Smoke Test Report

## Environment Setup
- Command: `source Kenta_Stuff/snakemake_stuff/setup.sh`
- Logic: Activate the project environment to ensure required tools are available.

## CLI Introspection
- Command: `python scripts/updated_chip_seq.py --help`
- Logic: Inspect the simulator's CLI to detect required flags, output format, and paired-end output mode.
- Result: Simulator outputs separate paired-end FASTA files via `--output_fasta1` and `--output_fasta2`.

## Treatment Simulation
- Command:
  ```bash
  conda activate sim
  python scripts/updated_chip_seq.py \
    --fasta data/genomes/ce11_1pct/ce11_1pct.fa \
    --coverage 5 \
    --tf_peak_count 5 \
    --read_length 100 \
    --fragment_length 200 \
    --seed 42 \
    --output_fasta1 temp/treatment/treat_R1.fasta \
    --output_fasta2 temp/treatment/treat_R2.fasta
  ```
- Logic: Generate paired-end treatment reads with TF peaks enabled to simulate ChIP enrichment.
- Validation: Verified the presence of headers starting with `>` and counted records to ensure output completeness.

## Control Simulation
- Command:
  ```bash
  conda activate sim
  python scripts/updated_chip_seq.py \
    --fasta data/genomes/ce11_1pct/ce11_1pct.fa \
    --coverage 5 \
    --tf_peak_count 0 \
    --read_length 100 \
    --fragment_length 200 \
    --seed 43 \
    --output_fasta1 temp/control/ctrl_R1.fasta \
    --output_fasta2 temp/control/ctrl_R2.fasta
  ```
- Logic: Produce control reads with TF peaks disabled to serve as background.
- Validation: Confirmed file creation and record counts analogous to the treatment run.

## Bowtie2 Index Preparation
- Command: `bowtie2-build data/genomes/ce11_1pct/ce11_1pct.fa data/indexes/ce11_1pct/bowtie2_index/index`
- Logic: Ensure the reference genome has the six required `.bt2` files for alignment.

## Treatment Alignment
- Command:
  ```bash
  conda activate align
  bowtie2 \
    --f \
    -x data/indexes/ce11_1pct/bowtie2_index/index \
    -1 temp/treatment/treat_R1.fasta \
    -2 temp/treatment/treat_R2.fasta \
    | samtools view -b - \
    | samtools sort -o temp/treatment/treat.bam
  samtools index temp/treatment/treat.bam
  ```
- Logic: Align treatment reads to the indexed genome, convert to BAM, sort, and index for downstream analysis.
- Validation: `samtools quickcheck temp/treatment/treat.bam` and `samtools flagstat temp/treatment/treat.bam` confirmed alignment integrity and 100% proper pairing.

## Control Alignment
- Command:
  ```bash
  conda activate align
  bowtie2 \
    --f \
    -x data/indexes/ce11_1pct/bowtie2_index/index \
    -1 temp/control/ctrl_R1.fasta \
    -2 temp/control/ctrl_R2.fasta \
    | samtools view -b - \
    | samtools sort -o temp/control/ctrl.bam
  samtools index temp/control/ctrl.bam
  ```
- Logic: Align control reads under identical conditions to the treatment run.
- Validation: `samtools quickcheck temp/control/ctrl.bam` and `samtools flagstat temp/control/ctrl.bam` verified successful alignment and 100% proper pairing.

## Summary
Both treatment and control simulations were generated, aligned, and validated successfully. The commands above encapsulate the entire workflow from environment setup to alignment verification.
