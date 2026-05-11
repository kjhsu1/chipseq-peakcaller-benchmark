# Real SKN-1 ENCODE/modERN rebuild

This directory contains a local rebuild of real `C. elegans` SKN-1 ChIP-seq signal
tracks from official ENCODE ce11 BAMs for experiment `ENCSR012BJM`.

Goal:
- create a pooled treatment-only BigWig for quantitative comparison
- optionally create a pooled control-normalized BigWig from the same rebuild
- keep the steps explicit and reproducible

Tracked provenance in git:
- this `README.md`
- `logs/*.flagstat.txt`
- `meta/input_manifest.tsv`
- `meta/output_manifest.tsv`
- `meta/tool_versions.txt`
- `../../scripts/rebuild_real_encode_skn1_bigwigs.sh`

Intentionally untracked artifacts:
- `downloads/*.bam`
- `bam/*.bam`
- `bam/*.bai`
- `bw/*.bw`

Chosen strategy:
- source BAMs from official ENCODE ce11 alignments
- use the replicate BAMs that ENCODE already labels as filtered alignments
- pool treatment replicates with `samtools merge`
- use the pooled input BAM for control-aware signal
- generate standardized BigWigs locally with deepTools

Files used:
- treatment replicate 1 BAM: `ENCFF775SIN`
- treatment replicate 2 BAM: `ENCFF927CRC`
- pooled input BAM: `ENCFF904DXM`

Signal products to build:
- pooled treatment-only CPM BigWig
- pooled treatment-vs-input log2 BigWig

Method choices:
- genome assembly: `ce11`
- signal bin size: `10 bp`
- treatment-only normalization: `CPM`
- control-aware normalization: `log2` ratio via `bamCompare`
- single-end fragment extension: `125 bp`

Why `125 bp`:
- ENCODE QC for the two treatment BAMs reports fragment lengths of `130` and `120`
- this rebuild uses `125 bp` as a simple pooled approximation

Caveats:
- this is a standardized local rebuild, not an attempt to exactly reproduce the
  original ENCODE BigWig byte-for-byte
- ENCODE portal metadata exposes the file provenance and step classes, but not
  the exact original signal-generation command line

Execution notes:
- the downloaded ENCODE BAMs were not coordinate-sorted as downloaded
- local coordinate-sorted copies were created before indexing or coverage
- the rebuild was executed inside a dedicated git worktree on branch
  `codex/real-encode-skn1-bigwig`
- the exact command sequence has been captured in
  `../../scripts/rebuild_real_encode_skn1_bigwigs.sh`

Executed workflow:
```bash
samtools sort -o bam/ENCFF775SIN.sorted.bam downloads/ENCFF775SIN.bam
samtools sort -o bam/ENCFF927CRC.sorted.bam downloads/ENCFF927CRC.bam
samtools sort -o bam/ENCFF904DXM.sorted.bam downloads/ENCFF904DXM.bam

samtools index bam/ENCFF775SIN.sorted.bam
samtools index bam/ENCFF927CRC.sorted.bam
samtools index bam/ENCFF904DXM.sorted.bam

samtools merge -f bam/ENCFF775SIN_ENCFF927CRC.pooled.bam \
  bam/ENCFF775SIN.sorted.bam \
  bam/ENCFF927CRC.sorted.bam
samtools index bam/ENCFF775SIN_ENCFF927CRC.pooled.bam

bamCoverage \
  -b bam/ENCFF775SIN_ENCFF927CRC.pooled.bam \
  -o bw/ENCFF775SIN_ENCFF927CRC.pooled.ce11.CPM.extend125.bin10.bw \
  --binSize 10 \
  --normalizeUsing CPM \
  --extendReads 125 \
  -p 2

bamCompare \
  -b1 bam/ENCFF775SIN_ENCFF927CRC.pooled.bam \
  -b2 bam/ENCFF904DXM.sorted.bam \
  -o bw/ENCFF775SIN_ENCFF927CRC_vs_ENCFF904DXM.ce11.log2ratio.extend125.bin10.bw \
  --binSize 10 \
  --extendReads 125 \
  --operation log2 \
  --pseudocount 1 1 \
  --scaleFactorsMethod readCount \
  -p 2
```

Built outputs:
- `bw/ENCFF775SIN_ENCFF927CRC.pooled.ce11.CPM.extend125.bin10.bw`
- `bw/ENCFF775SIN_ENCFF927CRC_vs_ENCFF904DXM.ce11.log2ratio.extend125.bin10.bw`

Supporting logs:
- `logs/ENCFF775SIN.sorted.flagstat.txt`
- `logs/ENCFF927CRC.sorted.flagstat.txt`
- `logs/ENCFF904DXM.sorted.flagstat.txt`
- `logs/ENCFF775SIN_ENCFF927CRC.pooled.flagstat.txt`
