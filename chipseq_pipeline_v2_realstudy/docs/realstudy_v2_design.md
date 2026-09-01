# Realstudy v2 Scientific Design

## Question and first principles

For two real *C. elegans* ChIP-seq experiments, how many eligible mapped
fragments from the experiment's matched input control are needed before MACS2
peak calls stop changing materially under one fixed workflow?

The treatment data define observed enrichment plus experimental biases. The
matched input captures observable background from that context. Realstudy v2
keeps one pooled treatment BAM fixed and only removes fragments from the pooled
empirical control. Every smaller control is an exact subset of every larger
control for the same seed. This isolates control sampling depth while retaining
the empirical fragment distribution in the prepared library.

This estimates what would have happened if fewer already-observed fragments had
been sequenced from the same library. It assumes that the parent's empirical
fragment distribution is the relevant fixed population. It does not model a
new preparation, batch, animal, antibody, or library chemistry.

## Locked studies and matrix

- Narrow: modERN ceh-27, ENCSR951SLK with ENCSR925LDI matched input;
  39,587,306 treatment and 58,638,263 control raw reads.
- Broad: GSE67028 adult `fer-1` H3K9me2 with SRR1917671 input;
  82,391,358 treatment and 64,271,618 control raw reads.
- SKN-1 remains in the audit but is excluded because 2,963,282 raw control
  reads cannot support the mapped-fragment range.

All FASTQs use the same ce11 Bowtie2/SAMtools path. Replicates are aligned and
QC'd separately, then pooled. With genome size 100,286,401 bp and 150-bp
normalization fragments, the locked targets are:

| Depth | Exact fragments |
|---:|---:|
| 0.5× | 334,288 |
| 1× | 668,576 |
| 2× | 1,337,152 |
| 4× | 2,674,304 |
| 8× | 5,348,608 |
| 16× | 10,697,216 |
| 32× | 21,394,432 |

Seeds 11, 23, and 37 produce 42 sampled calls; two full-control calls are
anchors. A control parent below 21,394,432 eligible fragments fails explicitly.

## Sampling and analysis boundary

Unmapped, secondary, supplementary, and QC-failed alignments are excluded;
duplicate-marked fragments remain and are measured. Fragment identity is read
group plus query name, so paired mates share one decision. SHA-256 ranks include
algorithm version, parent checksum, seed, read group, and query name. A
disk-backed SQLite ledger selects exact nested identifiers and all seven BAMs
are written from one coordinate-ordered pass.

The full-control peak set is an anchor, not truth. Primary comparisons use any
overlap; 10% and 50% reciprocal overlap are sensitivity checks. Metrics include
anchor retention, query concordance, base Jaccard, count/base ratios, score-rank
correlation, narrow summits, broad boundaries/widths, adjacent gain, seed
reproducibility, and ce11 genomic-context strata. Conclusions are limited to
these two experiments and this fixed workflow.
