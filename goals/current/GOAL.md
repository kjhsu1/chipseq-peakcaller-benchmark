# Goal: Realstudy v2 Empirical Control-Depth Study

Deliver a code-complete, audit-ready workflow that uses two deeply sequenced
real matched ChIP/input experiments to determine the lowest tested empirical
control depth at which MACS2 peak calls stabilize, without using ChIPs or making
biological-truth claims.

## Context

- Project: Realstudy v2 control-depth sufficiency paper workflow.
- Stack: Python, pandas, NumPy, SciPy, pysam, SQLite, matplotlib, Snakemake,
  Bowtie2, SAMtools, MACS2, Bash/Slurm.
- Branch/worktree: `codex/realstudy-v2-code` in
  `.codex_worktrees/realstudy-v2-code`.
- Inputs: modERN ceh-27/ENCODE matched input and GSE67028 fer-1 H3K9me2/input.
- Constraint for this wave: author and statically check code only. Do not
  download real FASTQs; do not run tests, Snakemake, alignment, sampling,
  MACS2, figures, scientific analysis, or HPC submission.

## Success criteria

1. Exactly 42 sampled-control calls plus two full-control anchors are encoded.
2. All seven FASTQ inputs and paper evidence sources are versioned with counts,
   checksums, URLs, provenance, and inclusion rationale.
3. Deterministic exact nested fragment sampling is disk-backed, pair-safe, and
   fails below 21,394,432 eligible control fragments.
4. Narrow/broad MACS2 calls hold treatment fixed and use `--keep-dup all`.
5. Anchor-relative metrics, reciprocal-overlap sensitivity, seed uncertainty,
   genomic strata, and the registered sufficiency decision are implemented.
6. A transactional 23-table SQLite database and checksummed CSV exports are
   implemented with constraints and foreign keys.
7. Six publication figures plus supplementary-source contracts are implemented
   in PDF/SVG/300-DPI PNG with source CSV and caption outputs.
8. Local/Slurm launchers, focused future tests, docs, and final integrity
   validation are present without executing the scientific workflow.

## Claim boundary

The eventual paper may report control-depth effects on MACS2 stability for these
two experiments under this workflow. It may not claim biological precision or
recall, new library/batch simulation, biological replication from seeds,
universal depth sufficiency, or biological truth from the full-control anchor.
