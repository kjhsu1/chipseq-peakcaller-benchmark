# Log

## 2026-09-01 — Realstudy v2 empirical-control implementation

- Assumption: the active finalized plan is the user-pasted empirical
  control-depth plan; it supersedes the older SKN-1/five-depth/ChIPs v2 scaffold.
- Preserved the committed legacy ChIPs pipeline and reverted only overlapping
  uncommitted scaffold edits in legacy files.
- Locked primary inputs to ENCSR951SLK + ENCSR925LDI (ceh-27) and GSE67028
  (fer-1 H3K9me2). SKN-1 is retained only as an excluded audit candidate.
- Corrected the Guo paper registry against PubMed: title is “Enrichment of
  H3K9me2 on Unsynapsed Chromatin in Caenorhabditis elegans Does Not Target de
  Novo Sites” and PMCID is PMC4555223.
- Important sampling assumption: downsampling estimates reduced sequencing from
  the same empirical parent library; the fragment distribution is held fixed.
  It does not create new experimental preparations or biological replicates.
- Important eligibility decision: duplicate-marked fragments remain eligible
  and are measured; unmapped, secondary, supplementary, and QC-failed records
  are excluded. Read-group + query-name identity keeps paired mates together.
- Important failure behavior: parent controls below 21,394,432 eligible mapped
  fragments fail; no depth is silently removed or relabeled.
- Important analysis decision: full-control calls are anchors, never truth.
  Precision/recall/F1 truth language was removed from the v2 path.
- Important database decision: runtime tables are consolidated before a
  temporary SQLite transaction; publication occurs only after foreign-key
  validation. Every table exports with a row count and SHA-256.
- Important figure decision: representative loci are computed from actual BAMs
  and peak artifacts; no placeholder panel is generated.
- Static validation performed: all new Python scripts/tests passed
  `python -m py_compile`; all new launchers passed `bash -n`; `git diff --check`
  passed.
- Scope boundary honored: no FASTQ download, test execution, Snakemake dry-run,
  scientific tool, figure render, local analysis, or HPC job was executed.
