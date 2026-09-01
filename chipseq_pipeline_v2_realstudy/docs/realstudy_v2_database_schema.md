# Realstudy v2 Database Schema

The future workflow creates
`analysis_outputs/realstudy_v2/database/realstudy_v2.sqlite` transactionally,
with foreign keys enabled, then atomically publishes it after validation.

- Provenance: `schema_metadata`, `publications`, `publication_claims`,
  `studies`, `experiments`, `samples`, `input_files`,
  `experiment_publications`, `parent_libraries`, `software_versions`.
- Design/execution: `sampling_blocks`, `control_subsamples`, `runs`,
  `artifacts`, `validation_events`.
- Analysis: `peaks`, `reference_regions`, `peak_overlaps`,
  `reference_region_recovery`, `run_metrics`, `seed_pair_metrics`,
  `stratified_metrics`, `enough_control_decisions`.

Stable keys, coordinate/count checks, legal statuses/types, unique
study/seed/depth combinations, complete checksums/anchor linkage, and foreign
keys are enforced in SQL. Indexes cover publication, study/depth/seed/run,
interval, overlap, and metric access. Every table is exported to a deterministic
CSV with row count and SHA-256 in `csv_exports_manifest.csv`.

The hierarchy is publication → experiment → sample → input file → pooled parent
→ exact subsample → MACS2 run → peak/overlap/metric → study decision. This lets
the paper trace every aggregate point back to the source FASTQ and rank cutoff.
