# Two-Run Realstudy Demo Report

## Scope

This report summarizes the two-run realstudy demo executed with the demo pairing config:

- output root: `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2_realstudy/chips_runs/chips_realsim_ontology_demo_20260528_181205`
- config: `chipseq_pipeline_v2_realstudy/configs/chips_demo_pairing.yaml`
- run table: `chipseq_pipeline_v2_realstudy/metadata/chips_demo_pairing_run_table.csv`
- Slurm job: `14561728`

This was a representative demo, not the full sweep. It included exactly two runs:

1. `encode_skn1_celegans_l1_t5_c4_s11_fl150_rl38_albowtie2_pcmacs2_mnarrow`
2. `geo_gse67028_celegans_h3k9me2_adult_t5_c4_s11_fl150_rl38_albowtie2_pcmacs2_mbroad`

## Final Status

Scientifically, both demo runs completed and produced alignment outputs, peak calls, ontology region tables, combined ontology summaries, summary plots, and validation artifacts.

Operationally, the top-level run is marked `RUN_FAILED`, but not because the scientific outputs were missing. The validation report shows:

- expected run count: `2`
- completed run count: `2`
- all peak calls present: `PASS`
- ontology scoring complete for all runs: `PASS`
- combined ontology table built: `PASS`
- final summary tables built: `PASS`
- final plots built: `PASS`

The overall validation failed for workflow-management reasons:

- `RUN_COMPLETE` marker was not present when validation ran
- ontology coverage threshold was intentionally not met by a two-run demo
- the reproducibility package expected production-scale files that the demo wrapper did not copy

So this demo is scientifically usable, but not production-validation clean.

## Studies Ingested

### Study 1: `encode_skn1_celegans_l1`

- source: `ENCODE`
- assembly: `ce11`
- signal class: `narrow`
- selected treatment inputs:
  - `ENCFF775SIN.bam`
  - `ENCFF927CRC.bam`
- selected control input:
  - `ENCFF904DXM.bam`
- ingest mode: processed BAM ingest
- needs alignment: `false`

### Study 2: `geo_gse67028_celegans_h3k9me2_adult`

- source: `GEO`
- assembly: `ce11`
- signal class: `broad`
- selected treatment inputs:
  - `SRR1917669.fastq.gz`
  - `SRR1917670.fastq.gz`
- selected control input:
  - `SRR1917671.fastq.gz`
- ingest mode: FASTQ ingest plus alignment
- needs alignment: `true`

## Demo Run Specifications

Shared settings:

- treatment coverage: `5x`
- control coverage: `4x`
- seed: `11`
- fragment length: `150`
- read length: `38`
- aligner: `bowtie2`
- peak caller: `MACS2`
- ontology enabled: `true`
- local raw-input requirement: `true`

Natural pairing applied:

- SKN-1 run used `MACS2 narrow`
- H3K9me2 run used `MACS2 broad`

## Workflow From Beginning To End

### 1. Study selection and ingest metadata

Purpose:
- define which real studies and files are eligible
- record source, assembly, file format, and whether alignment is needed

Main code:
- `chipseq_pipeline_v2_realstudy/metadata/data_manifest.csv`
- `chipseq_pipeline_v2_realstudy/rules/ingest_real_data.smk`
- `chipseq_pipeline_v2_realstudy/scripts/download_real_study_file.py`

Main file types:
- `.csv` manifest
- raw `.bam`
- raw `.fastq.gz`

### 2. Realstudy alignment and merged realstudy BAM generation

Purpose:
- align FASTQ-based studies
- merge selected treatment and control BAMs per study

Main code:
- `chipseq_pipeline_v2_realstudy/rules/peakcalling.smk`

Main file types:
- aligned `.bam`
- indexed `.bam.bai`

Main structure:
- `data/aligned/{study_id}/{role}/aligned.sorted.bam`
- `analysis_outputs/realstudy_bams/{study_id}/{group}.merged.bam`

### 3. Realstudy truth peak generation

Purpose:
- call MACS2 peaks on the ingested real studies
- use those peaks as the truth-template BED for downstream scoring

Main code:
- `chipseq_pipeline_v2_realstudy/rules/peakcalling.smk`

Main file types:
- `.bed`
- `.narrowPeak` or `.broadPeak`
- `.xls`
- summit files for narrow calls

Main structure:
- `analysis_outputs/realstudy_peakcalls/{study_id}/{study_id}_peaks.bed`

Important detail:
- the ontology analysis does not score genome-wide sliding windows
- it scores one row per truth-template realstudy peak interval

### 4. ChIPs simulation

Purpose:
- learn simulation behavior from real data
- generate simulated treatment and whole-cell-extract control reads

Main code:
- `chipseq_pipeline_v2_realstudy/rules/chips_simulation.smk`

Main file types:
- simulated `.fastq`
- model `.json`

Main structure:
- `analysis_outputs/chips_models/{study_id}/{study_id}.json`
- `results_chips/{run_id}/treat/reads_R1.fastq`
- `results_chips/{run_id}/treat/reads_R2.fastq`
- `results_chips/{run_id}/con/reads_R1.fastq`
- `results_chips/{run_id}/con/reads_R2.fastq`

### 5. Simulated-read alignment

Purpose:
- align the simulated treatment and control reads back to `ce11`

Main code:
- `chipseq_pipeline_v2_realstudy/rules/chips_simulation.smk`

Main file types:
- aligned `.bam`
- `.bam.bai`

Main structure:
- `results_chips/{run_id}/bowtie2/treat/aligned.sorted.bam`
- `results_chips/{run_id}/bowtie2/con/aligned.sorted.bam`

### 6. Simulated peak calling

Purpose:
- call peaks on the simulated treatment/control BAM pair

Main code:
- `chipseq_pipeline_v2_realstudy/rules/chips_simulation.smk`
- `chipseq_pipeline_v2_realstudy/rules/peakcalling.smk`

Main file types:
- `.bed`
- `.narrowPeak` or `.broadPeak`
- `.xls`
- `.bed` summit files for narrow calls

Main structure:
- `results_chips/{run_id}/peaks/macs2/{run_id}_peaks.bed`

### 7. Ontology region scoring

Purpose:
- compare the simulated peak set against the truth-template BED
- compute one scored row per truth-template region

Main code:
- `chipseq_pipeline_v2_realstudy/rules/ontology_analysis.smk`
- `chipseq_pipeline_v2_realstudy/scripts/score_chips_ontology_regions.py`

Main file types:
- `.csv`

Main structure:
- `analysis_outputs/chips_ontology/{run_id}/region_metrics.csv`

Columns include:
- genomic coordinates
- template score
- called overlap
- treatment and control read counts
- control variability
- treatment local peak z-score
- log2 enrichment
- precision, recall, and F1 under the truth-template overlap proxy

### 8. Ontology classification

Purpose:
- map each scored region into ontology labels

Main code:
- `chipseq_pipeline_v2_realstudy/scripts/classify_regions.py`
- `chipseq_pipeline_v2_realstudy/scripts/ontology_lib.py`

Main file types:
- `.csv`
- `.bed`
- `.json`

Main structure:
- `analysis_outputs/chips_ontology/{run_id}/classified.csv`
- `analysis_outputs/chips_ontology/{run_id}/classified_summary.csv`
- `analysis_outputs/chips_ontology/{run_id}/classified.bed`
- `analysis_outputs/chips_ontology/{run_id}/classified_definition.json`

### 9. Cross-run summary evaluation

Purpose:
- combine both runs
- summarize mean precision, recall, and F1 by run, ontology class, and failure mode

Main code:
- `chipseq_pipeline_v2_realstudy/scripts/combine_csv_tables.py`
- `chipseq_pipeline_v2_realstudy/scripts/evaluate_by_region_ontology.py`

Main file types:
- `.csv`
- `.png`
- `.md`
- `.txt`

Main structure:
- `analysis_outputs/chips_ontology/combined_region_metrics.csv`
- `analysis_outputs/chips_ontology/summary/per_run_overall_metrics.csv`
- `analysis_outputs/chips_ontology/summary/per_ontology_metrics.csv`
- `analysis_outputs/chips_ontology/summary/control_response_by_ontology.csv`
- `analysis_outputs/chips_ontology/summary/failure_mode_metrics.csv`
- `analysis_outputs/chips_ontology/summary/ontology_class_coverage.csv`
- `analysis_outputs/chips_ontology/summary/ontology_f1_heatmap.csv`
- `analysis_outputs/chips_ontology/summary/*.png`

## Ontology Definition

The ontology is a Cartesian product of:

- `4` background classes
- `4` foreground classes

Total possible ontology classes: `16`

### Background classes

Definitions come from `chipseq_pipeline_v2_realstudy/scripts/ontology_lib.py`.

1. `hilly_background`
   - criterion: `ctrl_bump_z >= 3.0`
   - interpretation: control contains a strong localized bump

2. `wavy_background`
   - criterion: `ctrl_bump_z < 3.0` and `ctrl_cv >= 0.25`
   - interpretation: control is variable or uneven, but not sharply hotspot-like

3. `flat_background`
   - criterion: `ctrl_bump_z < 3.0` and `ctrl_cv <= 0.25`
   - interpretation: control is comparatively even

4. `mixed_background`
   - criterion: fallback if none of the above match
   - interpretation: intermediate control structure

### Foreground classes

1. `plateau_like_signal`
   - criterion: `plateau_width_bp >= 1000` and `log2_enrichment >= 1.0`

2. `strong_peak`
   - criterion: not plateau-like, and `treat_peak_z >= 6.0`

3. `weak_peak`
   - criterion: not plateau-like, and `3.0 <= treat_peak_z < 6.0`

4. `background_only`
   - criterion: `treat_peak_z < 3.0` and not plateau-like

### Confounders and failure modes

Confounders:

- `hotspot_like_background`
  - assigned when `background_class == hilly_background` and foreground is `weak_peak` or `background_only`

- `plateau_width_conflict`
  - assigned when `plateau_width_bp >= 1000` but the region was not classified as `plateau_like_signal`

- `weak_enrichment`
  - assigned when `log2_enrichment < 1.0` but the region was classified as a non-background foreground

Failure modes:

- `background_dominant`
  - assigned when foreground is `background_only` and background is not `flat_background`

- `confounded_signal`
  - assigned when any confounder is present and the row is not already `background_dominant`

- `none`
  - assigned otherwise

## Ontology Classes Actually Observed

Observed classes across the two-run demo: `7 / 16`

Observed:

- `flat_background__background_only`
- `flat_background__plateau_like_signal`
- `hilly_background__background_only`
- `hilly_background__plateau_like_signal`
- `wavy_background__background_only`
- `wavy_background__plateau_like_signal`
- `wavy_background__weak_peak`

Not observed:

- all `mixed_background__*`
- all `*__strong_peak`
- `flat_background__weak_peak`
- `hilly_background__weak_peak`
- `flat_background__strong_peak`
- `hilly_background__strong_peak`
- `wavy_background__strong_peak`

## Per-Run Results

### Run 1: SKN-1 narrow

Run id:
- `encode_skn1_celegans_l1_t5_c4_s11_fl150_rl38_albowtie2_pcmacs2_mnarrow`

Outputs:

- simulated peak BED lines: `3511`
- ontology region rows: `3632`
- missed truth-template regions: `22`
- per-run precision: `0.9939427313`
- per-run recall: `0.9939427313`
- per-run F1: `0.9939427313`

Ontology distribution:

- `flat_background__background_only`: `1663` regions, `45.79%`
- `wavy_background__background_only`: `1962` regions, `54.02%`
- `wavy_background__plateau_like_signal`: `7` regions, `0.19%`

Main finding:
- this run was almost perfectly recovered under the truth-template overlap proxy
- its missed regions were only `22 / 3632`
- all misses were background-only:
  - `12` `wavy_background__background_only`
  - `10` `flat_background__background_only`

Interpretation:
- the narrow run recovered essentially all obvious truth-template regions
- the few failures were borderline, background-like intervals rather than strong signal classes

### Run 2: H3K9me2 broad

Run id:
- `geo_gse67028_celegans_h3k9me2_adult_t5_c4_s11_fl150_rl38_albowtie2_pcmacs2_mbroad`

Outputs:

- simulated peak BED lines: `6116`
- ontology region rows: `10358`
- missed truth-template regions: `5360`
- per-run precision: `0.4825255841`
- per-run recall: `0.4825255841`
- per-run F1: `0.4825255841`

Ontology distribution:

- `flat_background__background_only`: `1162` regions, `11.22%`
- `flat_background__plateau_like_signal`: `7` regions, `0.07%`
- `hilly_background__background_only`: `3` regions, `0.03%`
- `hilly_background__plateau_like_signal`: `128` regions, `1.24%`
- `wavy_background__background_only`: `6601` regions, `63.73%`
- `wavy_background__plateau_like_signal`: `2443` regions, `23.59%`
- `wavy_background__weak_peak`: `14` regions, `0.14%`

Missed-region ontology breakdown:

- `wavy_background__background_only`: `4259`
- `flat_background__background_only`: `969`
- `wavy_background__plateau_like_signal`: `116`
- `wavy_background__weak_peak`: `11`
- `hilly_background__background_only`: `3`
- `hilly_background__plateau_like_signal`: `1`
- `flat_background__plateau_like_signal`: `1`

Interpretation:
- the broad run is much harder overall than the narrow run under this overlap metric
- most failures came from background-only classes, especially `wavy_background__background_only`
- true plateau-like classes were recovered much better than the background-like classes

## Combined Summary Results

### Per-run performance

- SKN-1 narrow F1: `0.9939`
- H3K9me2 broad F1: `0.4825`

### Per-ontology mean F1

- `hilly_background__plateau_like_signal`: `0.9922`
- `wavy_background__plateau_like_signal`: `0.9527`
- `flat_background__plateau_like_signal`: `0.8571`
- `flat_background__background_only`: `0.6535`
- `wavy_background__background_only`: `0.5012`
- `wavy_background__weak_peak`: `0.2143`
- `hilly_background__background_only`: `0.0`

### Total ontology class counts across both runs

Combined region count: `13990`

- `wavy_background__background_only`: `8563`, `61.21%`
- `flat_background__background_only`: `2825`, `20.19%`
- `wavy_background__plateau_like_signal`: `2450`, `17.51%`
- `hilly_background__plateau_like_signal`: `128`, `0.92%`
- `wavy_background__weak_peak`: `14`, `0.10%`
- `flat_background__plateau_like_signal`: `7`, `0.05%`
- `hilly_background__background_only`: `3`, `0.02%`

### Failure mode summary

- `none`: mean F1 `0.7972`
- `background_dominant`: mean F1 `0.5011`
- `confounded_signal`: mean F1 `0.2143`

## Detailed Findings

1. The ontology space was dominated by background-only classes.
   - `81.4%` of all regions were either `wavy_background__background_only` or `flat_background__background_only`

2. Plateau-like classes were recovered much better than weak or background-only classes.
   - plateau-like classes all had mean F1 above `0.85`
   - `wavy_background__plateau_like_signal` was especially common and still had high recovery at `0.9527`

3. The narrow TF-like run behaved like a near-ideal recovery case.
   - almost perfect F1
   - only `22` misses
   - misses were all borderline background-like intervals

4. The broad histone-like run exposed the main difficulty of this demo setup.
   - overall F1 dropped to `0.4825`
   - the loss was driven mainly by `wavy_background__background_only` regions
   - plateau-like broad regions were much more stable than background-like broad regions

5. This two-run demo is enough to demonstrate pipeline behavior, but not enough to satisfy the production ontology coverage criterion.
   - only `2` runs were summarized
   - the production validation threshold expects at least `3` ontology classes with at least `10` runs each

## Practical Bottom Line

This demo shows that the realstudy pipeline can ingest two real worm studies, simulate matched ChIPs treatment/control reads, align them with Bowtie2, call peaks with MACS2, and classify truth-template regions into ontology categories.

The strongest scientific signal from the demo is:

- narrow SKN-1 is recovered almost perfectly
- broad H3K9me2 is much harder overall
- plateau-like broad regions recover well
- broad, wavy, background-dominated regions are the main failure mode
