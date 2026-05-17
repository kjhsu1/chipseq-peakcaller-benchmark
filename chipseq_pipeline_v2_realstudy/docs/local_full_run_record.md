# Local Full Realstudy ChIPs Run Record

## Objective

Prepare real-study ChIP-seq inputs for ChIPs, validate the local Snakemake path
end to end, estimate local scaling, and decide whether a roughly 500-run local
sweep should be launched.

## Boundary

All edits and generated records for this goal stay inside
`chipseq_pipeline_v2_realstudy/`.

## Required Real Inputs

For `chips learn`, each study needs:

- sorted and indexed treatment BAM
- matching treatment/control-derived peak BED
- score column for the BED file

Raw FASTQ studies must first be aligned. Processed BAM studies can be used
directly after download, sorting/index checking, and merge steps.

## Study Inputs

- Narrow study: `encode_skn1_celegans_l1`
  - treatment BAM: `ENCFF775SIN`
  - treatment BAM: `ENCFF927CRC`
  - control BAM: `ENCFF904DXM`
- Broad study: `geo_gse67028_celegans_h3k9me2_adult`
  - treatment FASTQ: `SRR1917669`
  - treatment FASTQ: `SRR1917670`
  - input FASTQ: `SRR1917671`

## Reference Inputs

- `data/local/references/ce11/genome.fa`
- `data/local/references/ce11/genome.fa.fai`
- `data/local/indexes/ce11/bowtie2/ce11*.bt2`

## Current Log

- Started local preparation goal.
- Confirmed only tiny smoke assets were present at goal start.
- Found ENCODE SKN-1 manifest URLs pointed at file pages rather than direct BAM
  download URLs; patching to direct `@@download/*.bam` URLs before download.
- 2026-05-16 22:39 PDT: Preserving existing repo-level `/goals` history. For
  this goal, continuing the durable run record here inside the realstudy
  directory to respect the active write boundary.
- 2026-05-16 22:39 PDT: Regenerated `metadata/data_manifest.csv` from
  `manifests/study_selection.yaml` and the patched
  `manifests/study_file_manifest.csv`.
- 2026-05-16 22:39 PDT: `python -m py_compile scripts/*.py` passed after the
  manifest URL patch.

## Validation Log

- `python scripts/fetch_real_study_data.py --study-selection manifests/study_selection.yaml --file-manifest manifests/study_file_manifest.csv --output-dir analysis_outputs/realstudy_ingest_prep`
  completed successfully.
- `python -m py_compile scripts/*.py` completed successfully.
- First focused download dry-run found a Snakemake wildcard ambiguity in the
  old marker pattern `download_{study_id}_{role}.done`; study IDs and roles both
  contain underscores.
- Patched download markers to
  `analysis_outputs/realstudy_ingest_prep/download_markers/{study_id}/{role}.done`
  and updated all workflow references inside the realstudy pipeline.
- Started six-file download workflow with three concurrent jobs.
- Broad input/control `SRR1917671.fastq.gz` completed: `3.0G`.
- Broad treatment rep 1 `SRR1917669.fastq.gz` completed: `2.9G`.
- ENCODE SKN-1 download failed on the script's `HEAD` preflight with HTTP `403`
  from the redirected S3 URL.
- Patched `download_real_study_file.py` so failed/unsupported `HEAD` requests
  fall back to streamed `GET`; this keeps ENCODE downloads usable while still
  using remote size/range information when available.
- Resumed the missing downloads after the downloader patch.
- Narrow SKN-1 treatment BAM `ENCFF775SIN.bam` completed: `120M`.
- Narrow SKN-1 treatment BAM `ENCFF927CRC.bam` completed: `174M`.
- Narrow SKN-1 control BAM `ENCFF904DXM.bam` completed: `137M`.
- Broad treatment rep 2 `SRR1917670.fastq.gz` completed: `864M`.
- All six selected real-study files are now staged under `data/raw/`.
- Attempted to stage full `ce11` reference FASTA and Bowtie2 index under
  `data/local/`, but the required network/escalated command was blocked by the
  current approval/usage limit. No workaround attempted.
- Decision: continue with SKN-1 narrow BAM preprocessing and ChIPs model
  learning, because that path needs BAM plus peak BED and does not require the
  reference FASTA/index until simulated reads are generated/aligned.
- SKN-1 `chips learn` dry-run was clean: two merge jobs, one MACS2 real-study
  peak-call job, and one ChIPs learn job.
- First SKN-1 execution failed during indexing of the merged ENCODE control BAM:
  `Chromosome blocks not continuous`. This confirms the downloaded ENCODE BAMs
  cannot be assumed coordinate-sorted/index-ready after simple merge.
- Patched `merge_selected_realstudy_bams` so merged real-study BAMs are always
  coordinate-sorted before indexing.
- SKN-1 merge and MACS2 peak calling then completed successfully.
- First ChIPs learn attempt failed with `std::invalid_argument: stof: no
  conversion` because the workflow passed `-c 4`. For MACS2 narrowPeak and
  broadPeak output, column 4 is peak name and column 5 is the numeric score.
- Patched ChIPs `peak_score_column` from `4` to `5` for the MACS2 peak BED
  files used by this workflow.
- Re-ran SKN-1 ChIPs learn after the score-column fix.
- SKN-1 ChIPs learn completed successfully.
- SKN-1 generated model:
  `analysis_outputs/chips_models/encode_skn1_celegans_l1/encode_skn1_celegans_l1.json`.
- SKN-1 real-study peak BED:
  `analysis_outputs/realstudy_peakcalls/encode_skn1_celegans_l1/encode_skn1_celegans_l1_peaks.bed`,
  `3632` peaks.
- SKN-1 merged treatment BAM: `250M`.
- SKN-1 merged control BAM: `123M`.
- Updated `configs/chips_local_full.yaml` so local full ChIPs runs include both
  `encode_skn1_celegans_l1` and `geo_gse67028_celegans_h3k9me2_adult`.
- Current full prototype matrix across the two studies resolves to `144` runs:
  `72` narrow SKN-1 runs and `72` broad worm H3K9me2 runs.
- `pytest tests` passed after the download/preprocessing fixes: `17 passed`.

## Current Blocker

Full local end-to-end simulation/alignment still needs a full `ce11` FASTA and
Bowtie2 index staged under:

- `data/local/references/ce11/genome.fa`
- `data/local/references/ce11/genome.fa.fai`
- `data/local/indexes/ce11/bowtie2/ce11*.bt2`

The attempted network/escalated reference download was blocked by the current
approval/usage limit. Resume from this point by staging those reference assets
inside `chipseq_pipeline_v2_realstudy/data/local/`.

## 2026-05-17 Resume Validation

- Resumed the local full-run goal on branch
  `codex/chips-realstudy-bugfix-wave`.
- Staged full `ce11` reference assets inside the realstudy directory:
  - `data/local/references/ce11/genome.fa`
  - `data/local/references/ce11/genome.fa.fai`
  - `data/local/indexes/ce11/bowtie2/ce11*.bt2`
- Full local dry-run now passes with the staged reference/index:
  `./run_chips_realsim_local.sh`.
- Dry-run DAG for the two-study ChIPs matrix:
  - `728` total jobs
  - `144` final simulated MACS2 peak-call runs
  - `288` simulated-read Bowtie2 alignments
  - `288` ChIPs simulated-read jobs
  - `3` real FASTQ alignments for the broad study
  - `1` real-study MACS2 call and `1` broad-study ChIPs learn job
- Broad worm H3K9me2 preprocessing completed locally:
  - treatment rep 1: `62,865,092` reads, `62.94%` alignment
  - treatment rep 2: `19,526,266` reads, `82.65%` alignment
  - input/control: `64,271,618` reads, `67.65%` alignment
  - broad real-study MACS2 peaks: `10,358`
  - broad ChIPs model:
    `analysis_outputs/chips_models/geo_gse67028_celegans_h3k9me2_adult/geo_gse67028_celegans_h3k9me2_adult.json`
- First broad model attempt failed because Snakemake was invoked with an
  absolute binary path but without the `background_project` environment on
  `PATH`, so `macs2` was not found. Re-running with
  `PATH=/opt/anaconda3/envs/background_project/bin:$PATH` fixed the issue.
- Tiny local end-to-end validation passed for one narrow and one broad run:
  - command runtime: `real 557.96s`
  - narrow target:
    `results_chips/encode_skn1_celegans_l1_t5_c0p5_s11_fl150_rl38_albowtie2_pcmacs2_mnarrow/peaks/macs2/encode_skn1_celegans_l1_t5_c0p5_s11_fl150_rl38_albowtie2_pcmacs2_mnarrow_peaks.bed`
  - broad target:
    `results_chips/geo_gse67028_celegans_h3k9me2_adult_t5_c0p5_s11_fl150_rl38_albowtie2_pcmacs2_mbroad/peaks/macs2/geo_gse67028_celegans_h3k9me2_adult_t5_c0p5_s11_fl150_rl38_albowtie2_pcmacs2_mbroad_peaks.bed`
  - narrow simulated peaks: `3266`
  - broad simulated peaks: `4407`
- Larger local timing validation passed for four additional final runs:
  - command runtime: `real 1468.12s`
  - narrow `t5_c1_s11` simulated peaks: `3430`
  - broad `t5_c1_s11` simulated peaks: `5864`
  - narrow `t10_c0p5_s11` simulated peaks: `3275`
  - broad `t10_c0p5_s11` simulated peaks: `4411`
- Found and fixed a run-ID consistency bug: the standalone params-table script
  parsed integer depths through `argparse` as floats and emitted IDs like
  `t5p0`, while the live Snakefile config emitted runnable IDs like `t5`.
  Added `format_run_id_value()` so integer-valued floats normalize to integer
  labels while true decimals still use `p`, for example `0.5 -> 0p5`.
- Regenerated `metadata/prototype_run_table.csv` and
  `metadata/prototype_run_table.summary.json` through Snakemake after the
  formatter fix. Verification:
  - `144` rows
  - `72` rows per selected study
  - no run IDs containing `_t5p0_`
  - run IDs containing `_t5_` are present
- Validation after the formatter fix:
  - `python -m py_compile scripts/*.py` passed
  - `pytest tests` passed: `18 passed`
- Final full local ChIPs dry-run after validations and the run-ID fix passed.
  Because six of the `144` final runs have already completed locally, the
  remaining dry-run DAG now resolves to `691` jobs:
  - `138` remaining final MACS2 peak-call runs
  - `276` remaining Bowtie2 alignments
  - `138` remaining ChIPs treatment simulations
  - `138` remaining ChIPs control simulations

## Sweep Decision

Do not start the approximately 500-run local sweep.

Evidence:

- Tiny validation lower-bound average: `557.96s / 2 = 278.98s` per final run.
- Larger validation average: `1468.12s / 4 = 367.03s` per final run.
- A 500-run sweep at the optimistic tiny-run average would take about
  `38.7` hours.
- A 500-run sweep at the larger-validation average would take about `51.0`
  hours.
- The full configured local matrix includes higher `20x` treatment depth and
  control depths up to `24x`, so those estimates are likely optimistic rather
  than conservative.
- The current configured two-study matrix has `144` final ChIPs peak-call
  runs. Even using the optimistic tiny-run average, `144` runs estimates to
  about `11.2` hours, above the requested 8-hour local threshold.

Conclusion: the local real-study ChIPs pipeline is now executable and validated
end-to-end on small and larger local batches, but the requested roughly
500-run local sweep should not be launched on this Mac under the 8-hour
constraint. The larger sweep should run on the cluster or be reduced further
before local execution.
