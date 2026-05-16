# Revised 11-Step Benchmark Plan

## First Principles

The controlled pipeline remains the truth-aware benchmark. It is intentionally
simplified so we can isolate how specific biases affect precision, recall, and
F1. The realistic pipeline should test whether those lessons survive a more
realistic experimental generator, but it should not pretend that one observed
real pileup is the true read-generation distribution.

The revised plan keeps real-study grounding and ontology windowing, but delegates
realistic read generation to ChIPs where possible.

## Steps

1. Stabilize the controlled baseline and treat old six TF-clean plots as
   comparison artifacts until rerun/audited under repaired logic.
2. Document controlled-simulator nomenclature, known limitations, and fixed
   metric definitions.
3. Use BAM-only realism scorecards as descriptive audits, not as truth claims.
4. Keep selected real-study manifests and downloaded worm broad-mark inputs as
   reference assets.
5. Learn ChIPs model parameters from selected real-study aligned BAMs and
   peak/truth-proxy regions.
6. Generate ChIPs treatment and WCE control FASTQs across coverage-depth sweeps.
7. Align ChIPs FASTQs with the existing Snakemake alignment structure.
8. Call peaks with MACS2 and, in a matched controlled comparison, EPIC2.
9. Build ontology windows from reference/truth-proxy regions, background
   properties, and simulator truth outputs when available.
10. Evaluate precision, recall, F1, false positives, and control response by
    ontology class.
11. Assemble final outputs that separate controlled truth-aware conclusions,
    ChIPs-based realism conclusions, unresolved limitations, and next steps.

## Ontology Handoff

The ontology remains the interpretive endpoint. Windows should be labeled using
controlled truth intervals in `chipseq_pipeline_v2`, and using ChIPs peak inputs,
learned/reference regions, and real-study annotations as truth-proxy evidence in
`chipseq_pipeline_v2_realstudy`. Real-study regions must be described as
reference or truth-proxy regions, not direct biological truth.

## Current Execution State

- Controlled bug fixes are prerequisite context for any new interpretation.
- Old six TF-clean plots remain baseline comparison artifacts until rerun or
  audited under the repaired logic.
- The first realistic workflow target is ChIPs model learning and paired FASTQ
  generation, then reuse of existing alignment and peak-calling patterns.
- EPIC2 is prepared as a matched controlled-pipeline follow-up sweep to test
  whether conclusions are MACS2-specific.

## What Changed From The Naive Plan

- Real pileups are no longer converted directly into "true" PMFs.
- ChIPs becomes the realistic read-generation engine because it models
  shearing, pulldown, PCR, sequencing, and can learn parameters from real data.
- Ontology evaluation remains central, but its truth inputs are explicitly
  labeled as controlled truth, simulator truth, or real-study truth-proxy.
- Cluster-only work is documented through Slurm scripts and local dry-runs rather
  than pretending local Mac execution can validate EPIC2.
