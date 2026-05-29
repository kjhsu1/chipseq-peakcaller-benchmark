# Progress

## Realstudy Production Validation Wave

- [x] New `/goal` objective loaded from `~/Downloads/prompt.txt`.
- [x] Existing `/goals` history preserved instead of deleted.
- [x] `GOAL.md` rewritten around the realstudy production objective.
- [x] `PLAN.md` rewritten around the 11 production success criteria.
- [x] Current realstudy workflow files re-inspected against the new objective.
- [x] Exact production output paths re-confirmed from live workflow code.
- [x] Latest realstudy cluster failure root cause re-confirmed in context of the
  new objective.
- [x] Production validation script added or upgraded.
- [x] Runtime/resource capture implemented in the production path.
- [x] Reproducibility packaging implemented in the production path.
- [x] Required summary plot list finalized.
- [x] Ontology class coverage validation implemented.
- [x] Compile checks rerun under the new validation wave.
- [x] Tests rerun under the new validation wave.
- [x] Production Snakemake dry-run rerun under the new validation wave.
- [x] Realstudy production launch path judged ready and launched on Slurm.
- [x] Live cluster fetch failure diagnosed as external DNS/network access.
- [x] Production ingest switched to pre-staged local BAM/FASTQ inputs.

## Live Production State

- [x] Production Slurm job submitted: `14516460`
- [x] Production output root identified:
  `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2_realstudy/chips_runs/chips_realsim_ontology_128cpu_2tb_20260526_034056`
- [ ] Production run completed successfully.
- [ ] Final validation report generated.
- [ ] Runtime/resource report generated.
- [ ] Reproducibility package generated.

## 11 Success Criteria Status

- [ ] 1. Successful end-to-end production run
- [ ] 2. Complete expected run set produced
- [ ] 3. All simulated alignment outputs present
- [ ] 4. All peak-call outputs present
- [ ] 5. Ontology scoring complete for all runs
- [ ] 6. Combined ontology table built
- [ ] 7. Final performance summary tables built
- [ ] 8. Final performance plots built
- [ ] 9. Ontology class coverage acceptable
- [ ] 10. Runtime/resource report captured
- [ ] 11. Reproducibility package complete

## Historical Controlled + Local Validation Context

## Current Center of Gravity

The current center of gravity is the new ChIPs realstudy bug-fix wave. The
Python tests are green, and the realstudy ChIPs Snakemake path is now local
dry-run clean. ChIPs `2.4` is now installed in `background_project`, and a tiny
WCE/control ChIPs smoke target runs locally. Full production execution is still
cluster-oriented because it requires full real-study downloads and production
reference/index assets.

## ChIPs Environment and Tiny Smoke Wave

- [x] Confirmed `background_project` did not expose `chips` at wave start.
- [x] Installed or exposed ChIPs inside `background_project`.
- [x] Added or selected a tiny local smoke config.
- [x] Re-ran realstudy compile checks.
- [x] Re-ran realstudy tests.
- [x] Ran a tiny ChIPs Snakemake dry-run.
- [x] Ran one tiny ChIPs output target, if safe.
- [x] Updated the local-runnability assessment.

## Local and Slurm ChIPs Execution Wave

- [x] Inspected current ChIPs launch/config files.
- [x] Added local full-run config for staged local assets.
- [x] Added local Snakemake launcher.
- [x] Updated local-vs-Slurm documentation.
- [x] Re-ran compile and tests.
- [x] Validated local launcher dry-run behavior.
- [x] Rechecked tiny smoke execution or output.
- [x] Logged final local and cluster run commands.

## Cluster Resume Wave

- [x] EPIC2 hilly peaks config completed on Slurm.
- [x] EPIC2 hilly plateaus config completed on Slurm.
- [x] EPIC2 hilly peaks metrics generated.
- [x] EPIC2 hilly plateaus metrics generated.
- [x] Six-config EPIC2 report README, attempt log, and per-config plots
  generated.
- [x] Realstudy failure root cause identified as package-import invocation in
  `prepare_realstudy_ingest_plan`, separate from the later Snakemake cleanup
  traceback.
- [x] Realstudy Snakemake invocation patched and local tests/dry-run rerun.
- [x] Realstudy ontology job resubmitted with the requested `64` CPUs, `50G`
  RAM, `publicgrp/low`, and `2-00:00:00` walltime.
- [ ] Realstudy resubmission finished successfully.
- [ ] Realstudy ontology metrics and plots collected.

## Broad/Narrow-Safe Peakcaller Retry Wave

- [x] New additive wave appended to all `/goals` tracking files.
- [x] Six-category shape/mode audit completed.
- [x] Shape/mode mismatch fixed in controlled configs and workflow logic.
- [x] Shape/mode mismatch documented so it cannot silently recur.
- [x] Candidate replacement peakcallers researched from primary sources.
- [x] Best-fit replacement caller selected and justified.
- [x] Replacement caller integrated into `chipseq_pipeline_v2`.
- [x] Broad-mode smoke validation passes locally.
- [x] Narrow-mode smoke validation passes locally.
- [x] Initial parameter-tuning loop completed and logged.
- [x] Pilot output namespacing added so one-run broad and narrow validations do
  not collide.
- [x] First staged batch configs created and representative narrow/broad runs
  completed locally.
- [x] Batch anomaly checks started for staged representative runs.
- [x] Six-category staged report builder updated to use canonical category labels.
- [x] Confirmed by dry-run that HOMER tagdir/peak rules are reachable for completed BAM pairs; no hidden workflow dependency bug found.
- [x] Flatearth narrow and broad 128-run HOMER staged categories completed and scored locally.
- [x] Staged HOMER combined-report wrapper added and validated on the four completed categories.
- [x] Staged HOMER headline-summary script added and validated on the four completed categories.
- [x] Staged HOMER decision-summary artifact added and validated on the four completed categories.
- [x] One-command staged HOMER artifact refresh script added and validated locally.
- [x] Staged HOMER headline wording clarified to define mean per-category best F1.
- [x] Matched MACS2 128-run control configs added for the six HOMER staged parameter sets.
- [x] Reusable staged HOMER-vs-MACS2 comparison helper added and validated on
  the currently completed paired categories.
- [x] One-command staged artifact refresh now includes the paired HOMER-vs-MACS2
  comparison outputs.
- [x] Representative MACS2 128-run control batches completed and scored for
  wavy narrow and wavy broad.
- [x] Matched wavy HOMER-vs-MACS2 control comparison completed; MACS2 is much
  stronger on the same parameter sets.
- [x] Remaining four matched MACS2 128-run control categories completed and
  scored locally, yielding the full six-category MACS2 control report.
- [x] Flatearth matched MACS2 128-run control categories completed and scored.
- [x] Matched flatearth HOMER-vs-MACS2 control comparison completed; MACS2 is
  much stronger on the same parameter sets.
- [x] Local strategy updated to treat MACS2 as the leading larger-sweep path
  pending the remaining hilly matched controls.
- [x] Long-run monitoring strategy tightened to milestone-based,
  heartbeat-first checks so `/goals` progress does not spend tokens on repeated
  near-term polling.
- [x] Remaining four 128-run HOMER staged categories completed/scored locally,
  yielding the first full six-category staged HOMER report.
- [x] Staged HOMER summary script fixed for the complete `0`-missing-category
  state and refreshed successfully.
- [x] Full six-category staged HOMER-vs-MACS2 paired comparison completed and
  refreshed from generated artifacts.
- [x] MACS2 rerun path checked for planted-shape/decoder mismatch.
- [x] Local sequential launcher added and dry-run validated for the next
  six-config `balanced_tfclean_*_288` MACS2 rerun stage.
- [~] Next larger local MACS2 rerun stage has started with a representative
  `288`-run category under the isolated launcher.
- [~] Local six-config/`1728` MACS2 wave is no longer just the original wavy
  pair: one bounded flatearth continuation root has now been launched as well.
- [x] Wavy local `288` MACS2 pair completed through scoring/report refresh and
  now contributes the first larger-scale local report outputs.
- [x] Flatearth peak local `288` MACS2 root completed through scoring/report
  refresh and now joins the larger-scale local report outputs.
- [x] Flatearth plateau local `288` MACS2 root completed through scoring/report
  refresh and now joins the larger-scale local report outputs.
- [~] Local six-config/`1728` MACS2 wave has now advanced again: the bounded
  flatearth plateau continuation root has been launched.
- [~] Local six-config/`1728` MACS2 wave has now advanced into the hilly pair
  as well: the bounded hilly peak continuation root has been launched.
- [x] Hilly peak local `288` MACS2 root completed through scoring/report
  refresh and now joins the larger-scale local report outputs.
- [x] Final hilly plateau local `288` MACS2 root completed through
  scoring/report refresh.
- [x] Authoritative all-six local `288` report and decision artifacts were
  refreshed after the final hilly plateau score pass.
- [x] Local six-config/`1728` MACS2 wave is fully scored across all six
  categories.
- [x] Guarded local scorer/report refresher added for completed
  `balanced_tfclean_*_288` MACS2 runs.
- [x] Pair-level and full-batch canonical selectors added for the local
  `balanced_tfclean_*_288` runner and refresher helpers.
- [x] Low-token local `288` progress summarizer added so heartbeat resumes can
  inspect one milestone artifact instead of repeated manual directory checks.
- [x] One-command local `288` progress refresh wrapper added so heartbeat
  checkpoints use the same selector-driven shell pattern as the runner and
  score refresher.
- [x] Full local six-config/`1728` wave status can now be refreshed from one
  artifact with explicit launch-state and next-action buckets.
- [x] Remaining not-started local `288` configs can now be selected and
  launched automatically without re-targeting active roots.
- [x] Remaining-launch helper now supports subset-driven staged continuation
  like `flatearth_pair` or `hilly_pair`.
- [x] Remaining-launch helper now supports bounded continuation with
  `--limit N`.
- [x] Progress artifact scoping fixed so subset refreshes no longer overwrite
  the authoritative `all_six` full-wave view.
- [x] Batch-level local-wave decision artifact added on top of the authoritative
  `all_six` progress table.
- [x] One-command decision refresh wrapper added for the authoritative local
  six-config wave checkpoint.
- [x] HPC realstudy continuation status refreshed in the logs.

## Completion Audit Snapshot

- [x] Staged `128` HOMER-vs-MACS2 comparison is complete and decision-quality.
- [x] Local MACS2 `288` rerun wave is complete across all six categories.
- [ ] HPC realstudy resubmission finished successfully.
- [ ] HPC realstudy ontology metrics and plots collected.
