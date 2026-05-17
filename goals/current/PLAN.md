# Plan

## ChIPs Ontology Analysis and Cluster Sweep Wave

1. Inspect current ChIPs realstudy workflow, ontology helper scripts, and latest
   successful Slurm submission patterns. Done.
2. Before implementation, write a scientific validity report that evaluates the
   experimental method statistically, logically, and biologically, with sources.
   Done.
3. Use the report to adjust the experimental plan if needed, especially around
   what ChIPs can and cannot provide as ground truth. Done: the downstream
   analysis is framed as empirical truth-template recovery rather than exact
   intrinsic-PMF recovery.
4. Add an ontology downstream step to the realstudy ChIPs Snakemake workflow.
   Done.
5. Add dry-run-testable targets for ontology outputs after ChIPs MACS2 peak
   calls. Done.
6. Add 128-core, 2TB-RAM Slurm and shell submission scripts for the cluster,
   matching the latest successful sweep script structure as closely as possible.
   Done.
7. Validate with compile checks, tests, and Snakemake dry-runs, then record
   proof and limitations. Done.

## ChIPs Realstudy Bug-Fix Wave

1. Reproduce and document the current ChIPs realstudy failure points without
   changing code. Done.
2. Repair the ingest-to-alignment handoff so downloaded FASTQs are real
   workflow outputs, not out-of-band assumptions.
3. Split processed-BAM studies from FASTQ-needs-alignment studies so ENCODE BAM
   inputs do not get forced through `data/aligned/...` paths meant for FASTQ
   ingest.
4. Remove stale parse-time metadata assumptions in `Snakefile.py` so regenerated
   manifests and run tables do not silently disagree with the DAG.
5. Unify local path configuration for references and Bowtie2 indexes so local
   dry-runs and local execution use the same override mechanism across ingest
   alignment and ChIPs simulation.
6. Make peak-calling genome size assembly-aware instead of hard-coding one worm
   scale value for every study.
7. Align ChIPs model-learning inputs with the real-study peak BED generation
   logic so the BAM used for `chips learn` matches the replicate structure used
   to define the training peaks.
8. Tighten manifest/local-file bookkeeping so `downloaded` and `local_exists`
   reflect the current checkout accurately and do not overstate local readiness.
9. Future-proof run IDs against parameter-axis growth so later expansions do not
   create silent collisions.
10. Re-run the realstudy test suite and local Snakemake dry-runs, then reassess
   whether the ChIPs path is truly locally runnable or still cluster-only.

1. Record the current goal and initialize progress/log tracking. Done.
2. Document the controlled-sweep bug-fix state as prerequisite context. Done.
3. Write simulator compatibility criteria, compare ChIPs/ChIPulate/isChIP/ChIPsim,
   and choose the simulator for the realistic benchmark path.
   Done.
4. Revise the 11-step plan around ChIPs and the ontology handoff. Done.
5. Add ChIPs workflow config/rules/Slurm launch assets in
   `chipseq_pipeline_v2_realstudy`. Done.
6. Add matched EPIC2 sweep config and Slurm/shell launch assets in
   `chipseq_pipeline_v2`. Done.
7. Add `/goals` operating instructions to root `AGENTS.md`. Done.
8. Prune low-utility active scripts and add a plain-language script index. Done.
9. Run fast checks, tests, and dry-runs; record proof and limitations. Done in
   the local `background_project` Conda environment.

## Next Cluster Step

Run these from the cluster once the branch is available there:

```bash
cd chipseq_pipeline_v2
source ../snakemake_stuff/setup.sh
snakemake -s Snakefile.py --configfile configs/epic2_tfclean_realistic_peaks_wavy_narrow_integrated_288.yaml --dry-run
bash slurm/submit_epic2_tfclean_1728_series.sh
```

```bash
cd chipseq_pipeline_v2_realstudy
source ../snakemake_stuff/setup.sh
snakemake -s Snakefile.py --config enable_chips_targets=true --dry-run
bash slurm/submit_chips_realsim.sh
```

## Local Validation Command Pattern

On this Mac, use:

```bash
eval "$(/opt/anaconda3/bin/conda shell.zsh hook)"
conda activate background_project
HOME=/private/tmp/chipseq_snakemake_home snakemake -s Snakefile.py --dry-run
```

## ChIPs Environment and Tiny Smoke Wave

1. Confirm whether `background_project` exposes `chips`. Done: it did not at
   the start of this wave.
2. Install or expose ChIPs inside `background_project` so workflow commands do
   not depend on another Conda environment. Done.
3. Add or reuse a tiny local smoke config that resolves one ChIPs run only, with
   very small coverage settings and no large sweep targets. Done.
4. Validate in a tight loop:
   - `python -m py_compile scripts/*.py`
   - `pytest tests`
   - Snakemake dry-run for the tiny ChIPs config
   - one tiny ChIPs-controlled output target if the installed binary is usable
   Done.
5. Record proof, limitations, and whether the realstudy pipeline is locally
   runnable beyond dry-run. Done.

## Local and Slurm ChIPs Execution Wave

1. Inspect the current ChIPs configs, docs, and launch scripts. Done.
2. Add a local execution config that uses local staged asset paths instead of
   cluster production paths. Done.
3. Add a small local launcher for dry-run and full-run modes without Slurm.
   Done.
4. Update docs so local Snakemake, tiny smoke, and cluster Slurm commands are
   all explicit and separate. Done.
5. Validate compile/tests plus local launcher dry-runs and tiny smoke execution.
   Done.
