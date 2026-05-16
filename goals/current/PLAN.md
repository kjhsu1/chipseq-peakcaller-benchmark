# Plan

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
