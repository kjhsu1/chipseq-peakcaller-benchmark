# Balanced 288 Config Report

This directory contains one subdirectory per completed balanced 288 config.
Each category folder includes:
- `per_run_stats.csv` produced by `scripts/peak_pr_stats.py`
- `group_summary.csv` produced by `scripts/peak_pr_stats.py`
- `plot_point_summary.csv` with one aggregated plot point per treatment/control pair
- `pr_recall_f1_vs_ctrl_coverage.png` with 3 panels and one curve per treatment coverage
- `data_info.md` summarizing the source data, swept parameters, and point counts

## Included Categories
- `balanced_tfclean_flatearth_peaks_broad_integrated_288` from `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/archived_results/balanced_tfclean_flatearth_peaks_broad_integrated_288_20260505_005251`
- `balanced_tfclean_flatearth_plateaus_broad_integrated_288` from `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/archived_results/balanced_tfclean_flatearth_plateaus_broad_integrated_288_20260505_005251`
- `balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288` from `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/archived_results/balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288_20260505_021306`
- `balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288` from `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/archived_results/balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288_20260505_005251`
- `balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288` from `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/archived_results/balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288_20260505_021306`
- `balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288` from `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/archived_results/balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288_20260505_005251`

If present, `attempt_history.log` records the known Slurm/log history for these six config runs.
