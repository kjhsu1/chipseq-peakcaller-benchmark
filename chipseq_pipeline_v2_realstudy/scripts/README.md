# Active Realstudy Scripts Index

This folder contains scripts used by the active realstudy/realism-oriented
benchmark. The workflow should stay simple: add a new script only when it makes
the Snakemake path easier to verify or avoids meaningful repeated logic.

## Manifests And Ingestion

- `realstudy_manifest_lib.py`: validates and normalizes study metadata.
- `fetch_real_study_data.py`: builds ingest planning tables and manifest rows.
- `download_real_study_file.py`: resumable file-level download helper.

## Reference Tracks And Sampling

- `make_binned_tracks.py`: creates binned treatment/control signal tracks.
- `fit_reference_intensity.py`: estimates reference intensity tracks and JSON
  sidecars without calling them truth.
- `sample_reads_from_intensity.py`: generates prototype run tables and
  sampling-depth labels.
- `realstudy_sampling_lib.py`: helpers for run-table generation and depth labels.

## Ontology And Evaluation

- `classify_regions.py`: classifies windows into ontology labels.
- `ontology_lib.py`: ontology threshold helpers.
- `evaluate_by_region_ontology.py`: aggregates performance by ontology class.
- `build_final_professor_summary.py`: assembles the final professor-facing
  summary directory.
