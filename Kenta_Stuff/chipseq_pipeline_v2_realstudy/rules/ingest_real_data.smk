rule prepare_realstudy_ingest_plan:
    input:
        "manifests/study_file_manifest.csv"
    output:
        "analysis_outputs/realstudy_ingest_prep/ingest_plan.csv"
    shell:
        """
        python scripts/fetch_real_study_data.py \
          --file-manifest {input} \
          --output-dir analysis_outputs/realstudy_ingest_prep
        """
