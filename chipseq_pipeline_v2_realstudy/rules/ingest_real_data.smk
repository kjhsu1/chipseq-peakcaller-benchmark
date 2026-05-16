rule prepare_realstudy_ingest_plan:
    input:
        selection="manifests/study_selection.yaml",
        files="manifests/study_file_manifest.csv"
    output:
        plan="analysis_outputs/realstudy_ingest_prep/ingest_plan.csv",
        summary="analysis_outputs/realstudy_ingest_prep/ingest_summary.csv",
        metadata="metadata/data_manifest.csv"
    shell:
        """
        python scripts/fetch_real_study_data.py \
          --study-selection {input.selection} \
          --file-manifest {input.files} \
          --output-dir analysis_outputs/realstudy_ingest_prep
        """


rule download_realstudy_file:
    input:
        manifest="metadata/data_manifest.csv"
    output:
        touch("analysis_outputs/realstudy_ingest_prep/download_{study_id}_{role}.done")
    shell:
        """
        python scripts/download_real_study_file.py \
          --manifest-csv {input.manifest} \
          --study-id {wildcards.study_id} \
          --role {wildcards.role}
        touch {output}
        """
