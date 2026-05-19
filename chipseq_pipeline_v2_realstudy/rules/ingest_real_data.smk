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
        python -m scripts.fetch_real_study_data \
          --study-selection {input.selection} \
          --file-manifest {input.files} \
          --output-dir analysis_outputs/realstudy_ingest_prep
        """


def download_target_path(wc):
    return manifest_file_path(find_manifest_row(wc.study_id, wc.role))


rule download_realstudy_file:
    input:
        manifest="metadata/data_manifest.csv"
    output:
        done="analysis_outputs/realstudy_ingest_prep/download_markers/{study_id}/{role}.done"
    params:
        output_path=download_target_path
    shell:
        """
        python -m scripts.download_real_study_file \
          --manifest-csv {input.manifest} \
          --study-id {wildcards.study_id} \
          --role {wildcards.role} \
          --output-path {params.output_path}
        touch {output.done}
        """
