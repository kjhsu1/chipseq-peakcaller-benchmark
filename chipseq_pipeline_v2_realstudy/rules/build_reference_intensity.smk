rule make_binned_tracks:
    input:
        treat=lambda wc: f"data/reference/{wc.study_id}/treatment.bam",
        ctrl=lambda wc: f"data/reference/{wc.study_id}/control.bam"
    output:
        csv="analysis_outputs/reference_tracks/{study_id}/reference_binned.csv",
        meta="analysis_outputs/reference_tracks/{study_id}/reference_binned.json"
    shell:
        """
        python scripts/make_binned_tracks.py \
          --treatment-bam {input.treat} \
          --control-bam {input.ctrl} \
          --bin-size 25 \
          --study-id {wildcards.study_id} \
          --output-prefix analysis_outputs/reference_tracks/{wildcards.study_id}/reference_binned
        """


rule fit_reference_intensity:
    input:
        "analysis_outputs/reference_tracks/{study_id}/reference_binned.csv"
    output:
        csv="analysis_outputs/reference_tracks/{study_id}/reference_intensity.csv",
        meta="analysis_outputs/reference_tracks/{study_id}/reference_intensity.json"
    shell:
        """
        python scripts/fit_reference_intensity.py \
          --binned-track-csv {input} \
          --study-id {wildcards.study_id} \
          --output-prefix analysis_outputs/reference_tracks/{wildcards.study_id}/reference_intensity
        """
