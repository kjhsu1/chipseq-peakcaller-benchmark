rule make_binned_tracks:
    input:
        treat="data/reference/treatment.bam",
        ctrl="data/reference/control.bam"
    output:
        "analysis_outputs/reference_tracks/reference_binned.csv"
    shell:
        """
        python scripts/make_binned_tracks.py \
          --treatment-bam {input.treat} \
          --control-bam {input.ctrl} \
          --bin-size 25 \
          --output-prefix analysis_outputs/reference_tracks/reference_binned
        """


rule fit_reference_intensity:
    input:
        "analysis_outputs/reference_tracks/reference_binned.csv"
    output:
        "analysis_outputs/reference_tracks/reference_intensity.csv"
    shell:
        """
        python scripts/fit_reference_intensity.py \
          --binned-track-csv {input} \
          --output-prefix analysis_outputs/reference_tracks/reference_intensity
        """
