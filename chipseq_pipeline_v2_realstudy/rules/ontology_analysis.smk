def chips_ontology_cfg(name, default=None):
    return config.get("chips", {}).get("ontology", {}).get(name, default)


def chips_ontology_score_targets():
    return [
        f"analysis_outputs/chips_ontology/{row['run_id']}/region_metrics.csv"
        for row in chips_runs()
    ]


def chips_ontology_classified_targets():
    return [
        f"analysis_outputs/chips_ontology/{row['run_id']}/classified.csv"
        for row in chips_runs()
    ]


def chips_ontology_summary_targets():
    return [
        "analysis_outputs/chips_ontology/combined_region_metrics.csv",
        "analysis_outputs/chips_ontology/summary/per_run_overall_metrics.csv",
        "analysis_outputs/chips_ontology/summary/per_ontology_metrics.csv",
        "analysis_outputs/chips_ontology/summary/control_response_by_ontology.csv",
        "analysis_outputs/chips_ontology/summary/failure_mode_metrics.csv",
        "analysis_outputs/chips_ontology/summary/ontology_f1_heatmap.csv",
        "analysis_outputs/chips_ontology/summary/ontology_f1_heatmap.png",
        "analysis_outputs/chips_ontology/summary/failure_mode_summary.md",
    ]


def chips_ontology_targets():
    return chips_ontology_classified_targets() + chips_ontology_summary_targets()


def chips_truth_bed_for_run(wc):
    study_id = find_row(wc.run_id)["study_id"]
    return f"analysis_outputs/realstudy_peakcalls/{study_id}/{study_id}_peaks.bed"


rule score_chips_ontology_regions:
    input:
        truth_bed=chips_truth_bed_for_run,
        called_bed="results_chips/{run_id}/peaks/macs2/{run_id}_peaks.bed",
        treat_bam="results_chips/{run_id}/bowtie2/treat/aligned.sorted.bam",
        treat_bai="results_chips/{run_id}/bowtie2/treat/aligned.sorted.bam.bai",
        control_bam="results_chips/{run_id}/bowtie2/con/aligned.sorted.bam",
        control_bai="results_chips/{run_id}/bowtie2/con/aligned.sorted.bam.bai",
    output:
        csv="analysis_outputs/chips_ontology/{run_id}/region_metrics.csv",
    params:
        run_table=lambda wc: config["params_table"],
        score_col=lambda wc: chips_ontology_cfg("peak_score_column", chips_cfg("peak_score_column", 5)),
    shell:
        r"""
        python scripts/score_chips_ontology_regions.py \
          --run-id {wildcards.run_id} \
          --run-table-csv {params.run_table} \
          --truth-bed {input.truth_bed} \
          --called-bed {input.called_bed} \
          --treat-bam {input.treat_bam} \
          --control-bam {input.control_bam} \
          --score-column {params.score_col} \
          --output-csv {output.csv}
        """


rule classify_chips_ontology_regions:
    input:
        csv="analysis_outputs/chips_ontology/{run_id}/region_metrics.csv",
    output:
        csv="analysis_outputs/chips_ontology/{run_id}/classified.csv",
        summary="analysis_outputs/chips_ontology/{run_id}/classified_summary.csv",
        bed="analysis_outputs/chips_ontology/{run_id}/classified.bed",
        definition="analysis_outputs/chips_ontology/{run_id}/classified_definition.json",
    params:
        prefix=lambda wc: f"analysis_outputs/chips_ontology/{wc.run_id}/classified",
    shell:
        r"""
        python -m scripts.classify_regions \
          --input-csv {input.csv} \
          --output-prefix {params.prefix}
        """


rule combine_chips_ontology_tables:
    input:
        chips_ontology_classified_targets(),
    output:
        csv="analysis_outputs/chips_ontology/combined_region_metrics.csv",
    shell:
        r"""
        python scripts/combine_csv_tables.py \
          --output-csv {output.csv} \
          {input}
        """


rule evaluate_chips_ontology:
    input:
        csv="analysis_outputs/chips_ontology/combined_region_metrics.csv",
    output:
        per_run="analysis_outputs/chips_ontology/summary/per_run_overall_metrics.csv",
        per_ontology="analysis_outputs/chips_ontology/summary/per_ontology_metrics.csv",
        control_response="analysis_outputs/chips_ontology/summary/control_response_by_ontology.csv",
        failure_modes="analysis_outputs/chips_ontology/summary/failure_mode_metrics.csv",
        heatmap_csv="analysis_outputs/chips_ontology/summary/ontology_f1_heatmap.csv",
        heatmap_png="analysis_outputs/chips_ontology/summary/ontology_f1_heatmap.png",
        summary_md="analysis_outputs/chips_ontology/summary/failure_mode_summary.md",
    params:
        outdir="analysis_outputs/chips_ontology/summary",
    shell:
        r"""
        python scripts/evaluate_by_region_ontology.py \
          --input-csv {input.csv} \
          --output-dir {params.outdir}
        """


rule chips_ontology_done:
    input:
        chips_ontology_targets()
