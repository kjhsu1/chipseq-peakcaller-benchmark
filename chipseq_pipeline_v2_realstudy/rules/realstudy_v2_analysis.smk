def v2_anchor_peak_for_run(wc):
    study_id = V2_RUN_BY_ID[wc.run_id]["study_id"]
    return f"{v2_root()}/peaks/{study_id}__full_control_anchor/peaks.bed"


def v2_comparison_paths(filename):
    return [f"{v2_root()}/comparisons/{row['run_id']}/{filename}" for row in v2_sampled_runs()]


def v2_parent_qc_paths():
    return [
        v2_parent_qc(study_id, role)
        for study_id in V2_CONFIG.get("studies", {})
        for role in ["treatment", "control"]
    ]


def v2_stage_markers():
    return [v2_stage_marker(file_id) for file_id in v2_file_ids()]


def v2_prepare_outputs():
    root = f"{v2_root()}/tables/prepared"
    return {
        name: f"{root}/{name}.csv"
        for name in [
            "peaks", "reference_regions", "peak_overlaps", "reference_region_recovery",
            "run_metrics", "parent_libraries", "control_subsamples", "sampling_blocks",
            "software_versions", "artifacts", "validation_events",
        ]
    }


def v2_locus_inputs(wc):
    study_id = wc.study_id
    return [
        v2_parent_bam(study_id, "treatment"),
        v2_parent_bam(study_id, "control"),
        v2_control_subset(study_id, 11, 0.5),
        v2_control_subset(study_id, 11, 4),
        v2_control_subset(study_id, 11, 32),
    ]


def v2_locus_track_args(wc):
    paths = v2_locus_inputs(wc)
    labels = ["fixed_treatment", "full_control", "control_0p5x_seed11", "control_4x_seed11", "control_32x_seed11"]
    return " ".join(f"--track {label}={path}" for label, path in zip(labels, paths))


def v2_locus_peak_inputs(wc):
    study_id = wc.study_id
    run_ids = [
        f"{study_id}__full_control_anchor",
        f"{study_id}__control_0p5x__seed_11",
        f"{study_id}__control_4x__seed_11",
        f"{study_id}__control_32x__seed_11",
    ]
    return [f"{v2_root()}/peaks/{run_id}/peaks.bed" for run_id in run_ids]


def v2_locus_peak_args(wc):
    labels = ["peaks_full_control", "peaks_0p5x_seed11", "peaks_4x_seed11", "peaks_32x_seed11"]
    return " ".join(f"--peak-track {label}={path}" for label, path in zip(labels, v2_locus_peak_inputs(wc)))


rule compare_realstudy_v2_peaksets:
    input:
        query=f"{v2_root()}/peaks/{{run_id}}/peaks.bed",
        anchor=v2_anchor_peak_for_run
    output:
        peaks=f"{v2_root()}/comparisons/{{run_id}}/peaks.csv",
        reference=f"{v2_root()}/comparisons/{{run_id}}/reference_regions.csv",
        overlaps=f"{v2_root()}/comparisons/{{run_id}}/peak_overlaps.csv",
        recovery=f"{v2_root()}/comparisons/{{run_id}}/reference_region_recovery.csv",
        metrics=f"{v2_root()}/comparisons/{{run_id}}/run_metrics.csv"
    params:
        study=lambda wc: V2_RUN_BY_ID[wc.run_id]["study_id"],
        mode=lambda wc: V2_RUN_BY_ID[wc.run_id]["macs2_mode"],
        thresholds=lambda wc: " ".join(f"--threshold {value}" for value in V2_CONFIG["overlap_thresholds"]),
        outdir=lambda wc: f"{v2_root()}/comparisons/{wc.run_id}"
    shell:
        r"""
        python scripts/compare_realstudy_peaksets.py \
          --run-id {wildcards.run_id} --study-id {params.study} --mode {params.mode} \
          --query-peaks {input.query} --anchor-peaks {input.anchor} \
          {params.thresholds} --output-dir {params.outdir}
        """


rule prepare_realstudy_v2_database_tables:
    input:
        runs=f"{v2_root()}/design/run_table.csv",
        files=lambda wc: V2_CONFIG["manifests"]["files"],
        staged=v2_stage_markers(),
        metrics=v2_comparison_paths("run_metrics.csv"),
        peaks=v2_comparison_paths("peaks.csv"),
        reference=v2_comparison_paths("reference_regions.csv"),
        overlaps=v2_comparison_paths("peak_overlaps.csv"),
        recovery=v2_comparison_paths("reference_region_recovery.csv"),
        parents=v2_parent_qc_paths(),
        samplers=v2_sampling_manifests(),
        peak_calls=v2_peak_targets(),
        reference_assets=[
            path for path in [V2_CONFIG.get("reference_fasta", ""), V2_CONFIG.get("genome_annotation_gff", "")] + v2_bowtie2_index_files()
            if path
        ]
    output:
        peaks=f"{v2_root()}/tables/prepared/peaks.csv",
        files=f"{v2_root()}/tables/prepared/input_files.csv",
        reference=f"{v2_root()}/tables/prepared/reference_regions.csv",
        overlaps=f"{v2_root()}/tables/prepared/peak_overlaps.csv",
        recovery=f"{v2_root()}/tables/prepared/reference_region_recovery.csv",
        metrics=f"{v2_root()}/tables/prepared/run_metrics.csv",
        parents=f"{v2_root()}/tables/prepared/parent_libraries.csv",
        subsamples=f"{v2_root()}/tables/prepared/control_subsamples.csv",
        blocks=f"{v2_root()}/tables/prepared/sampling_blocks.csv",
        versions=f"{v2_root()}/tables/prepared/software_versions.csv",
        artifacts=f"{v2_root()}/tables/prepared/artifacts.csv",
        validation=f"{v2_root()}/tables/prepared/validation_events.csv"
    params:
        parent_args=lambda wc, input: " ".join(f"--parent {path}" for path in input.parents),
        sampler_args=lambda wc, input: " ".join(f"--sampler-manifest {path}" for path in input.samplers),
        reference_args=lambda wc, input: " ".join(
            f"--reference-asset asset_{index}={path}" for index, path in enumerate(input.reference_assets)
        ),
        root=v2_root(),
        outdir=f"{v2_root()}/tables/prepared"
    shell:
        r"""
        python scripts/prepare_realstudy_v2_database_tables.py \
          --run-table {input.runs} --files-manifest {input.files} --comparison-root {params.root}/comparisons \
          --peak-root {params.root}/peaks {params.parent_args} {params.sampler_args} {params.reference_args} \
          --output-dir {params.outdir}
        """


rule annotate_realstudy_v2_regions:
    input:
        reference=f"{v2_root()}/tables/prepared/reference_regions.csv",
        recovery=f"{v2_root()}/tables/prepared/reference_region_recovery.csv",
        runs=f"{v2_root()}/design/run_table.csv",
        gff=lambda wc: V2_CONFIG["genome_annotation_gff"]
    output:
        reference=f"{v2_root()}/tables/analysis/reference_regions.csv",
        stratified=f"{v2_root()}/tables/analysis/stratified_metrics.csv"
    shell:
        r"""
        python scripts/annotate_realstudy_v2_regions.py \
          --reference-regions {input.reference} --recovery {input.recovery} \
          --run-table {input.runs} --gff {input.gff} \
          --annotated-regions {output.reference} --stratified-metrics {output.stratified}
        """


rule compare_realstudy_v2_seed_pairs:
    input:
        runs=f"{v2_root()}/design/run_table.csv",
        peaks=v2_peak_targets()
    output:
        metrics=f"{v2_root()}/tables/analysis/seed_pair_metrics.csv"
    params:
        root=v2_root()
    shell:
        r"""
        python scripts/compare_realstudy_v2_seed_pairs.py \
          --run-table {input.runs} --peak-root {params.root}/peaks --output {output.metrics}
        """


rule analyze_realstudy_v2_sufficiency:
    input:
        config="configs/realstudy_v2.yaml",
        runs=f"{v2_root()}/design/run_table.csv",
        metrics=f"{v2_root()}/tables/prepared/run_metrics.csv",
        parents=f"{v2_root()}/tables/prepared/parent_libraries.csv"
    output:
        seed=f"{v2_root()}/tables/analysis/seed_level_metrics.csv",
        depth=f"{v2_root()}/tables/analysis/depth_summary.csv",
        adjacent=f"{v2_root()}/tables/analysis/adjacent_depth_gains.csv",
        decisions=f"{v2_root()}/tables/analysis/enough_control_decisions.csv",
        overlap_sensitivity=f"{v2_root()}/tables/analysis/overlap_sensitivity_summary.csv",
        policy_sensitivity=f"{v2_root()}/tables/analysis/policy_sensitivity.csv",
        augmented_metrics=f"{v2_root()}/tables/analysis/run_metrics_augmented.csv"
    params:
        outdir=f"{v2_root()}/tables/analysis"
    shell:
        r"""
        python scripts/realstudy_v2_analysis.py \
          --config {input.config} --run-table {input.runs} --run-metrics {input.metrics} \
          --parent-libraries {input.parents} --output-dir {params.outdir}
        """


rule build_realstudy_v2_locus_tracks:
    input:
        reference=f"{v2_root()}/tables/analysis/reference_regions.csv",
        tracks=v2_locus_inputs,
        peaks=v2_locus_peak_inputs
    output:
        csv=f"{v2_root()}/tables/analysis/loci/{{study_id}}.csv"
    params:
        track_args=v2_locus_track_args,
        peak_args=v2_locus_peak_args
    shell:
        r"""
        python scripts/build_realstudy_v2_locus_tracks.py \
          --study-id {wildcards.study_id} --reference-regions {input.reference} \
          {params.track_args} {params.peak_args} --output {output.csv}
        """


rule combine_realstudy_v2_locus_tracks:
    input:
        [f"{v2_root()}/tables/analysis/loci/{study_id}.csv" for study_id in V2_CONFIG.get("studies", {})]
    output:
        csv=f"{v2_root()}/tables/analysis/representative_locus_tracks.csv"
    shell:
        r"""
        python scripts/combine_csv_tables.py --output-csv {output.csv} {input}
        """


rule build_realstudy_v2_database:
    input:
        publications=lambda wc: V2_CONFIG["manifests"]["publications"],
        experiments=lambda wc: V2_CONFIG["manifests"]["experiments"],
        files=f"{v2_root()}/tables/prepared/input_files.csv",
        runs=f"{v2_root()}/design/run_table.csv",
        parents=f"{v2_root()}/tables/prepared/parent_libraries.csv",
        versions=f"{v2_root()}/tables/prepared/software_versions.csv",
        blocks=f"{v2_root()}/tables/prepared/sampling_blocks.csv",
        subsamples=f"{v2_root()}/tables/prepared/control_subsamples.csv",
        artifacts=f"{v2_root()}/tables/prepared/artifacts.csv",
        validation=f"{v2_root()}/tables/prepared/validation_events.csv",
        peaks=f"{v2_root()}/tables/prepared/peaks.csv",
        reference=f"{v2_root()}/tables/analysis/reference_regions.csv",
        overlaps=f"{v2_root()}/tables/prepared/peak_overlaps.csv",
        recovery=f"{v2_root()}/tables/prepared/reference_region_recovery.csv",
        metrics=f"{v2_root()}/tables/analysis/run_metrics_augmented.csv",
        seed_pairs=f"{v2_root()}/tables/analysis/seed_pair_metrics.csv",
        overlap_sensitivity=f"{v2_root()}/tables/analysis/overlap_sensitivity_summary.csv",
        policy_sensitivity=f"{v2_root()}/tables/analysis/policy_sensitivity.csv",
        stratified=f"{v2_root()}/tables/analysis/stratified_metrics.csv",
        decisions=f"{v2_root()}/tables/analysis/enough_control_decisions.csv",
        figures=f"{v2_root()}/figures/.all_six_complete"
    output:
        database=f"{v2_root()}/database/realstudy_v2.sqlite",
        manifest=f"{v2_root()}/database/csv_exports_manifest.csv"
    params:
        exports=f"{v2_root()}/database/csv_exports",
        figures=f"{v2_root()}/figures"
    shell:
        r"""
        python scripts/build_realstudy_v2_database.py \
          --publications {input.publications} --experiments {input.experiments} --files {input.files} \
          --run-table {input.runs} --mark-runs-complete \
          --table parent_libraries={input.parents} --table software_versions={input.versions} \
          --table sampling_blocks={input.blocks} --table control_subsamples={input.subsamples} \
          --table artifacts={input.artifacts} --table validation_events={input.validation} \
          --table peaks={input.peaks} --table reference_regions={input.reference} \
          --table peak_overlaps={input.overlaps} --table reference_region_recovery={input.recovery} \
          --table run_metrics={input.metrics} --table seed_pair_metrics={input.seed_pairs} \
          --table stratified_metrics={input.stratified} --table enough_control_decisions={input.decisions} \
          --artifact-root {params.figures} \
          --output {output.database} --export-dir {params.exports} --export-manifest {output.manifest}
        """


rule render_realstudy_v2_figures:
    input:
        config="configs/realstudy_v2.yaml",
        audit=lambda wc: V2_CONFIG["manifests"]["selection_audit"],
        files=lambda wc: V2_CONFIG["manifests"]["files"],
        seed=f"{v2_root()}/tables/analysis/seed_level_metrics.csv",
        depth=f"{v2_root()}/tables/analysis/depth_summary.csv",
        adjacent=f"{v2_root()}/tables/analysis/adjacent_depth_gains.csv",
        seed_pairs=f"{v2_root()}/tables/analysis/seed_pair_metrics.csv",
        overlap_sensitivity=f"{v2_root()}/tables/analysis/overlap_sensitivity_summary.csv",
        policy_sensitivity=f"{v2_root()}/tables/analysis/policy_sensitivity.csv",
        stratified=f"{v2_root()}/tables/analysis/stratified_metrics.csv",
        loci=f"{v2_root()}/tables/analysis/representative_locus_tracks.csv"
    output:
        touch(f"{v2_root()}/figures/.all_six_complete")
    params:
        outdir=f"{v2_root()}/figures"
    shell:
        r"""
        python scripts/render_realstudy_v2_figures.py \
          --config {input.config} --selection-audit {input.audit} --files {input.files} \
          --seed-level {input.seed} --depth-summary {input.depth} --adjacent-gains {input.adjacent} --seed-pairs {input.seed_pairs} \
          --overlap-sensitivity {input.overlap_sensitivity} --policy-sensitivity {input.policy_sensitivity} \
          --stratified-metrics {input.stratified} --locus-tracks {input.loci} --output-dir {params.outdir}
        """


rule validate_realstudy_v2_integrity:
    input:
        database=f"{v2_root()}/database/realstudy_v2.sqlite",
        manifest=f"{v2_root()}/database/csv_exports_manifest.csv",
        runs=f"{v2_root()}/design/run_table.csv",
        figures=f"{v2_root()}/figures/.all_six_complete"
    output:
        report=f"{v2_root()}/validation/final_integrity_report.json"
    params:
        figures=f"{v2_root()}/figures"
    shell:
        r"""
        python scripts/validate_realstudy_v2.py \
          --database {input.database} --run-table {input.runs} \
          --figure-dir {params.figures} --export-manifest {input.manifest} \
          --report {output.report}
        """
