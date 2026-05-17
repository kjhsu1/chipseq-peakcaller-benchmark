configfile: "config.yaml"

import pandas as pd
from pathlib import Path

from scripts.realstudy_manifest_lib import (
    build_data_manifest,
    load_file_manifest,
    load_study_selection,
    validate_selection_table,
)
from scripts.realstudy_sampling_lib import build_run_table_rows, load_observed_depths


SELECTION_PATH = Path("manifests/study_selection.yaml")
FILE_MANIFEST_PATH = Path("manifests/study_file_manifest.csv")


def config_bool(name, default):
    value = config.get(name, default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"1", "true", "yes", "on"}:
            return True
        if normalized in {"0", "false", "no", "off"}:
            return False
    return bool(value)


def load_runs():
    observed_by_study = load_observed_studies()
    return build_run_table_rows(
        study_ids=list(config["study_ids"]),
        observed_by_study=observed_by_study,
        coverage_treat=list(config["coverage_treat"]),
        coverage_ctrl=list(config["coverage_ctrl"]),
        seeds=list(config["seed"]),
        fragment_lengths=list(config["fragment_length"]),
        read_lengths=list(config["read_length"]),
        aligners=list(config["aligners"]),
        peakcallers=list(config["peakcaller_list"]),
        macs2_modes=list(config["macs2_mode"]),
    )


def load_data_manifest():
    if not SELECTION_PATH.exists() or not FILE_MANIFEST_PATH.exists():
        return []
    selection_df = validate_selection_table(load_study_selection(SELECTION_PATH))
    file_df = load_file_manifest(FILE_MANIFEST_PATH)
    manifest_df = build_data_manifest(selection_df, file_df, [Path.cwd()])
    return manifest_df.fillna("").to_dict(orient="records")


def load_observed_studies():
    return load_observed_depths(Path(config["study_depths_csv"]))

RUNS = load_runs()
RUN_BY_ID = {row["run_id"]: row for row in RUNS}
DATA_MANIFEST = load_data_manifest()
OBSERVED_STUDIES = load_observed_studies()


def find_row(run_id):
    return RUN_BY_ID[run_id]


def selected_manifest_rows():
    return [row for row in DATA_MANIFEST if str(row.get("selection_status", "")).strip() == "selected"]


def selected_fastq_rows():
    rows = []
    for row in selected_manifest_rows():
        if str(row.get("file_format", "")).strip().lower() not in {"fastq", "fq"}:
            continue
        if str(row.get("needs_alignment", "")).strip().lower() != "true":
            continue
        rows.append(row)
    return rows


def selected_bam_rows():
    rows = []
    for row in selected_manifest_rows():
        if str(row.get("file_format", "")).strip().lower() != "bam":
            continue
        if str(row.get("needs_alignment", "")).strip().lower() == "true":
            continue
        rows.append(row)
    return rows


def selected_study_ids():
    return sorted({row["study_id"] for row in selected_manifest_rows()})


def find_manifest_row(study_id, role):
    for row in DATA_MANIFEST:
        if row.get("study_id") == study_id and row.get("role") == role:
            return row
    raise KeyError(f"Missing manifest row for {study_id} / {role}")


def observed_study(study_id):
    return OBSERVED_STUDIES.get(study_id, {})


def config_mapping(name):
    value = config.get(name, {})
    return value if isinstance(value, dict) else {}


def chips_mapping(name):
    value = config.get("chips", {}).get(name, {})
    return value if isinstance(value, dict) else {}


def reference_fasta_for_assembly(assembly):
    return (
        config_mapping("reference_fasta_by_assembly").get(assembly)
        or chips_mapping("reference_fasta_by_assembly").get(assembly)
        or f"references/{assembly}/genome.fa"
    )


def bowtie2_index_for_assembly(assembly):
    return (
        config_mapping("bowtie2_index_by_assembly").get(assembly)
        or chips_mapping("bowtie2_index_by_assembly").get(assembly)
        or f"indexes/{assembly}/bowtie2_index"
    )


def genome_size_for_assembly(assembly):
    return str(
        config_mapping("genome_size_by_assembly").get(assembly)
        or chips_mapping("genome_size_by_assembly").get(assembly)
        or 1000000
    )


def normalized_local_path(row):
    local_path = str(row.get("local_path", "")).strip()
    if local_path:
        return local_path
    remote_url = str(row.get("remote_url", "")).strip()
    if remote_url:
        remote_name = Path(remote_url).name
        if remote_name:
            return str(Path("data/raw") / str(row["study_id"]) / remote_name)
    return ""


def manifest_file_path(row):
    return normalized_local_path(row)


def manifest_bam_path(row):
    if str(row.get("needs_alignment", "")).strip().lower() == "true":
        return f"data/aligned/{row['study_id']}/{row['role']}/aligned.sorted.bam"
    return manifest_file_path(row)


def selected_rows_for_prefix(study_id, prefix):
    rows = []
    for row in selected_manifest_rows():
        if row["study_id"] != study_id:
            continue
        if str(row.get("role", "")).startswith(prefix):
            rows.append(row)
    return rows


def study_assembly(study_id):
    observed = observed_study(study_id)
    assembly = str(observed.get("assembly", "")).strip()
    if assembly:
        return assembly
    for row in selected_manifest_rows():
        if row["study_id"] == study_id:
            assembly = str(row.get("assembly", "")).strip()
            if assembly:
                return assembly
    return "unknown"


include: "rules/ingest_real_data.smk"
include: "rules/build_reference_intensity.smk"
include: "rules/simulation.smk"
include: "rules/alignment.smk"
include: "rules/peakcalling.smk"
include: "rules/chips_simulation.smk"


def default_targets():
    targets = [
        "analysis_outputs/realstudy_ingest_prep/ingest_plan.csv",
        config["params_table"],
        "metadata/prototype_run_table.summary.json",
    ]
    if config_bool("enable_alignment_peakcalling_targets", False):
        targets.extend(peaks_all())
    if config_bool("enable_ingest_peakcalling_targets", False):
        targets.extend(ingested_peak_targets())
    if config_bool("enable_chips_targets", False):
        targets.extend(chips_peak_targets())
    return targets


rule all:
    input:
        default_targets()
    default_target: True
