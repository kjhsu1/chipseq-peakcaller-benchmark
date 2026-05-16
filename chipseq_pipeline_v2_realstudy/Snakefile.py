configfile: "config.yaml"

import pandas as pd
from pathlib import Path


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
    path = config["params_table"]
    if not Path(path).exists():
        return []
    return pd.read_csv(path).to_dict(orient="records")


def load_data_manifest():
    path = Path("metadata/data_manifest.csv")
    if not path.exists():
        return []
    return pd.read_csv(path).fillna("").to_dict(orient="records")


def load_observed_studies():
    path = Path(config["study_depths_csv"])
    if not path.exists():
        return {}
    df = pd.read_csv(path).fillna("")
    return {row["study_id"]: row for row in df.to_dict(orient="records")}

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


def selected_study_ids():
    return sorted({row["study_id"] for row in selected_manifest_rows()})


def find_manifest_row(study_id, role):
    for row in DATA_MANIFEST:
        if row.get("study_id") == study_id and row.get("role") == role:
            return row
    raise KeyError(f"Missing manifest row for {study_id} / {role}")


def observed_study(study_id):
    return OBSERVED_STUDIES.get(study_id, {})


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
