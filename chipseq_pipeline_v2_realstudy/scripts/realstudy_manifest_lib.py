"""Helpers for real-study manifest normalization and validation."""

"""Imports"""

from pathlib import Path
from typing import Dict, Iterable, List
from urllib.parse import urlparse

import pandas as pd


"""Functions"""


SELECTION_COLUMNS = [
    "study_id",
    "source",
    "accession",
    "control_accession",
    "signal_class",
    "has_treatment",
    "has_control",
    "processed_bam_available",
    "usable_for_realism_scorecard",
    "usable_for_realstudy_simulation",
    "selection_status",
    "selection_reason",
]

FILE_COLUMNS = [
    "study_id",
    "role",
    "accession",
    "file_format",
    "assembly",
    "local_path",
    "remote_url",
    "needs_alignment",
    "selection_status",
    "selection_reason",
]

FASTQ_FORMATS = {"fastq", "fq"}


def normalize_scalar(value: str) -> str:
    """Normalize simple YAML scalar text."""
    cleaned = value.strip()
    if cleaned.startswith(("'", '"')) and cleaned.endswith(("'", '"')) and len(cleaned) >= 2:
        cleaned = cleaned[1:-1]
    return cleaned


def filename_from_remote_url(remote_url: str) -> str:
    """Extract a filename from a remote URL path."""
    parsed = urlparse(str(remote_url).strip())
    return Path(parsed.path).name


def parse_simple_yaml_list(path: Path) -> List[Dict[str, str]]:
    """Parse the limited study-selection YAML format used in this prototype."""
    studies: List[Dict[str, str]] = []
    current: Dict[str, str] | None = None
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.rstrip()
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped == "studies:":
            continue
        if stripped.startswith("- "):
            if current is not None:
                studies.append(current)
            current = {}
            key, value = stripped[2:].split(":", 1)
            current[key.strip()] = normalize_scalar(value)
            continue
        if current is None:
            continue
        key, value = stripped.split(":", 1)
        current[key.strip()] = normalize_scalar(value)
    if current is not None:
        studies.append(current)
    return studies


def load_study_selection(path: Path) -> pd.DataFrame:
    """Load and normalize the study selection manifest."""
    df = pd.DataFrame(parse_simple_yaml_list(path))
    for column in SELECTION_COLUMNS:
        if column not in df.columns:
            df[column] = ""
    return df[SELECTION_COLUMNS].copy()


def normalize_bool_text(value: object) -> str:
    """Normalize booleans and unknowns to explicit lowercase strings."""
    normalized = str(value).strip().lower()
    if normalized in {"true", "false", "unknown"}:
        return normalized
    if normalized in {"1", "yes", "y"}:
        return "true"
    if normalized in {"0", "no", "n"}:
        return "false"
    return "unknown" if not normalized else normalized


def validate_selection_table(df: pd.DataFrame) -> pd.DataFrame:
    """Add explicit gating flags for study eligibility."""
    validated = df.copy()
    if "selection_status" not in validated.columns:
        validated["selection_status"] = ""
    for column in [
        "has_treatment",
        "has_control",
        "processed_bam_available",
        "usable_for_realism_scorecard",
        "usable_for_realstudy_simulation",
    ]:
        validated[column] = validated[column].map(normalize_bool_text)
    validated["eligible_pairing"] = validated["has_treatment"].eq("true") & validated["has_control"].eq("true")
    validated["eligible_files"] = validated["processed_bam_available"].eq("true")
    validated["eligible_realism_reference"] = validated["usable_for_realism_scorecard"].eq("true")
    validated["eligible_realstudy_simulation"] = validated["usable_for_realstudy_simulation"].eq("true")
    validated["selection_status"] = validated["selection_status"].fillna("").astype(str)
    return validated


def load_file_manifest(path: Path) -> pd.DataFrame:
    """Load and normalize the file manifest."""
    df = pd.read_csv(path).fillna("")
    for column in FILE_COLUMNS:
        if column not in df.columns:
            df[column] = ""
    df["needs_alignment"] = df["needs_alignment"].map(normalize_bool_text)
    return df[FILE_COLUMNS].copy()


def build_data_manifest(
    selection_df: pd.DataFrame,
    file_df: pd.DataFrame,
    local_search_roots: Iterable[Path] | None = None,
) -> pd.DataFrame:
    """Build a resumable file-level ingest manifest."""
    merged = file_df.merge(
        selection_df[
            [
                "study_id",
                "source",
                "signal_class",
                "selection_status",
                "selection_reason",
                "eligible_pairing",
                "eligible_files",
                "eligible_realism_reference",
                "eligible_realstudy_simulation",
            ]
        ],
        on="study_id",
        how="left",
        suffixes=("", "_study"),
    )
    merged["download_status"] = "not_started"
    merged["metadata_status"] = merged["selection_status"].replace(
        {
            "selected": "locally_grounded",
            "pending_external_fetch": "needs_live_resolution",
            "fallback": "fallback_not_promoted",
        }
    ).fillna("unknown")
    def infer_ingest_mode(row: pd.Series) -> str:
        file_format = str(row.get("file_format", "")).lower()
        remote_url = str(row.get("remote_url", "")).strip()
        if file_format == "bam":
            return "download_processed_bam"
        if file_format in FASTQ_FORMATS and remote_url:
            return "download_fastq"
        return "resolve_metadata_first"

    merged["ingest_mode"] = merged.apply(infer_ingest_mode, axis=1)
    merged["failure_scope"] = "per_file"
    merged["local_exists"] = False
    roots = list(local_search_roots or [])
    if roots:
        for idx, row in merged.iterrows():
            local_path = str(row.get("local_path", "")).strip()
            if not local_path:
                remote_name = filename_from_remote_url(str(row.get("remote_url", "")))
                if remote_name:
                    local_path = str(Path("data/raw") / str(row["study_id"]) / remote_name)
                    merged.at[idx, "local_path"] = local_path
            if local_path:
                local_exists = any(
                    (root / local_path).exists() for root in roots if not Path(local_path).is_absolute()
                )
                if Path(local_path).is_absolute():
                    local_exists = Path(local_path).exists()
                merged.at[idx, "local_exists"] = local_exists
                if local_exists:
                    merged.at[idx, "download_status"] = "downloaded"
    return merged
