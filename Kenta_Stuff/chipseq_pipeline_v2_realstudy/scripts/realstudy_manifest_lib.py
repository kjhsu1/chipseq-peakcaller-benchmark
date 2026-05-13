"""Helpers for real-study manifest normalization and validation."""

"""Imports"""

from pathlib import Path
from typing import Dict, List

import pandas as pd


"""Functions"""


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
            current[key.strip()] = value.strip()
            continue
        if current is None:
            continue
        key, value = stripped.split(":", 1)
        current[key.strip()] = value.strip()
    if current is not None:
        studies.append(current)
    return studies


def validate_selection_table(df: pd.DataFrame) -> pd.DataFrame:
    """Add explicit gating flags for study eligibility."""
    validated = df.copy()
    validated["eligible_pairing"] = validated["has_treatment"].astype(str).eq("true") & validated["has_control"].astype(str).eq("true")
    validated["eligible_files"] = validated["processed_bam_available"].astype(str).eq("true")
    validated["usable_for_realism_scorecard"] = validated["usable_for_realism_scorecard"].astype(str)
    validated["usable_for_realstudy_simulation"] = validated["usable_for_realstudy_simulation"].astype(str)
    return validated
