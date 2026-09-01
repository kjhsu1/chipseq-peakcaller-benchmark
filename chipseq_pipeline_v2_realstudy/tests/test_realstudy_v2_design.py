"""Tests for the locked empirical Realstudy v2 design and registries."""

"""Imports"""

from pathlib import Path

import pandas as pd
import yaml

from scripts.realstudy_sampling_lib import (
    audit_realstudy_v2_design,
    build_realstudy_v2_run_rows,
    fragment_rank_digest,
    target_fragments,
)


"""Functions"""


def design() -> dict:
    """Load the versioned empirical-control design."""
    root = Path(__file__).resolve().parents[1]
    return yaml.safe_load((root / "configs/realstudy_v2.yaml").read_text())["realstudy_v2"]


def test_design_expands_to_exactly_forty_four_calls():
    """The matrix must contain 42 sampled calls and two full-control anchors."""
    rows = build_realstudy_v2_run_rows(design())
    audit = audit_realstudy_v2_design(rows, design())
    assert len(rows) == 44
    assert audit["sampled_peak_calls"] == 42
    assert audit["anchor_peak_calls"] == 2
    assert len({row["run_id"] for row in rows}) == 44
    assert {row["macs2_mode"] for row in rows} == {"narrow", "broad"}


def test_locked_fragment_targets_are_exact():
    """ce11 coverage conversion must reproduce every preregistered target."""
    expected = [334288, 668576, 1337152, 2674304, 5348608, 10697216, 21394432]
    observed = [target_fragments(100286401, 150, coverage) for coverage in [0.5, 1, 2, 4, 8, 16, 32]]
    assert observed == expected


def test_rank_is_reproducible_and_seed_specific():
    """Identical rank inputs repeat exactly while a changed seed changes rank."""
    first = fragment_rank_digest("sha256-rg-qname-v1", "a" * 64, 11, "RG1\x1fq1")
    assert first == fragment_rank_digest("sha256-rg-qname-v1", "a" * 64, 11, "RG1\x1fq1")
    assert first != fragment_rank_digest("sha256-rg-qname-v1", "a" * 64, 23, "RG1\x1fq1")


def test_input_registry_matches_locked_accessions_counts_and_checksums():
    """All seven locked FASTQs must have positive counts and MD5 checksums."""
    root = Path(__file__).resolve().parents[1]
    files = pd.read_csv(root / "manifests/realstudy_v2_files.csv", dtype=str).fillna("")
    assert set(files["file_id"]) == {
        "ENCFF077SJJ", "ENCFF577WJY", "ENCFF800SHI", "ENCFF364NJX",
        "SRR1917669", "SRR1917670", "SRR1917671",
    }
    assert pd.to_numeric(files["raw_reads"]).gt(0).all()
    assert pd.to_numeric(files["compressed_bytes"]).gt(0).all()
    assert files["md5"].str.fullmatch(r"[0-9a-f]{32}").all()


def test_every_registered_publication_has_a_supported_claim():
    """No paper claim source may enter the registry without an explicit role."""
    root = Path(__file__).resolve().parents[1]
    publications = pd.read_csv(root / "manifests/realstudy_v2_publications.csv", dtype=str).fillna("")
    assert publications["citation_id"].is_unique
    assert publications["supported_claim_or_method"].str.strip().ne("").all()
    assert publications["url"].str.startswith("https://").all()
