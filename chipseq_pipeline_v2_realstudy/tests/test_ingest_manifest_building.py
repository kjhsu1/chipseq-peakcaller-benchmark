"""Tests for realstudy ingest bookkeeping."""

"""Imports"""

from pathlib import Path

import pandas as pd

from scripts.realstudy_manifest_lib import build_data_manifest


def test_build_data_manifest_adds_bookkeeping_columns(tmp_path: Path):
    """Data manifest rows should get resumable ingest bookkeeping fields."""
    selection = pd.DataFrame(
        [
            {
                "study_id": "study1",
                "source": "ENCODE",
                "signal_class": "narrow",
                "selection_status": "selected",
                "selection_reason": "grounded",
                "eligible_pairing": True,
                "eligible_files": True,
                "eligible_realism_reference": True,
                "eligible_realstudy_simulation": True,
            }
        ]
    )
    files = pd.DataFrame(
        [
            {
                "study_id": "study1",
                "role": "treatment",
                "accession": "ENCFF1",
                "file_format": "bam",
                "assembly": "ce11",
                "local_path": "",
                "remote_url": "https://example.org",
                "needs_alignment": "false",
                "selection_status": "selected",
                "selection_reason": "grounded",
            }
        ]
    )
    manifest = build_data_manifest(selection, files, [tmp_path])
    assert manifest.iloc[0]["download_status"] == "not_started"
    assert manifest.iloc[0]["ingest_mode"] == "download_processed_bam"
    assert manifest.iloc[0]["failure_scope"] == "per_file"


def test_build_data_manifest_marks_fastq_with_remote_url_as_downloadable(tmp_path: Path):
    """FASTQ rows with remote URLs should become download-ready."""
    selection = pd.DataFrame(
        [
            {
                "study_id": "study1",
                "source": "GEO",
                "signal_class": "broad",
                "selection_status": "selected",
                "selection_reason": "grounded",
                "eligible_pairing": True,
                "eligible_files": False,
                "eligible_realism_reference": True,
                "eligible_realstudy_simulation": True,
            }
        ]
    )
    files = pd.DataFrame(
        [
            {
                "study_id": "study1",
                "role": "treatment",
                "accession": "SRR1",
                "file_format": "fastq",
                "assembly": "ce11",
                "local_path": "",
                "remote_url": "https://example.org/run.fastq.gz",
                "needs_alignment": "true",
                "selection_status": "selected",
                "selection_reason": "grounded",
            }
        ]
    )
    manifest = build_data_manifest(selection, files, [tmp_path])
    assert manifest.iloc[0]["ingest_mode"] == "download_fastq"


def test_build_data_manifest_infers_local_path_and_downloaded_status(tmp_path: Path):
    """A downloaded file under data/raw should be discovered automatically."""
    selection = pd.DataFrame(
        [
            {
                "study_id": "study1",
                "source": "GEO",
                "signal_class": "broad",
                "selection_status": "selected",
                "selection_reason": "grounded",
                "eligible_pairing": True,
                "eligible_files": False,
                "eligible_realism_reference": True,
                "eligible_realstudy_simulation": True,
            }
        ]
    )
    study_dir = tmp_path / "data" / "raw" / "study1"
    study_dir.mkdir(parents=True)
    downloaded = study_dir / "run.fastq.gz"
    downloaded.write_bytes(b"abc")
    files = pd.DataFrame(
        [
            {
                "study_id": "study1",
                "role": "treatment",
                "accession": "SRR1",
                "file_format": "fastq",
                "assembly": "ce11",
                "local_path": "",
                "remote_url": "https://example.org/run.fastq.gz",
                "needs_alignment": "true",
                "selection_status": "selected",
                "selection_reason": "grounded",
            }
        ]
    )
    manifest = build_data_manifest(selection, files, [tmp_path])
    assert manifest.iloc[0]["local_path"] == "data/raw/study1/run.fastq.gz"
    assert bool(manifest.iloc[0]["local_exists"]) is True
    assert manifest.iloc[0]["download_status"] == "downloaded"
