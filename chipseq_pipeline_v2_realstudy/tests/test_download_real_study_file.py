"""Tests for real-study file downloading helpers."""

from scripts.download_real_study_file import (
    choose_transfer_mode,
    ensure_required_local_file,
    filename_from_url,
    should_trust_existing_file,
)


def test_filename_from_url_uses_path_tail():
    """The output filename should come from the URL path."""
    assert (
        filename_from_url("https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR191/000/SRR1917670/SRR1917670.fastq.gz")
        == "SRR1917670.fastq.gz"
    )


def test_choose_transfer_mode_skips_matching_complete_file():
    """A complete local file should be reused."""
    assert choose_transfer_mode(100, 100, "bytes", False) == "skip"


def test_choose_transfer_mode_resumes_partial_range_capable_file():
    """A partial local file should resume when the server supports byte ranges."""
    assert choose_transfer_mode(25, 100, "bytes", False) == "resume"


def test_choose_transfer_mode_restarts_when_server_lacks_ranges():
    """A partial file should restart if resume is not supported."""
    assert choose_transfer_mode(25, 100, "none", False) == "restart"


def test_should_trust_existing_file_when_remote_metadata_is_unavailable():
    """A non-empty local file should be reused when the remote cannot be queried."""
    assert should_trust_existing_file(25, None, False) is True


def test_should_not_trust_existing_file_when_overwriting():
    """Overwrite requests should ignore an existing local file."""
    assert should_trust_existing_file(25, None, True) is False


def test_ensure_required_local_file_fails_on_missing_file(tmp_path):
    """Local-only mode should fail before trying a network transfer."""
    missing_path = tmp_path / "missing.fastq.gz"
    try:
        ensure_required_local_file(missing_path, 0)
    except SystemExit as exc:
        assert "Missing required pre-staged input" in str(exc)
    else:
        raise AssertionError("Expected missing local input to raise SystemExit")
