"""Tests for real-study file downloading helpers."""

from scripts.download_real_study_file import choose_transfer_mode, filename_from_url


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
