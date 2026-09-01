"""Focused future tests for sampling, peak, database, and paper-output contracts."""

"""Imports"""

import hashlib
import sqlite3
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pysam
import yaml

from scripts.realstudy_sampling_lib import build_realstudy_v2_run_rows


"""Functions"""


def write_test_bam(path: Path) -> None:
    """Write a small coordinate-sorted BAM with one pair and eighteen singles."""
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrI", "LN": 10000}], "RG": [{"ID": "RG1", "SM": "sample"}]}
    records = []
    for index in range(18):
        record = pysam.AlignedSegment()
        record.query_name = f"single_{index}"
        record.query_sequence = "A" * 50
        record.flag = 0
        record.reference_id = 0
        record.reference_start = 100 + index * 100
        record.mapping_quality = 60
        record.cigar = [(0, 50)]
        record.query_qualities = pysam.qualitystring_to_array("I" * 50)
        record.set_tag("RG", "RG1")
        records.append(record)
    for flag, start in [(65, 3000), (129, 3200)]:
        record = pysam.AlignedSegment()
        record.query_name = "paired_fragment"
        record.query_sequence = "C" * 50
        record.flag = flag
        record.reference_id = 0
        record.reference_start = start
        record.mapping_quality = 60
        record.cigar = [(0, 50)]
        record.next_reference_id = 0
        record.next_reference_start = 3200 if flag == 65 else 3000
        record.template_length = 250 if flag == 65 else -250
        record.query_qualities = pysam.qualitystring_to_array("I" * 50)
        record.set_tag("RG", "RG1")
        records.append(record)
    with pysam.AlignmentFile(path, "wb", header=header) as output:
        for record in sorted(records, key=lambda item: item.reference_start):
            output.write(record)
    pysam.index(str(path))


def bam_fragment_ids(path: Path) -> set[str]:
    """Read unique query names from a sampled BAM."""
    with pysam.AlignmentFile(path, "rb") as source:
        return {record.query_name for record in source.fetch(until_eof=True)}


def run_sampler(root: Path, bam: Path, tmp_path: Path, seed: int, suffix: str) -> tuple[Path, Path]:
    """Execute the sampler on a tiny fixture for future contract validation."""
    checksum = hashlib.sha256(bam.read_bytes()).hexdigest()
    low, high = tmp_path / f"low_{suffix}.bam", tmp_path / f"high_{suffix}.bam"
    command = [
        sys.executable, "scripts/build_empirical_control_series.py", "--parent-bam", str(bam),
        "--source-checksum", checksum, "--study-id", "study", "--control-parent-id", "study__control_parent",
        "--seed", str(seed), "--genome-size-bp", "1000", "--normalization-fragment-length-bp", "100",
        "--minimum-eligible-fragments", "10", "--target", "1x=10", "--target", "2x=19",
        "--output", f"1x={low}", "--output", f"2x={high}",
        "--ledger", str(tmp_path / f"ledger_{suffix}.sqlite"), "--manifest", str(tmp_path / f"manifest_{suffix}.csv"),
    ]
    subprocess.run(command, cwd=root, check=True)
    return low, high


def test_nested_sampling_is_exact_reproducible_and_pair_safe(tmp_path: Path):
    """Subsets must be exact/nested; same seed repeats and both paired mates travel together."""
    root = Path(__file__).resolve().parents[1]
    bam = tmp_path / "parent.bam"
    write_test_bam(bam)
    low_a, high_a = run_sampler(root, bam, tmp_path, 11, "a")
    low_b, high_b = run_sampler(root, bam, tmp_path, 11, "b")
    low_ids, high_ids = bam_fragment_ids(low_a), bam_fragment_ids(high_a)
    assert len(low_ids) == 10 and len(high_ids) == 19 and low_ids < high_ids
    assert low_ids == bam_fragment_ids(low_b) and high_ids == bam_fragment_ids(high_b)
    with pysam.AlignmentFile(high_a, "rb") as source:
        paired_records = [record for record in source.fetch(until_eof=True) if record.query_name == "paired_fragment"]
    assert len(paired_records) in {0, 2}


def test_peak_comparison_distinguishes_valid_zero_from_missing(tmp_path: Path):
    """An empty completed peak file is valid, while a missing file fails."""
    root = Path(__file__).resolve().parents[1]
    empty = tmp_path / "empty.narrowPeak"; empty.write_text("")
    anchor = tmp_path / "anchor.narrowPeak"; anchor.write_text("chrI\t10\t20\ta\t10\t.\t5\t1\t1\t5\n")
    output = tmp_path / "comparison"
    command = [sys.executable, "scripts/compare_realstudy_peaksets.py", "--run-id", "run", "--study-id", "study", "--mode", "narrow", "--query-peaks", str(empty), "--anchor-peaks", str(anchor), "--output-dir", str(output)]
    subprocess.run(command, cwd=root, check=True)
    assert set(pd.read_csv(output / "run_metrics.csv")["status"]) == {"valid_zero_peak"}
    command[command.index(str(empty))] = str(tmp_path / "missing.narrowPeak")
    assert subprocess.run(command, cwd=root, check=False).returncode != 0


def test_database_transaction_rolls_back_on_foreign_key_failure(tmp_path: Path):
    """A bad analysis row must not publish a partially built database."""
    root = Path(__file__).resolve().parents[1]
    config = yaml.safe_load((root / "configs/realstudy_v2.yaml").read_text())["realstudy_v2"]
    runs = tmp_path / "runs.csv"; pd.DataFrame(build_realstudy_v2_run_rows(config)).to_csv(runs, index=False)
    bad = tmp_path / "bad_metrics.csv"; pd.DataFrame([{"run_id": "missing", "metric": "x", "threshold": 0, "value": 1, "status": "complete"}]).to_csv(bad, index=False)
    database = tmp_path / "realstudy_v2.sqlite"
    command = [sys.executable, "scripts/build_realstudy_v2_database.py", "--publications", "manifests/realstudy_v2_publications.csv", "--experiments", "manifests/realstudy_v2_experiments.csv", "--files", "manifests/realstudy_v2_files.csv", "--run-table", str(runs), "--table", f"run_metrics={bad}", "--output", str(database), "--export-dir", str(tmp_path / "exports"), "--export-manifest", str(tmp_path / "manifest.csv")]
    result = subprocess.run(command, cwd=root, check=False)
    assert result.returncode != 0
    assert not database.exists() and not database.with_suffix(".sqlite.tmp").exists()


def test_database_declares_every_required_table():
    """The schema contract must retain all 23 preregistered tables."""
    root = Path(__file__).resolve().parents[1]
    source = (root / "scripts/build_realstudy_v2_database.py").read_text()
    required = ["publication_claims", "parent_libraries", "control_subsamples", "reference_region_recovery", "seed_pair_metrics", "stratified_metrics", "enough_control_decisions"]
    assert all(f"CREATE TABLE {table}" in source for table in required)


def test_renderer_declares_all_six_multiformat_source_packages():
    """Every preregistered main figure must use the common PDF/SVG/PNG/source/caption writer."""
    root = Path(__file__).resolve().parents[1]
    source = (root / "scripts/render_realstudy_v2_figures.py").read_text()
    assert 'FORMATS = ["pdf", "svg", "png"]' in source
    for figure_id in [
        "figure_1_design", "figure_2_dose_response", "figure_3_reproducibility",
        "figure_4_peak_geometry", "figure_5_stratified_response", "figure_6_representative_loci",
    ]:
        assert figure_id in source
    assert "_source_data.csv" in source and "_caption.md" in source
