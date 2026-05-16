"""End-to-end tests for seed standardization in the simulator."""

"""Imports"""

import csv
import subprocess
import sys
from pathlib import Path

from scripts.seed_helpers import derive_condition_seed


"""Functions"""


def write_fasta(path: Path) -> None:
    """Write a small multi-chromosome FASTA for simulator tests."""
    path.write_text(
        ">chr1\n"
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
        ">chr2\n"
        "TTGCAATGCCGATTAACCGGTTGCAATGCCGATTAACCGGTTGCAATGCCGA\n",
        encoding="utf-8",
    )


def read_lines(path: Path) -> list[str]:
    """Return stripped non-empty lines from a text file."""
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def read_pmf_rows(path: Path) -> list[dict[str, str]]:
    """Return PMF CSV rows."""
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def run_simulation(
    tmp_path: Path,
    run_name: str,
    extra_args: list[str],
) -> tuple[Path, Path]:
    """Run the simulator and return planted BED and PMF CSV paths."""
    repo_root = Path(__file__).resolve().parents[1]
    fasta = tmp_path / "toy.fa"
    if not fasta.exists():
        write_fasta(fasta)

    out_dir = tmp_path / run_name
    out_dir.mkdir()
    r1 = out_dir / "reads_R1.fasta"
    r2 = out_dir / "reads_R2.fasta"
    planted = out_dir / "planted_peaks.bed"
    planted_intervals = out_dir / "planted_peak_intervals.bed"
    pmf_csv = out_dir / "pmf.csv"

    defaults = {
        "--seed": "101",
        "--tf_peak_count": "3",
    }
    for idx in range(0, len(extra_args), 2):
        defaults[extra_args[idx]] = extra_args[idx + 1]

    cmd = [
        sys.executable,
        "-m",
        "scripts.updated_chip_seq",
        "--fasta", str(fasta),
        "--coverage", "0",
        "--fragment_length", "10",
        "--read_length", "5",
        "--tf_sigma", "1.5",
        "--tf_enrich", "4",
        "--nb_k", "1000",
        "--output_fasta1", str(r1),
        "--output_fasta2", str(r2),
        "--pmf_csv", str(pmf_csv),
        "--planted_peaks_bed", str(planted),
        "--planted_peak_intervals_bed", str(planted_intervals),
    ]
    for key, value in defaults.items():
        cmd.extend([key, value])
    subprocess.run(cmd, check=True, cwd=repo_root)
    return planted, pmf_csv


def test_seed_fallback_keeps_tf_centers_fixed_when_mappability_changes(tmp_path):
    """Shared seed should keep TF centers fixed even when mappability differs."""
    planted_a, _ = run_simulation(
        tmp_path,
        "seed_only_map0",
        ["--map_coverage_pct", "0"],
    )
    planted_b, _ = run_simulation(
        tmp_path,
        "seed_only_map20",
        ["--map_coverage_pct", "20", "--map_sigma", "1.5", "--map_enrich", "5"],
    )
    assert read_lines(planted_a) == read_lines(planted_b)


def test_map_seed_keeps_mappability_fixed_when_general_seed_changes(tmp_path):
    """Explicit map_seed should stabilize mappability PMFs across runs."""
    _, pmf_a = run_simulation(
        tmp_path,
        "map_seed_a",
        [
            "--seed", "101",
            "--tf_seed", "11",
            "--map_seed", "23",
            "--tf_peak_count", "0",
            "--map_coverage_pct", "20",
            "--map_sigma", "1.5",
            "--map_enrich", "5",
        ],
    )
    _, pmf_b = run_simulation(
        tmp_path,
        "map_seed_b",
        [
            "--seed", "999",
            "--tf_seed", "77",
            "--map_seed", "23",
            "--tf_peak_count", "0",
            "--map_coverage_pct", "20",
            "--map_sigma", "1.5",
            "--map_enrich", "5",
        ],
    )
    assert read_pmf_rows(pmf_a) == read_pmf_rows(pmf_b)


def test_condition_seed_derivation_is_distinct():
    """Condition-specific seed derivation should separate control and treatment RNGs."""
    base_seed = 101
    assert derive_condition_seed(base_seed, 31) != derive_condition_seed(base_seed, 32)
    assert derive_condition_seed(base_seed, 21) != derive_condition_seed(base_seed, 22)


def test_read_seed_changes_reads_without_changing_truth_or_pmf(tmp_path):
    """Read-seed changes should affect reads, not planted peaks or PMFs."""
    repo_root = Path(__file__).resolve().parents[1]
    fasta = tmp_path / "toy.fa"
    write_fasta(fasta)

    def run_named(name: str, read_seed: str) -> tuple[list[str], list[str], list[dict[str, str]]]:
        out_dir = tmp_path / name
        out_dir.mkdir()
        r1 = out_dir / "reads_R1.fasta"
        r2 = out_dir / "reads_R2.fasta"
        planted = out_dir / "planted_peaks.bed"
        planted_intervals = out_dir / "planted_peak_intervals.bed"
        pmf_csv = out_dir / "pmf.csv"
        cmd = [
            sys.executable,
            "-m",
            "scripts.updated_chip_seq",
            "--fasta", str(fasta),
            "--coverage", "1",
            "--fragment_length", "10",
            "--read_length", "5",
            "--tf_sigma", "1.5",
            "--tf_enrich", "4",
            "--nb_k", "1000",
            "--seed", "101",
            "--tf_seed", "11",
            "--map_seed", "23",
            "--read_seed", read_seed,
            "--map_coverage_pct", "20",
            "--map_sigma", "1.5",
            "--map_enrich", "5",
            "--output_fasta1", str(r1),
            "--output_fasta2", str(r2),
            "--pmf_csv", str(pmf_csv),
            "--planted_peaks_bed", str(planted),
            "--planted_peak_intervals_bed", str(planted_intervals),
        ]
        subprocess.run(cmd, check=True, cwd=repo_root)
        return read_lines(planted), read_lines(r1), read_pmf_rows(pmf_csv)

    planted_a, reads_a, pmf_a = run_named("read_seed_a", "31")
    planted_b, reads_b, pmf_b = run_named("read_seed_b", "32")
    assert planted_a == planted_b
    assert pmf_a == pmf_b
    assert reads_a != reads_b
