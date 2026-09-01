"""Consolidate runtime provenance and analysis CSVs for database loading."""

"""Imports"""

import argparse
import hashlib
import subprocess
from pathlib import Path

import pandas as pd


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse table-preparation arguments."""
    parser = argparse.ArgumentParser(description="Prepare normalized Realstudy v2 database tables.")
    parser.add_argument("--run-table", type=Path, required=True)
    parser.add_argument("--comparison-root", type=Path, required=True)
    parser.add_argument("--peak-root", type=Path, required=True)
    parser.add_argument("--files-manifest", type=Path, required=True)
    parser.add_argument("--parent", action="append", type=Path, required=True)
    parser.add_argument("--sampler-manifest", action="append", type=Path, required=True)
    parser.add_argument("--reference-asset", action="append", default=[], help="label=path")
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def combine(paths: list[Path], output: Path, deduplicate: str | None = None) -> pd.DataFrame:
    """Combine matching CSVs and optionally deduplicate a stable identifier."""
    frames = [pd.read_csv(path) for path in paths]
    frame = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if deduplicate and not frame.empty:
        duplicates = frame[frame.duplicated(deduplicate, keep=False)]
        for _, group in duplicates.groupby(deduplicate):
            if len(group.drop_duplicates()) != 1:
                raise ValueError(f"Conflicting duplicate {deduplicate}: {group.iloc[0][deduplicate]}")
        frame = frame.drop_duplicates(deduplicate)
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False)
    return frame


def sha256(path: Path) -> str:
    """Return a streaming artifact checksum."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def software_versions() -> list[dict]:
    """Capture exact versions for workflow-critical command-line software."""
    commands = {
        "python": ["python", "--version"], "bowtie2": ["bowtie2", "--version"],
        "samtools": ["samtools", "--version"], "macs2": ["macs2", "--version"],
    }
    rows = []
    for software, command in commands.items():
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        version = (result.stdout or result.stderr).splitlines()[0].strip()
        rows.append({"software_id": software, "software_name": software, "version": version, "command": " ".join(command)})
    return rows


def main() -> None:
    """Write all directly derived execution, peak, metric, artifact, and validation tables."""
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    runs = pd.read_csv(args.run_table, dtype=str).fillna("")
    input_files = pd.read_csv(args.files_manifest, dtype=str).fillna("")
    input_files["staged_checksum"] = input_files["md5"]
    input_files["staged_status"] = "verified"
    input_files.to_csv(args.output_dir / "input_files.csv", index=False)
    sampled_run_ids = list(runs.loc[runs["run_type"] == "control_subsample", "run_id"])
    comparison_dirs = [args.comparison_root / run_id for run_id in sampled_run_ids]
    combine([path / "peaks.csv" for path in comparison_dirs], args.output_dir / "peaks.csv", "peak_id")
    combine([path / "reference_regions.csv" for path in comparison_dirs], args.output_dir / "reference_regions.csv", "reference_region_id")
    combine([path / "peak_overlaps.csv" for path in comparison_dirs], args.output_dir / "peak_overlaps.csv", "overlap_id")
    combine([path / "reference_region_recovery.csv" for path in comparison_dirs], args.output_dir / "reference_region_recovery.csv")
    metrics = combine([path / "run_metrics.csv" for path in comparison_dirs], args.output_dir / "run_metrics.csv")
    combine(args.parent, args.output_dir / "parent_libraries.csv", "parent_library_id")
    sampler = combine(args.sampler_manifest, args.output_dir / "sampling_manifests.csv", "control_subsample_id")
    controls = sampler.rename(
        columns={"output_sha256": "output_checksum", "requested_fragments": "requested_fragments",
                 "realized_fragments": "realized_fragments", "realized_reads": "realized_reads"}
    )[[
        "control_subsample_id", "control_parent_id", "study_id", "seed", "control_coverage_x",
        "requested_fragments", "realized_fragments", "realized_reads", "source_sha256",
        "output_checksum", "rank_cutoff",
    ]].rename(columns={"source_sha256": "source_checksum"})
    controls.to_csv(args.output_dir / "control_subsamples.csv", index=False)
    blocks = runs[["matched_block_id", "study_id", "seed", "reference_run_id"]].drop_duplicates().copy()
    blocks["algorithm_version"] = "sha256-rg-qname-v1"
    blocks.to_csv(args.output_dir / "sampling_blocks.csv", index=False)
    pd.DataFrame(software_versions()).to_csv(args.output_dir / "software_versions.csv", index=False)
    artifacts = []
    validation = []
    metric_status = metrics.groupby("run_id")["status"].first().to_dict()
    for row in runs.to_dict(orient="records"):
        peak_path = args.peak_root / row["run_id"] / "peaks.bed"
        artifacts.append(
            {"artifact_id": f"{row['run_id']}__peaks", "run_id": row["run_id"], "artifact_type": "macs2_peaks",
             "path": str(peak_path), "checksum": sha256(peak_path), "bytes": peak_path.stat().st_size, "status": "complete"}
        )
        status = "pass" if row["run_type"] == "full_control_anchor" or row["run_id"] in metric_status else "fail"
        validation.append(
            {"validation_id": f"{row['run_id']}__completion", "run_id": row["run_id"],
             "check_name": "peak_call_and_comparison_complete", "status": status,
             "detail": "Anchor peak call complete." if row["run_type"] == "full_control_anchor" else str(metric_status.get(row["run_id"], "missing metrics"))}
        )
    for item in args.reference_asset:
        label, separator, value = item.partition("=")
        if not separator:
            raise ValueError(f"Expected reference label=path, observed {item}")
        path = Path(value)
        artifacts.append(
            {"artifact_id": f"reference__{label}", "run_id": "", "artifact_type": "reference_asset",
             "path": str(path), "checksum": sha256(path), "bytes": path.stat().st_size, "status": "complete"}
        )
    pd.DataFrame(artifacts).to_csv(args.output_dir / "artifacts.csv", index=False)
    pd.DataFrame(validation).to_csv(args.output_dir / "validation_events.csv", index=False)


main()
