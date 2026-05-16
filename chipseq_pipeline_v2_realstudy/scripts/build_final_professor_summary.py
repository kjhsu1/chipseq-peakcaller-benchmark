"""Assemble a professor-facing final summary directory."""

"""Imports"""

import argparse
from datetime import datetime
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Assemble final benchmark summary notes.")
    parser.add_argument(
        "--controlled-root",
        type=Path,
        default=Path("../chipseq_pipeline_v2/analysis_outputs"),
    )
    parser.add_argument(
        "--realstudy-root",
        type=Path,
        default=Path("analysis_outputs"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis_outputs") / f"final_professor_summary_{datetime.now().strftime('%Y%m%d')}",
    )
    return parser.parse_args()


def main() -> None:
    """Write a concise final summary directory index."""
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Final Professor-Facing Summary",
        "",
        "## Controlled Truth-Aware Conclusions",
        f"- controlled analysis root: `{args.controlled_root}`",
        "- use the TF-clean six-category outputs as the controlled baseline",
        "",
        "## Realstudy Truth-Proxy Conclusions",
        f"- realstudy analysis root: `{args.realstudy_root}`",
        "- use estimated/reference/truth-proxy language only",
        "",
        "## Remaining Scientific Limitations",
        "- broad-study live metadata resolution is still required",
        "- ontology outputs should be interpreted as prototype truth-proxy summaries",
        "",
        "## Next Steps",
        "- finalize live study manifests",
        "- run smoke-scale realstudy execution before the 144-run prototype",
    ]
    (args.output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


main()
