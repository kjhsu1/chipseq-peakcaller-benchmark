"""Re-render cleaned realistic-category stratified plots with a fixed layout."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir
    manifest = pd.read_csv(input_dir / "manifest.csv")

    colors = {5: "tab:blue", 10: "tab:orange", 15: "tab:green", 20: "tab:red"}
    metrics = [("precision", "Precision"), ("recall", "Recall"), ("f1", "F1")]

    for row in manifest.itertuples(index=False):
        csv_path = Path(row.output_csv)
        png_path = Path(row.output_png)
        pdf_path = Path(row.output_pdf)
        data = pd.read_csv(csv_path)

        fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharex=True)
        for ax, (metric, title) in zip(axes, metrics):
            for treat in sorted(data["coverage_treat"].unique()):
                line = data[data["coverage_treat"] == treat].sort_values("ratio")
                ax.plot(
                    line["ratio"],
                    line[metric],
                    marker="o",
                    linewidth=1.8,
                    color=colors.get(int(treat)),
                    label=f"treat={int(treat)}",
                )
            ax.set_xscale("log")
            ax.set_ylim(0, 1.05)
            ax.grid(True, alpha=0.25)
            ax.set_title(title)
            ax.set_xlabel("control:treat ratio")

        axes[0].set_ylabel("metric value")

        handles, labels = axes[0].get_legend_handles_labels()
        # Keep the legend centered under the suptitle and above the panels.
        fig.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.985),
            ncol=4,
            frameon=False,
            columnspacing=1.5,
            handletextpad=0.6,
        )

        subtitle = "all consistent families"
        if row.version == "drop_tfexp2_tf50":
            subtitle = "drop weak family (tf_exp=2, tf_enrich=50)"

        fig.suptitle(
            f"{row.category}\nStratified by fixed treatment coverage; {subtitle}",
            y=1.065,
            fontsize=13,
        )
        fig.tight_layout(rect=[0, 0, 1, 0.88])
        fig.savefig(png_path, dpi=220, bbox_inches="tight")
        fig.savefig(pdf_path, bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
