"""
Plot PMF by chromosome.

Usage:
  python -m scripts.plot_pmf --pmf_csv path/to/pmf.csv --out_dir path/to/plots
"""

import argparse
import os

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description="Plot PMF per chromosome")
    parser.add_argument("--pmf_csv", required=True, help="Path to pmf.csv with chrom column")
    parser.add_argument("--out_dir", required=True, help="Directory to write plots")
    parser.add_argument("--format", default="png", help="Image format (png, pdf, svg)")
    parser.add_argument("--dpi", type=int, default=200, help="Image DPI")
    return parser.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.pmf_csv)

    if "chrom" not in df.columns:
        raise ValueError("pmf.csv missing required 'chrom' column")

    os.makedirs(args.out_dir, exist_ok=True)

    for chrom, sub in df.groupby("chrom", sort=False):
        sub = sub.sort_values("bin_idx")
        plt.figure(figsize=(10, 3))
        plt.plot(sub["bin_idx"], sub["pmf"], linewidth=0.8)
        plt.title(f"PMF - {chrom}")
        plt.xlabel("bin_idx")
        plt.ylabel("pmf")
        plt.tight_layout()
        out_path = os.path.join(args.out_dir, f"pmf_{chrom}.{args.format}")
        plt.savefig(out_path, dpi=args.dpi)
        plt.close()


if __name__ == "__main__":
    main()
