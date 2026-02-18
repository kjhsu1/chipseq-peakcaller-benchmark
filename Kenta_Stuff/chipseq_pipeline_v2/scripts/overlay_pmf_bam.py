#!/usr/bin/env python3
"""Plot windowed overlays of expected (PMF-derived) vs observed (BAM depth) signal.

The script is intentionally standalone and lightweight:
- no Snakemake edits required
- works from existing run outputs
- supports condition filtering, peak index selection, chromosome override
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


@dataclass(frozen=True)
class Window:
    """0-based, half-open window."""

    name: str
    chrom: str
    start: int
    end: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Overlay PMF-derived expected coverage and BAM observed coverage in selected windows."
    )
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument("--params-csv", type=Path, default=None)
    parser.add_argument("--run-id", required=True, help="Run ID (e.g. 0037).")
    parser.add_argument(
        "--aligner",
        required=True,
        choices=("bowtie2", "bwa-mem"),
        help="Aligner folder name under results/{run_id}/",
    )
    parser.add_argument(
        "--window-bp",
        type=int,
        default=4000,
        help="Total window width in bp (default: 4000).",
    )
    parser.add_argument(
        "--cond",
        default="both",
        choices=("con", "treat", "both"),
        help="Which condition(s) to plot.",
    )
    parser.add_argument(
        "--peak-index",
        type=int,
        default=1,
        help="1-based planted peak index used for peak-centered window (default: 1).",
    )
    parser.add_argument("--chrom", default=None, help="Optional chromosome override.")
    parser.add_argument(
        "--center", type=int, default=None, help="Optional explicit center for peak window."
    )
    parser.add_argument(
        "--bg-center",
        type=int,
        default=None,
        help="Optional explicit center for background window.",
    )
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument(
        "--show-nb-band",
        action="store_true",
        help="Show approximate NB variance band around expected curve.",
    )
    return parser.parse_args()


def load_run_row(params_csv: Path, run_id: str) -> pd.Series:
    params = pd.read_csv(params_csv, dtype={"run_id": str})
    row = params.loc[params["run_id"] == run_id]
    if row.empty:
        raise ValueError(f"run_id={run_id} not found in {params_csv}")
    return row.iloc[0]


def load_planted_peaks(path: Path) -> List[Tuple[str, int]]:
    peaks: List[Tuple[str, int]] = []
    if not path.exists():
        return peaks
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                continue
            peaks.append((fields[0], int(fields[1])))
    return peaks


def choose_peak_center(
    peaks: Sequence[Tuple[str, int]],
    peak_index_1based: int,
    chrom_override: str | None,
    center_override: int | None,
) -> Tuple[str, int]:
    if center_override is not None:
        if chrom_override is None:
            raise ValueError("--center requires --chrom.")
        return chrom_override, int(center_override)

    if not peaks:
        raise ValueError(
            "No planted peaks found. Provide --chrom and --center explicitly for peak window."
        )

    if chrom_override is not None:
        peaks = [p for p in peaks if p[0] == chrom_override]
        if not peaks:
            raise ValueError(f"No planted peaks on chromosome {chrom_override}.")

    idx = peak_index_1based - 1
    if idx < 0 or idx >= len(peaks):
        raise ValueError(
            f"--peak-index {peak_index_1based} out of range for available peaks ({len(peaks)})."
        )
    return peaks[idx]


def make_window(chrom: str, center: int, width: int, chrom_len: int, name: str) -> Window:
    half = width // 2
    start = max(0, center - half)
    end = min(chrom_len, start + width)
    start = max(0, end - width)
    if end <= start:
        raise ValueError(f"Invalid window {chrom}:{start}-{end}")
    return Window(name=name, chrom=chrom, start=start, end=end)


def select_background_center(
    chrom_len: int,
    width: int,
    peak_centers: Iterable[int],
    forced_center: int | None,
) -> int:
    if forced_center is not None:
        return int(forced_center)

    peak_set = set(int(p) for p in peak_centers)
    candidate = chrom_len // 2
    step = max(1, width // 4)
    for direction in (1, -1):
        c = candidate
        while 0 <= c < chrom_len:
            if all(abs(c - p) > width for p in peak_set):
                return c
            c += direction * step
    return max(0, min(chrom_len - 1, candidate))


def load_pmf_curve(pmf_csv: Path, chrom: str) -> np.ndarray:
    df = pd.read_csv(pmf_csv)
    # Some archived PMFs store chrom labels like "I 1..150725"; normalize to first token.
    canon = df["chrom"].astype(str).str.split().str[0]
    df = df.loc[canon == chrom].sort_values("bin_idx")
    if df.empty:
        raise ValueError(f"No PMF entries for chromosome {chrom} in {pmf_csv}")
    return df["pmf"].to_numpy(dtype=float)


def expected_coverage_from_pmf(pmf: np.ndarray, frag_len: int, n_frag: float) -> np.ndarray:
    mu_start = np.maximum(0.0, pmf * n_frag)
    kernel = np.ones(int(frag_len), dtype=float)
    return np.convolve(mu_start, kernel, mode="full")


def expected_sd_from_nb(pmf: np.ndarray, frag_len: int, n_frag: float, nb_k: float) -> np.ndarray:
    mu_start = np.maximum(0.0, pmf * n_frag)
    # NB variance: Var = mu + mu^2 / k
    var_start = mu_start + (mu_start * mu_start) / max(nb_k, 1e-9)
    kernel = np.ones(int(frag_len), dtype=float)
    var_cov = np.convolve(var_start, kernel, mode="full")
    return np.sqrt(np.maximum(0.0, var_cov))


def extract_depth_window(bam: Path, chrom: str, start: int, end: int) -> np.ndarray:
    length = end - start
    depth = np.zeros(length, dtype=float)
    region = f"{chrom}:{start + 1}-{end}"
    cmd = ["samtools", "depth", "-aa", "-r", region, str(bam)]
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        c, pos1, d = line.split("\t")
        if c != chrom:
            continue
        idx = int(pos1) - 1 - start
        if 0 <= idx < length:
            depth[idx] = float(d)
    return depth


def scale_expected_to_observed(expected: np.ndarray, observed: np.ndarray) -> Tuple[np.ndarray, float]:
    esum = float(expected.sum())
    osum = float(observed.sum())
    if esum <= 0 or osum <= 0:
        return expected.copy(), 1.0
    factor = osum / esum
    return expected * factor, factor


def metrics(expected: np.ndarray, observed: np.ndarray) -> Tuple[float, float, float]:
    if expected.size == 0 or observed.size == 0:
        return float("nan"), float("nan"), float("nan")
    rmse = float(np.sqrt(np.mean((expected - observed) ** 2)))
    def corr_manual(a: np.ndarray, b: np.ndarray) -> float:
        a_centered = a - np.nanmean(a)
        b_centered = b - np.nanmean(b)
        denom = float(np.sqrt(np.nansum(a_centered**2) * np.nansum(b_centered**2)))
        if not np.isfinite(denom) or denom <= 1e-12:
            return float("nan")
        return float(np.nansum(a_centered * b_centered) / denom)

    pearson = corr_manual(expected, observed)
    exp_rank = pd.Series(expected).rank(method="average").to_numpy(dtype=float)
    obs_rank = pd.Series(observed).rank(method="average").to_numpy(dtype=float)
    spearman = corr_manual(exp_rank, obs_rank)
    return pearson, spearman, rmse


def conds_from_arg(cond: str) -> List[str]:
    if cond == "both":
        return ["con", "treat"]
    return [cond]


def plot_condition(
    run_id: str,
    cond: str,
    windows: List[Window],
    exp_curve: np.ndarray,
    exp_sd_curve: np.ndarray | None,
    bam_path: Path,
    out_png: Path,
) -> List[dict]:
    rows: List[dict] = []
    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=False)
    if len(windows) == 1:
        axes = [axes]

    for ax, w in zip(axes, windows):
        obs = extract_depth_window(bam_path, w.chrom, w.start, w.end)
        exp = exp_curve[w.start : w.end]
        exp_scaled, factor = scale_expected_to_observed(exp, obs)

        x = np.arange(w.start, w.end)
        ax.plot(x, obs, label="observed (BAM depth)", linewidth=1.2)
        ax.plot(x, exp_scaled, label="expected (PMF-derived)", linewidth=1.2)

        if exp_sd_curve is not None:
            sd = exp_sd_curve[w.start : w.end] * factor
            lo = np.maximum(0.0, exp_scaled - 1.96 * sd)
            hi = exp_scaled + 1.96 * sd
            ax.fill_between(x, lo, hi, alpha=0.2, label="expected +/- 1.96 SD")

        pearson, spearman, rmse = metrics(exp_scaled, obs)
        ax.set_title(
            f"{w.name}: {w.chrom}:{w.start}-{w.end} | r={pearson:.3f}, rho={spearman:.3f}, rmse={rmse:.3f}"
        )
        ax.set_xlabel("genomic position (bp)")
        ax.set_ylabel("coverage")
        ax.legend(loc="upper right", fontsize=8)

        rows.append(
            {
                "run_id": run_id,
                "cond": cond,
                "window_type": w.name,
                "chrom": w.chrom,
                "start": w.start,
                "end": w.end,
                "pearson_r": pearson,
                "spearman_r": spearman,
                "rmse": rmse,
            }
        )

    fig.suptitle(f"Run {run_id} | Condition: {cond}")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180)
    plt.close(fig)
    return rows


def main() -> None:
    args = parse_args()
    params_csv = args.params_csv or (args.results_dir / "params" / "run_params.csv")
    row = load_run_row(params_csv, args.run_id)

    frag_len = int(row["fragment_length"])
    nb_k = float(row["nb_k"])
    cov_treat = float(row["coverage_treat"])
    cov_ctrl = float(row["coverage_ctrl"])

    peaks_bed = args.results_dir / args.run_id / "treat" / "planted_peaks.bed"
    planted = load_planted_peaks(peaks_bed)
    peak_chrom, peak_center = choose_peak_center(
        peaks=planted,
        peak_index_1based=args.peak_index,
        chrom_override=args.chrom,
        center_override=args.center,
    )

    pmf_csv_treat = args.results_dir / args.run_id / "treat" / "pmf.csv"
    pmf_treat = load_pmf_curve(pmf_csv_treat, peak_chrom)
    n_frag_treat = cov_treat * len(pmf_treat) / max(frag_len, 1)
    exp_cov = expected_coverage_from_pmf(pmf_treat, frag_len, n_frag_treat)
    chrom_len = exp_cov.size

    peak_window = make_window(
        chrom=peak_chrom,
        center=peak_center,
        width=args.window_bp,
        chrom_len=chrom_len,
        name="peak",
    )

    peak_centers_same_chrom = [p for c, p in planted if c == peak_chrom]
    bg_center = select_background_center(
        chrom_len=chrom_len,
        width=args.window_bp,
        peak_centers=peak_centers_same_chrom,
        forced_center=args.bg_center,
    )
    bg_window = make_window(
        chrom=peak_chrom,
        center=bg_center,
        width=args.window_bp,
        chrom_len=chrom_len,
        name="background",
    )
    windows = [peak_window, bg_window]

    out_dir = args.out_dir or (args.results_dir / args.run_id / "viz")
    out_dir.mkdir(parents=True, exist_ok=True)

    all_rows: List[dict] = []
    for cond in conds_from_arg(args.cond):
        pmf_csv = args.results_dir / args.run_id / cond / "pmf.csv"
        if cond == "treat":
            n_frag = cov_treat * chrom_len / max(frag_len, 1)
        else:
            n_frag = cov_ctrl * chrom_len / max(frag_len, 1)
        pmf = load_pmf_curve(pmf_csv, peak_chrom)
        exp_curve = expected_coverage_from_pmf(pmf, frag_len, n_frag)
        exp_sd = (
            expected_sd_from_nb(pmf, frag_len, n_frag, nb_k)
            if args.show_nb_band
            else None
        )

        bam = (
            args.results_dir
            / args.run_id
            / args.aligner
            / cond
            / "aligned.sorted.bam"
        )
        if not bam.exists():
            raise FileNotFoundError(f"BAM not found: {bam}")

        out_png = out_dir / f"overlay_{args.run_id}_{cond}.png"
        rows = plot_condition(
            run_id=args.run_id,
            cond=cond,
            windows=windows,
            exp_curve=exp_curve,
            exp_sd_curve=exp_sd,
            bam_path=bam,
            out_png=out_png,
        )
        all_rows.extend(rows)
        print(f"Wrote {out_png}")

    metrics_csv = out_dir / f"overlay_metrics_{args.run_id}.csv"
    pd.DataFrame(all_rows).to_csv(metrics_csv, index=False)
    print(f"Wrote {metrics_csv}")


if __name__ == "__main__":
    main()
