"""
Summarize false-negative and false-positive patterns from peak recovery outputs.
"""

"""Imports"""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


"""Functions"""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Summarize seed, parameter, and locus patterns in peak recovery FP/FN outputs."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing per_run_recovery_summary.csv, missed_planted_peaks.csv, and false_positive_called_peaks.csv.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for summary outputs. Defaults to --input-dir.",
    )
    parser.add_argument(
        "--top-loci",
        type=int,
        default=25,
        help="Number of top repeated loci to include per category.",
    )
    return parser.parse_args()


def read_csv(path: Path) -> List[Dict[str, str]]:
    """Read a CSV into dictionaries."""
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: List[Dict[str, str]], fieldnames: List[str]) -> None:
    """Write dictionaries to a CSV."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def numeric_sort_key(value: str) -> Tuple[int, float, str]:
    """Sort numeric strings numerically and other strings lexically."""
    try:
        return (0, float(value), value)
    except ValueError:
        return (1, 0.0, value)


def distribution(values: Iterable[int]) -> str:
    """Return a compact count distribution string."""
    counts = Counter(values)
    return "; ".join(f"{value}:{count}" for value, count in sorted(counts.items()))


def summarize_recall_by_seed(per_run: List[Dict[str, str]]) -> List[Dict[str, str]]:
    """Summarize perfect recall and missed planted counts by category and seed."""
    grouped: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    for row in per_run:
        grouped[(row["category"], row["seed"])].append(row)

    rows: List[Dict[str, str]] = []
    for (category, seed), group in sorted(grouped.items(), key=lambda item: (item[0][0], numeric_sort_key(item[0][1]))):
        missed_counts = [int(float(row["missed_planted"])) for row in group]
        rows.append(
            {
                "category": category,
                "seed": seed,
                "runs": str(len(group)),
                "perfect_recall_runs": str(sum(count == 0 for count in missed_counts)),
                "imperfect_recall_runs": str(sum(count > 0 for count in missed_counts)),
                "total_missed_planted": str(sum(missed_counts)),
                "missed_per_run_distribution": distribution(missed_counts),
            }
        )
    return rows


def summarize_precision_by_seed(per_run: List[Dict[str, str]]) -> List[Dict[str, str]]:
    """Summarize perfect precision and false-positive counts by category and seed."""
    grouped: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    for row in per_run:
        grouped[(row["category"], row["seed"])].append(row)

    rows: List[Dict[str, str]] = []
    for (category, seed), group in sorted(grouped.items(), key=lambda item: (item[0][0], numeric_sort_key(item[0][1]))):
        fp_counts = [int(float(row["false_positive_called"])) for row in group]
        rows.append(
            {
                "category": category,
                "seed": seed,
                "runs": str(len(group)),
                "perfect_precision_runs": str(sum(count == 0 for count in fp_counts)),
                "imperfect_precision_runs": str(sum(count > 0 for count in fp_counts)),
                "total_false_positive_called": str(sum(fp_counts)),
                "false_positive_per_run_distribution": distribution(fp_counts),
            }
        )
    return rows


def summarize_missed_loci(missed: List[Dict[str, str]], top_loci: int) -> List[Dict[str, str]]:
    """Summarize repeated missed planted centers."""
    counts = Counter(
        (
            row["category"],
            row["seed"],
            row["planted_chrom"],
            row["planted_center"],
        )
        for row in missed
    )
    category_counts: Dict[str, int] = defaultdict(int)
    rows: List[Dict[str, str]] = []
    for (category, seed, chrom, center), count in counts.most_common():
        if category_counts[category] >= top_loci:
            continue
        category_counts[category] += 1
        rows.append(
            {
                "category": category,
                "seed": seed,
                "planted_chrom": chrom,
                "planted_center": center,
                "missed_count": str(count),
            }
        )
    return rows


def summarize_fp_loci(false_positives: List[Dict[str, str]], top_loci: int) -> List[Dict[str, str]]:
    """Summarize repeated false-positive called loci using 100 bp midpoint bins."""
    exact_counts = Counter(
        (
            row["category"],
            row["called_chrom"],
            row["called_start"],
            row["called_end"],
        )
        for row in false_positives
    )
    bin_counts: Counter = Counter()
    examples: Dict[Tuple[str, str, str], Tuple[str, str]] = {}
    for row in false_positives:
        start = int(float(row["called_start"]))
        end = int(float(row["called_end"]))
        midpoint_bin = str(((start + end) // 2 // 100) * 100)
        key = (row["category"], row["called_chrom"], midpoint_bin)
        bin_counts[key] += 1
        exact_key = (row["category"], row["called_chrom"], row["called_start"], row["called_end"])
        if key not in examples or exact_counts[exact_key] > exact_counts[(row["category"], row["called_chrom"], examples[key][0], examples[key][1])]:
            examples[key] = (row["called_start"], row["called_end"])

    category_counts: Dict[str, int] = defaultdict(int)
    rows: List[Dict[str, str]] = []
    for (category, chrom, midpoint_bin), count in bin_counts.most_common():
        if category_counts[category] >= top_loci:
            continue
        category_counts[category] += 1
        example_start, example_end = examples[(category, chrom, midpoint_bin)]
        rows.append(
            {
                "category": category,
                "called_chrom": chrom,
                "called_midpoint_100bp_bin": midpoint_bin,
                "false_positive_count": str(count),
                "representative_called_start": example_start,
                "representative_called_end": example_end,
            }
        )
    return rows


def summarize_parameter_grid(per_run: List[Dict[str, str]], metric: str) -> List[Dict[str, str]]:
    """Summarize perfect metric counts by swept parameter."""
    flag = f"has_perfect_{metric}"
    rows: List[Dict[str, str]] = []
    for category in sorted({row["category"] for row in per_run}):
        category_rows = [row for row in per_run if row["category"] == category]
        for parameter in ["seed", "coverage_treat", "coverage_ctrl", "tf_enrich", "tf_sigma", "tf_peak_count_treat"]:
            grouped: Dict[str, List[Dict[str, str]]] = defaultdict(list)
            for row in category_rows:
                grouped[row[parameter]].append(row)
            for value, group in sorted(grouped.items(), key=lambda item: numeric_sort_key(item[0])):
                rows.append(
                    {
                        "category": category,
                        "metric": metric,
                        "parameter": parameter,
                        "value": value,
                        "runs": str(len(group)),
                        "perfect_runs": str(sum(row[flag] == "True" for row in group)),
                    }
                )
    return rows


def render_markdown(
    per_run: List[Dict[str, str]],
    recall_seed: List[Dict[str, str]],
    precision_seed: List[Dict[str, str]],
    missed_loci: List[Dict[str, str]],
    fp_loci: List[Dict[str, str]],
) -> str:
    """Render a concise markdown interpretation."""
    categories = sorted({row["category"] for row in per_run})
    lines = [
        "# FP/FN Pattern Summary",
        "",
        "Matching criterion: a called peak is a true positive when it contains the planted 1-bp center. No padding is used.",
        "",
        "## Recall",
        "",
        "The 96/288 perfect-recall pattern is seed-driven, not coverage-driven. In each affected category, exactly two of six seeds have 48/48 perfect-recall runs; the other four seeds have at least one missed planted peak in every swept coverage/enrichment combination.",
        "",
        "| category | perfect-recall seeds | imperfect-recall seeds |",
        "| --- | --- | --- |",
    ]
    for category in categories:
        category_rows = [row for row in recall_seed if row["category"] == category]
        perfect = [row["seed"] for row in category_rows if row["perfect_recall_runs"] == row["runs"]]
        imperfect = [row["seed"] for row in category_rows if row["imperfect_recall_runs"] != "0"]
        lines.append(f"| {category} | {', '.join(perfect) or 'none'} | {', '.join(imperfect) or 'none'} |")

    lines.extend(
        [
            "",
            "The repeated missed centers confirm this is mostly a planted-location effect. A center missed 48 times means the same planted location failed across all 48 parameter combinations for that seed.",
            "",
            "## Precision",
            "",
            "Perfect precision is rare because almost every seed repeatedly produces at least one extra called interval. The only common perfect-precision seed is seed 37 in the two narrow/peak categories, but those runs have poor recall, so they are precise mostly because MACS2 calls very few peaks there.",
            "",
            "| category | perfect-precision runs | perfect-precision seeds/conditions |",
            "| --- | ---: | --- |",
        ]
    )
    for category in categories:
        category_rows = [row for row in per_run if row["category"] == category]
        perfect_rows = [row for row in category_rows if row["has_perfect_precision"] == "True"]
        if not perfect_rows:
            condition = "none"
        else:
            seeds = sorted({row["seed"] for row in perfect_rows}, key=numeric_sort_key)
            coverage_treat = sorted({row["coverage_treat"] for row in perfect_rows}, key=numeric_sort_key)
            tf_enrich = sorted({row["tf_enrich"] for row in perfect_rows}, key=numeric_sort_key)
            condition = f"seeds={','.join(seeds)}; treat_cov={','.join(coverage_treat)}; tf_enrich={','.join(tf_enrich)}"
        lines.append(f"| {category} | {len(perfect_rows)} | {condition} |")

    lines.extend(
        [
            "",
            "False positives also recur at specific genomic hotspots. A 100 bp bin with count 48 means the same approximate locus was falsely called across every parameter combination for one seed.",
            "",
            "## Output Tables",
            "",
            "- `recall_seed_summary.csv`: missed planted counts by category and seed.",
            "- `precision_seed_summary.csv`: false-positive called counts by category and seed.",
            "- `missed_locus_summary.csv`: top repeated missed planted centers.",
            "- `false_positive_locus_summary.csv`: top repeated false-positive called loci binned by called-interval midpoint.",
            "- `parameter_perfect_metric_summary.csv`: perfect recall/precision counts by swept parameter.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    """Run the summary workflow."""
    args = parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir or input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    per_run = read_csv(input_dir / "per_run_recovery_summary.csv")
    missed = read_csv(input_dir / "missed_planted_peaks.csv")
    false_positives = read_csv(input_dir / "false_positive_called_peaks.csv")

    recall_seed = summarize_recall_by_seed(per_run)
    precision_seed = summarize_precision_by_seed(per_run)
    missed_loci = summarize_missed_loci(missed, args.top_loci)
    fp_loci = summarize_fp_loci(false_positives, args.top_loci)
    parameter_grid = summarize_parameter_grid(per_run, "recall") + summarize_parameter_grid(per_run, "precision")

    write_csv(
        output_dir / "recall_seed_summary.csv",
        recall_seed,
        ["category", "seed", "runs", "perfect_recall_runs", "imperfect_recall_runs", "total_missed_planted", "missed_per_run_distribution"],
    )
    write_csv(
        output_dir / "precision_seed_summary.csv",
        precision_seed,
        ["category", "seed", "runs", "perfect_precision_runs", "imperfect_precision_runs", "total_false_positive_called", "false_positive_per_run_distribution"],
    )
    write_csv(
        output_dir / "missed_locus_summary.csv",
        missed_loci,
        ["category", "seed", "planted_chrom", "planted_center", "missed_count"],
    )
    write_csv(
        output_dir / "false_positive_locus_summary.csv",
        fp_loci,
        ["category", "called_chrom", "called_midpoint_100bp_bin", "false_positive_count", "representative_called_start", "representative_called_end"],
    )
    write_csv(
        output_dir / "parameter_perfect_metric_summary.csv",
        parameter_grid,
        ["category", "metric", "parameter", "value", "runs", "perfect_runs"],
    )
    (output_dir / "fp_fn_pattern_summary.md").write_text(
        render_markdown(per_run, recall_seed, precision_seed, missed_loci, fp_loci),
        encoding="utf-8",
    )


main()
