"""Helpers for hilly-density tuning summaries."""

"""Imports"""

import pandas as pd


"""Functions"""


def build_recommendation_table(map_qc_df: pd.DataFrame) -> pd.DataFrame:
    """Rank candidate mappability settings by mildness, bumpiness, and FP burden."""
    ranked = map_qc_df.copy()
    had_bumpiness_delta = "bumpiness_delta_vs_wavy" in ranked.columns
    for column, default in {
        "p95_over_nonzero_median": 0.0,
        "p99_over_nonzero_median": 0.0,
        "max_over_nonzero_median": 0.0,
        "fp_burden_ratio_vs_wavy": 1.0,
        "bumpiness_delta_vs_wavy": 0.0,
    }.items():
        if column not in ranked.columns:
            ranked[column] = default

    ranked["rejected_for_fp_burden"] = ranked["fp_burden_ratio_vs_wavy"] > 1.5
    ranked["mildness_score"] = (
        ranked["p95_over_nonzero_median"]
        + 0.5 * ranked["p99_over_nonzero_median"]
        + 0.1 * ranked["max_over_nonzero_median"]
    )
    ranked["distinct_from_wavy"] = (
        ranked["bumpiness_delta_vs_wavy"] > 0 if had_bumpiness_delta else True
    )
    ranked["rank_score"] = (
        ranked["mildness_score"]
        + ranked["rejected_for_fp_burden"].astype(float) * 1000.0
        - ranked["distinct_from_wavy"].astype(float) * 0.25
    )
    ranked = ranked.sort_values(
        ["rejected_for_fp_burden", "distinct_from_wavy", "rank_score", "map_enrich"],
        ascending=[True, False, True, True],
    ).reset_index(drop=True)
    ranked["recommendation_rank"] = ranked.index + 1
    ranked["recommendation"] = "candidate"
    preferred_mask = (~ranked["rejected_for_fp_burden"]) & ranked["distinct_from_wavy"]
    if preferred_mask.any():
        preferred_index = ranked[preferred_mask].index[0]
        ranked.loc[preferred_index, "recommendation"] = "preferred_mild_candidate"
    elif not ranked.empty:
        ranked.loc[ranked.index[0], "recommendation"] = "fallback_candidate"
    return ranked
