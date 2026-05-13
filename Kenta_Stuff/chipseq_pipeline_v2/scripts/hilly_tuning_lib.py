"""Helpers for hilly-density tuning summaries."""

"""Imports"""

import pandas as pd


"""Functions"""


def build_recommendation_table(map_qc_df: pd.DataFrame) -> pd.DataFrame:
    """Rank candidate mappability settings by mildness and bumpiness."""
    ranked = map_qc_df.copy()
    ranked["rank_score"] = (
        ranked["p95_over_nonzero_median"]
        + 0.5 * ranked["p99_over_nonzero_median"]
        + 0.1 * ranked["max_over_nonzero_median"]
    )
    ranked["recommendation_rank"] = ranked["rank_score"].rank(method="dense")
    ranked["recommendation"] = ranked["map_enrich"].apply(
        lambda value: "preferred_mild_candidate" if float(value) == float(ranked["map_enrich"].min()) else "candidate"
    )
    return ranked.sort_values(["recommendation_rank", "map_enrich"])
