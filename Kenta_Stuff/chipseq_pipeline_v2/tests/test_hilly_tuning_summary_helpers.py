"""Tests for hilly tuning summary helpers."""

"""Imports"""

import pandas as pd

from scripts.hilly_tuning_lib import build_recommendation_table


"""Functions"""


def test_build_recommendation_table_prefers_mildest_candidate():
    """The mildest tested enrichment should be the default preferred candidate."""
    df = pd.DataFrame(
        [
            {"map_enrich": 10, "p95_over_nonzero_median": 2.0, "p99_over_nonzero_median": 2.6, "max_over_nonzero_median": 5.2},
            {"map_enrich": 50, "p95_over_nonzero_median": 2.2, "p99_over_nonzero_median": 3.0, "max_over_nonzero_median": 7.4},
        ]
    )
    ranked = build_recommendation_table(df)
    assert ranked.iloc[0]["map_enrich"] == 10
    assert ranked.iloc[0]["recommendation"] == "preferred_mild_candidate"
