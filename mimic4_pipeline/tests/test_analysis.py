"""Unit tests for the guarded statistical comparison in task 7 -- the
brief's explicit requirement that this "runs end-to-end on the demo even
when group sizes are zero or one ... degrade gracefully with a clear
message, not crash or emit a spurious p-value."

Effect size (mean difference + 95% CI) is the primary result as of the
2026-08-23 revision; the Mann-Whitney p-value is secondary. Every "no
spurious p-value" guarantee below applies equally to the effect-size CI --
neither is ever populated on a status other than "ok".
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from mimic4_pipeline.analysis import (
    CLINICAL_THRESHOLDS,
    MIN_GROUP_SIZE_FOR_TEST,
    PRIMARY_ANALYTES,
    adjusted_odds_ratio,
    safe_two_group_comparison,
)


def _assert_nothing_populated(result):
    assert result["mean_diff_a_minus_b"] is None
    assert result["mean_diff_ci_low"] is None
    assert result["mean_diff_ci_high"] is None
    assert result["p_value_mannwhitney_secondary"] is None


def test_both_groups_empty_does_not_crash_and_reports_nothing():
    result = safe_two_group_comparison([], [], "test")
    _assert_nothing_populated(result)
    assert result["status"].startswith("insufficient_data")


def test_one_group_empty_does_not_crash_and_reports_nothing():
    result = safe_two_group_comparison([], [1.0, 2.0, 3.0], "test")
    _assert_nothing_populated(result)
    assert "group a is empty" in result["status"]

    result2 = safe_two_group_comparison([1.0, 2.0, 3.0], [], "test")
    _assert_nothing_populated(result2)
    assert "group b is empty" in result2["status"]


def test_single_observation_per_group_does_not_crash_and_reports_nothing():
    assert MIN_GROUP_SIZE_FOR_TEST >= 2  # sanity: the threshold this test relies on
    result = safe_two_group_comparison([1.0], [2.0], "test")
    _assert_nothing_populated(result)
    assert result["status"].startswith("insufficient_data")


def test_nan_and_none_values_are_dropped_not_crashed_on():
    result = safe_two_group_comparison(
        [1.0, float("nan"), 2.0, None, 3.0], [4.0, 5.0, None, 6.0], "test"
    )
    assert result["n_group_a"] == 3
    assert result["n_group_b"] == 3
    assert result["status"] == "ok"
    assert result["mean_diff_a_minus_b"] is not None
    assert result["p_value_mannwhitney_secondary"] is not None


def test_degenerate_identical_constant_groups_does_not_crash():
    result = safe_two_group_comparison([5.0, 5.0, 5.0], [5.0, 5.0, 5.0], "test")
    assert result["status"] == "degenerate: both groups constant and equal -- no variation to test"
    _assert_nothing_populated(result)


def test_sufficient_data_produces_a_real_effect_size_with_ci_and_a_secondary_p_value():
    result = safe_two_group_comparison([1.0, 2.0, 3.0, 4.0], [10.0, 11.0, 12.0, 13.0], "test")
    assert result["status"] == "ok"
    assert result["mean_diff_a_minus_b"] == -9.0
    assert result["mean_diff_ci_low"] < result["mean_diff_a_minus_b"] < result["mean_diff_ci_high"]
    p = result["p_value_mannwhitney_secondary"]
    assert p is not None and 0.0 <= p <= 1.0 and not math.isnan(p)


def test_result_always_reports_group_sizes_even_when_skipped():
    result = safe_two_group_comparison([], [1.0], "test")
    assert result["n_group_a"] == 0
    assert result["n_group_b"] == 1


# --- pre-specified clinical thresholds (must exist for every primary analyte,
# and must be genuinely pre-specified -- i.e. present regardless of what any
# comparison run finds) ---------------------------------------------------


def test_clinical_thresholds_defined_for_every_primary_analyte():
    assert set(CLINICAL_THRESHOLDS) == set(PRIMARY_ANALYTES)


def test_hepatic_analytes_use_multiple_of_uln_metric():
    for analyte in ["ALT", "AST", "Total Bilirubin"]:
        assert CLINICAL_THRESHOLDS[analyte]["metric"] == "extreme_x_uln"


def test_sodium_uses_absolute_value_metric_not_uln():
    # No "multiple of ULN" convention applies to sodium -- it's assessed on
    # the absolute mEq/L scale (moderate hyponatremia <=130 mEq/L).
    assert CLINICAL_THRESHOLDS["Sodium"]["metric"] == "extreme_value"
    assert CLINICAL_THRESHOLDS["Sodium"]["op"] == "<="


# --- adjusted_odds_ratio (severity-adjusted comparison core) ------------
# DB-independent by design -- see the module docstring on why the
# clinical-threshold-crossing indicator, not the raw delta, is what gets
# adjusted -- so these run on plain synthetic DataFrames.


def _synthetic_logit_df(n_treated=200, n_untreated=200, treated_rate=0.5, untreated_rate=0.1, seed=0):
    rng = np.random.default_rng(seed)
    treated = pd.DataFrame(
        {
            "is_tricyclic": 1,
            "crosses": rng.binomial(1, treated_rate, n_treated),
            "charlson_comorbidity_index": rng.integers(0, 6, n_treated),
        }
    )
    untreated = pd.DataFrame(
        {
            "is_tricyclic": 0,
            "crosses": rng.binomial(1, untreated_rate, n_untreated),
            "charlson_comorbidity_index": rng.integers(0, 6, n_untreated),
        }
    )
    return pd.concat([treated, untreated], ignore_index=True)


def test_adjusted_odds_ratio_insufficient_data_does_not_crash():
    df = pd.DataFrame({"crosses": [1], "is_tricyclic": [1], "charlson_comorbidity_index": [2]})
    result = adjusted_odds_ratio(df, "crosses", "is_tricyclic", ["charlson_comorbidity_index"])
    assert result["status"].startswith("insufficient_data")
    assert result["adjusted_or"] is None


def test_adjusted_odds_ratio_degenerate_no_outcome_variation_does_not_crash():
    df = pd.DataFrame(
        {
            "crosses": [0] * 10,
            "is_tricyclic": [1, 0] * 5,
            "charlson_comorbidity_index": [1, 2] * 5,
        }
    )
    result = adjusted_odds_ratio(df, "crosses", "is_tricyclic", ["charlson_comorbidity_index"])
    assert result["status"] == "degenerate: no variation in the outcome (all events or none)"
    assert result["adjusted_or"] is None


def test_adjusted_odds_ratio_separation_within_one_group_does_not_crash():
    # Every event falls in the treated group; untreated has zero events --
    # a real (quasi-)separation scenario a naive fit would mishandle.
    df = pd.DataFrame(
        {
            "crosses": [1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            "is_tricyclic": [1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
            "charlson_comorbidity_index": [1, 2, 3, 1, 2, 1, 2, 3, 1, 2],
        }
    )
    result = adjusted_odds_ratio(df, "crosses", "is_tricyclic", ["charlson_comorbidity_index"])
    assert result["status"].startswith("degenerate: no outcome variation within")
    assert result["adjusted_or"] is None


def test_adjusted_odds_ratio_valid_fit_reports_or_with_sane_ci():
    df = _synthetic_logit_df(treated_rate=0.5, untreated_rate=0.1)
    result = adjusted_odds_ratio(df, "crosses", "is_tricyclic", ["charlson_comorbidity_index"])
    assert result["status"] == "ok"
    assert result["adjusted_or"] > 1.0  # treated group has the higher event rate by construction
    assert result["adjusted_or_ci_low"] < result["adjusted_or"] < result["adjusted_or_ci_high"]
    assert result["adjusted_p_value"] is not None and 0.0 <= result["adjusted_p_value"] <= 1.0
    assert result["n_treated"] == 200 and result["n_untreated"] == 200


def test_adjusted_odds_ratio_equal_rates_gives_or_near_one():
    df = _synthetic_logit_df(treated_rate=0.2, untreated_rate=0.2, seed=1)
    result = adjusted_odds_ratio(df, "crosses", "is_tricyclic", ["charlson_comorbidity_index"])
    assert result["status"] == "ok"
    assert 0.5 < result["adjusted_or"] < 2.0  # no true effect -- should land near 1, not blow up
