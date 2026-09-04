"""Unit tests for the propensity-score-matched sensitivity analysis --
DB-independent by design (see propensity.py's module docstring), so these
run on plain synthetic DataFrames.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from mimic4_pipeline.propensity import (
    fit_propensity_scores,
    match_on_propensity,
    standardized_mean_diff,
)


# --- standardized_mean_diff ----------------------------------------------


def test_smd_zero_when_groups_identical():
    df = pd.DataFrame({"is_tricyclic": [1, 1, 0, 0], "age": [50, 60, 50, 60]})
    result = standardized_mean_diff(df, "is_tricyclic", ["age"])
    assert result.loc[0, "smd"] == 0.0


def test_smd_positive_when_treated_group_higher():
    df = pd.DataFrame({"is_tricyclic": [1, 1, 1, 0, 0, 0], "age": [70, 72, 74, 40, 42, 44]})
    result = standardized_mean_diff(df, "is_tricyclic", ["age"])
    assert result.loc[0, "smd"] > 0
    assert result.loc[0, "mean_treated"] == 72.0
    assert result.loc[0, "mean_control"] == 42.0


def test_smd_does_not_crash_on_zero_variance_covariate():
    # Every admission has the same value (e.g. a comorbidity flag that
    # happens to be 0 for everyone in this slice) -- pooled SD is 0, and
    # dividing by it must not raise or return inf/nan.
    df = pd.DataFrame({"is_tricyclic": [1, 1, 0, 0], "aids": [0, 0, 0, 0]})
    result = standardized_mean_diff(df, "is_tricyclic", ["aids"])
    assert result.loc[0, "smd"] == 0.0


def test_smd_covers_every_requested_covariate():
    df = pd.DataFrame(
        {"is_tricyclic": [1, 0, 1, 0], "age": [70, 40, 72, 42], "renal_disease": [1, 0, 0, 1]}
    )
    result = standardized_mean_diff(df, "is_tricyclic", ["age", "renal_disease"])
    assert set(result["covariate"]) == {"age", "renal_disease"}


# --- fit_propensity_scores ------------------------------------------------


def test_fit_propensity_scores_returns_probabilities_in_unit_interval():
    rng = np.random.default_rng(0)
    n = 200
    age = rng.normal(60, 10, n)
    is_tricyclic = rng.binomial(1, 1 / (1 + np.exp(-(age - 60) / 10)))
    df = pd.DataFrame({"is_tricyclic": is_tricyclic, "age": age})
    result = fit_propensity_scores(df, "is_tricyclic", ["age"])
    assert (result["propensity_score"] >= 0).all()
    assert (result["propensity_score"] <= 1).all()
    assert len(result) == n


def test_fit_propensity_scores_drops_rows_with_missing_covariates():
    # Some age overlap between groups (not perfectly separable) so the tiny
    # synthetic fit converges cleanly rather than just emitting a warning.
    df = pd.DataFrame(
        {
            "is_tricyclic": [1, 0, 1, 0, 1, 0, 1, 0],
            "age": [70, 40, 72, 55, None, 60, 68, 41],
        }
    )
    result = fit_propensity_scores(df, "is_tricyclic", ["age"])
    assert len(result) == 7  # one row dropped for missing age
    assert result["age"].notna().all()


# --- match_on_propensity --------------------------------------------------


def _ps_df(ids, treatment, ps):
    return pd.DataFrame({"hadm_id": ids, "is_tricyclic": treatment, "propensity_score": ps})


def test_match_on_propensity_no_duplicate_control_usage():
    # 5 treated units all clustered near 3 available controls -- matching
    # without replacement must never reuse a control.
    df = _ps_df(
        ids=list(range(8)),
        treatment=[1, 1, 1, 1, 1, 0, 0, 0],
        ps=[0.50, 0.51, 0.52, 0.53, 0.54, 0.50, 0.51, 0.52],
    )
    matched = match_on_propensity(df, "is_tricyclic", "propensity_score", caliper=1.0)
    assert matched["control_hadm_id"].is_unique
    assert matched["treated_hadm_id"].is_unique
    assert len(matched) <= 3  # only 3 controls exist


def test_match_on_propensity_respects_caliper():
    # One treated unit's propensity is far from every control -- with a
    # tight caliper it must go unmatched rather than force-matched.
    df = _ps_df(
        ids=[1, 2, 3, 4],
        treatment=[1, 1, 0, 0],
        ps=[0.5, 0.001, 0.5, 0.5],  # unit 2 is far from both controls (logit scale)
    )
    matched = match_on_propensity(df, "is_tricyclic", "propensity_score", caliper=0.1)
    assert 2 not in matched["treated_hadm_id"].tolist()
    assert 1 in matched["treated_hadm_id"].tolist()


def test_match_on_propensity_all_close_treated_units_get_matched():
    # Abundant, overlapping controls -- the regime the function is actually
    # used in (~7k treated vs ~80k controls, an ~11x ratio here) -- should
    # match nearly everyone. A 1:1 ratio is a harder case for the
    # batch-collision approximation and is not what this asserts.
    rng = np.random.default_rng(0)
    n_treated, n_control = 100, 1000
    ps_treated = rng.uniform(0.4, 0.6, n_treated)
    ps_control = rng.uniform(0.4, 0.6, n_control)
    df = _ps_df(
        ids=list(range(n_treated + n_control)),
        treatment=[1] * n_treated + [0] * n_control,
        ps=np.concatenate([ps_treated, ps_control]),
    )
    matched = match_on_propensity(df, "is_tricyclic", "propensity_score", caliper=2.0)
    assert len(matched) >= n_treated * 0.9


def test_match_on_propensity_empty_when_no_controls():
    df = _ps_df(ids=[1, 2], treatment=[1, 1], ps=[0.5, 0.5])
    matched = match_on_propensity(df, "is_tricyclic", "propensity_score")
    assert len(matched) == 0


def test_match_on_propensity_is_deterministic_given_a_seed():
    df = _ps_df(
        ids=list(range(20)),
        treatment=[1] * 10 + [0] * 10,
        ps=list(np.linspace(0.3, 0.7, 10)) + list(np.linspace(0.3, 0.7, 10)),
    )
    m1 = match_on_propensity(df, "is_tricyclic", "propensity_score", caliper=1.0, seed=42)
    m2 = match_on_propensity(df, "is_tricyclic", "propensity_score", caliper=1.0, seed=42)
    pd.testing.assert_frame_equal(m1.sort_values("treated_hadm_id").reset_index(drop=True),
                                   m2.sort_values("treated_hadm_id").reset_index(drop=True))
