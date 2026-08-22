from __future__ import annotations

import math

import pandas as pd
import pytest

from mimic4_pipeline.dose import (
    DOSE_THRESHOLDS_MG_PER_DAY,
    class_is_dose_ambiguous_low,
    classify_dose_stratum,
    estimate_daily_dose_mg,
    parse_dose_mg,
)


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("50", 50.0),
        ("7.5", 7.5),
        ("25-50", 50.0),  # range -> upper bound (conservative)
        ("300", 300.0),
        (None, None),
        ("", None),
        ("   ", None),
        ("QID", None),
        (float("nan"), None),
    ],
)
def test_parse_dose_mg(raw, expected):
    result = parse_dose_mg(raw)
    if expected is None:
        assert result is None
    else:
        assert result == expected


def test_estimate_daily_dose_mg_multiplies_by_frequency():
    assert estimate_daily_dose_mg("50", 2.0) == 100.0


@pytest.mark.parametrize("doses_per_24_hrs", [None, 0, 0.0, float("nan")])
def test_estimate_daily_dose_mg_treats_missing_frequency_as_once_daily(doses_per_24_hrs):
    assert estimate_daily_dose_mg("50", doses_per_24_hrs) == 50.0


def test_estimate_daily_dose_mg_unparseable_dose_is_none():
    assert estimate_daily_dose_mg("QID", 1.0) is None


def test_dose_thresholds_cover_the_three_ambiguous_agents():
    assert set(DOSE_THRESHOLDS_MG_PER_DAY) == {"amitriptyline", "doxepin", "trazodone"}


def test_classify_dose_stratum_not_applicable_for_unlisted_agent():
    assert classify_dose_stratum("sertraline", [200.0]) == "not_applicable"


def test_classify_dose_stratum_low_dose():
    assert classify_dose_stratum("trazodone", [25.0, 50.0]) == "low_dose"


def test_classify_dose_stratum_standard_dose_uses_max_over_admission():
    # Dose escalated from 25 (low) to 300 (standard) during the stay --
    # the admission is classified by the maximum, i.e. the strongest
    # evidence it reached antidepressant-range dosing.
    assert classify_dose_stratum("trazodone", [25.0, 300.0]) == "standard_dose"


def test_classify_dose_stratum_unknown_when_nothing_parses():
    assert classify_dose_stratum("amitriptyline", [None, None]) == "unknown_dose"


def test_classify_dose_stratum_unknown_when_values_are_nan_not_none():
    # Regression: pandas silently turns None into float NaN inside a float
    # column (e.g. after DataFrame.groupby(...)[...].tolist()) -- a naive
    # `is not None` filter lets NaN through, and max([nan, ...]) is
    # unreliable, both of which previously let an unparseable dose slip
    # through as a real (and always-False-threshold-comparison, hence
    # spuriously "low_dose") value.
    assert classify_dose_stratum("amitriptyline", [float("nan"), float("nan")]) == "unknown_dose"
    assert classify_dose_stratum("amitriptyline", pd.Series([float("nan")]).tolist()) == "unknown_dose"


def test_classify_dose_stratum_boundary_is_inclusive_of_standard():
    assert classify_dose_stratum("amitriptyline", [75.0]) == "standard_dose"
    assert classify_dose_stratum("amitriptyline", [74.9]) == "low_dose"


# --- class_is_dose_ambiguous_low ---------------------------------------


def test_class_dose_ambiguous_low_true_when_sole_agent_is_low_dose():
    assert class_is_dose_ambiguous_low(["trazodone"], {"trazodone": "low_dose"}) is True


def test_class_dose_ambiguous_low_false_when_agent_is_standard_dose():
    assert class_is_dose_ambiguous_low(["trazodone"], {"trazodone": "standard_dose"}) is False


def test_class_dose_ambiguous_low_false_when_dose_unknown():
    # Unknown dose is not evidence of low dose -- conservative default.
    assert class_is_dose_ambiguous_low(["amitriptyline"], {"amitriptyline": "unknown_dose"}) is False


def test_class_dose_ambiguous_low_false_when_a_non_ambiguous_agent_present():
    # bupropion has no dose ambiguity -- its presence in a class is always
    # real antidepressant-dose evidence, even alongside low-dose trazodone.
    assert (
        class_is_dose_ambiguous_low(
            ["trazodone", "bupropion"], {"trazodone": "low_dose", "bupropion": "not_applicable"}
        )
        is False
    )


def test_class_dose_ambiguous_low_false_when_empty():
    assert class_is_dose_ambiguous_low([], {}) is False


def test_class_dose_ambiguous_low_true_when_multiple_ambiguous_agents_all_low():
    assert (
        class_is_dose_ambiguous_low(
            ["amitriptyline", "doxepin"], {"amitriptyline": "low_dose", "doxepin": "low_dose"}
        )
        is True
    )
