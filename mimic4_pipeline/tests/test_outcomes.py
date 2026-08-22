"""Integration tests for outcome extraction against the built demo
database. Skipped if the database hasn't been built yet (see
scripts/build_database.sh) rather than failing -- this suite is about
validating extraction logic against real data shapes, not something that
should block a clean checkout with no data.
"""

from __future__ import annotations

import duckdb
import pytest

from mimic4_pipeline.config import DB_PATH
from mimic4_pipeline.outcomes import (
    ANALYTE_EXTREME_DIRECTION,
    BASELINE_WINDOW_HOURS,
    extract_outcomes,
    pivot_baseline,
    pivot_trajectory,
    resolve_lab_itemids,
)

pytestmark = pytest.mark.skipif(not DB_PATH.exists(), reason="demo database not built (run scripts/build_database.sh)")


@pytest.fixture(scope="module")
def con():
    c = duckdb.connect(str(DB_PATH), read_only=True)
    yield c
    c.close()


@pytest.fixture(scope="module")
def resolved(con):
    return resolve_lab_itemids(con, verbose=False)


def test_resolves_all_four_analytes_to_a_nonempty_itemid_list(resolved):
    assert set(resolved) == set(ANALYTE_EXTREME_DIRECTION)
    for analyte, itemids in resolved.items():
        assert itemids, f"{analyte} resolved to no itemids"


def test_sodium_resolves_to_blood_chemistry_not_whole_blood_gas(con, resolved):
    # The demo has two sodium specimen/method variants (see the printed
    # resolution report); blood-chemistry sodium is the one with data and
    # is the clinically standard serum sodium measurement.
    sql = "select label, fluid, category from mimiciv_hosp.d_labitems where itemid in ({})".format(
        ",".join(str(i) for i in resolved["Sodium"])
    )
    rows = con.execute(sql).fetchdf()
    assert (rows["category"] == "Chemistry").all()


def test_extract_outcomes_baseline_window_is_within_24h(con, resolved):
    long_df = extract_outcomes(con, resolved, verbose=False)
    with_baseline = long_df.dropna(subset=["baseline_hours_from_admit"])
    assert (with_baseline["baseline_hours_from_admit"] >= 0).all()
    assert (with_baseline["baseline_hours_from_admit"] <= BASELINE_WINDOW_HOURS).all()


def test_extract_outcomes_delta_only_computed_when_baseline_present(con, resolved):
    long_df = extract_outcomes(con, resolved, verbose=False)
    no_baseline = long_df[long_df["baseline_value"].isna()]
    assert no_baseline["delta_from_baseline"].isna().all()
    has_baseline = long_df[long_df["baseline_value"].notna()]
    assert has_baseline["delta_from_baseline"].notna().all()


def test_extract_outcomes_delta_sign_matches_extreme_minus_baseline(con, resolved):
    long_df = extract_outcomes(con, resolved, verbose=False)
    has_baseline = long_df[long_df["baseline_value"].notna()]
    expected = has_baseline["extreme_value"] - has_baseline["baseline_value"]
    assert (has_baseline["delta_from_baseline"].round(6) == expected.round(6)).all()


def test_sodium_extreme_direction_is_min_others_are_max():
    assert ANALYTE_EXTREME_DIRECTION["Sodium"] == "min"
    for analyte in ["ALT", "AST", "Total Bilirubin"]:
        assert ANALYTE_EXTREME_DIRECTION[analyte] == "max"


def test_extreme_x_uln_matches_extreme_value_over_its_own_ref_range(con, resolved):
    long_df = extract_outcomes(con, resolved, verbose=False)
    hepatic = long_df[long_df["analyte"].isin(["ALT", "AST", "Total Bilirubin"])]
    with_uln = hepatic.dropna(subset=["extreme_x_uln", "extreme_ref_range_upper"])
    assert len(with_uln) > 0
    expected = with_uln["extreme_value"] / with_uln["extreme_ref_range_upper"]
    assert (with_uln["extreme_x_uln"].round(6) == expected.round(6)).all()


def test_pivot_baseline_and_trajectory_are_one_row_per_hadm_id(con, resolved):
    long_df = extract_outcomes(con, resolved, verbose=False)
    baseline_wide = pivot_baseline(long_df)
    trajectory_wide = pivot_trajectory(long_df)
    assert baseline_wide["hadm_id"].is_unique
    assert trajectory_wide["hadm_id"].is_unique
    assert set(baseline_wide["hadm_id"]) == set(long_df["hadm_id"])
