"""Task 9: reviewer-requested sensitivity analyses.

Four checks on the primary Tricyclic-vs-(SSRI+SNRI) comparison, each
re-using the pre-specified thresholds (:data:`CLINICAL_THRESHOLDS`) and the
Charlson-adjusted logistic model (:func:`adjusted_odds_ratio`) unchanged:

  1. **Exposure-before-outcome.** The task-4 trajectory takes the extreme
     value over the *whole* admission, so part of it can precede the first
     antidepressant order. Here the extreme is re-derived from lab draws
     charted at or after the admission's first systemic order of its
     assigned exposure class ("post-exposure"), and again restricted to
     admissions not already past the threshold before that order
     ("incident").
  2. **Low-dose tricyclic exclusion.** Low-dose amitriptyline/doxepin
     (pain/sleep use) stays in the single-class Tricyclic arm under the
     task-3 rule; this drops Tricyclic admissions whose dose stratum is
     ``low_dose``.
  3. **Matched-cohort threshold odds ratios**, from the task-8 matched
     pairs.
  4. **Baseline sodium by group**, to make the sign of the sodium
     delta-from-baseline interpretable alongside threshold crossing.
"""

from __future__ import annotations

import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mimic4_pipeline.analysis import (
    CLINICAL_THRESHOLDS,
    PRIMARY_ANALYTES,
    _crosses_threshold,
    adjusted_odds_ratio,
)
from mimic4_pipeline.config import DB_PATH, OUTPUT_DIR
from mimic4_pipeline.outcomes import ANALYTE_EXTREME_DIRECTION, resolve_lab_itemids

SENSITIVITY_CSV = OUTPUT_DIR / "analysis_reviewer_sensitivity.csv"
TIMING_CSV = OUTPUT_DIR / "analysis_extreme_before_exposure.csv"
BASELINE_CSV = OUTPUT_DIR / "analysis_baseline_sodium.csv"

CLASS_LABEL = {"Tricyclic": "Tricyclic", "SSRI": "SSRI", "SNRI": "SNRI"}


def _base_cohort(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    cohort = con.execute(
        "select c.hadm_id, c.exposure_class, c.dose_stratum, c.admittime, s.charlson_comorbidity_index "
        "from mimiciv_analysis.cohort c left join mimiciv_analysis.severity s using (hadm_id) "
        "where c.exposure_class in ('Tricyclic', 'SSRI', 'SNRI')"
    ).fetchdf()
    cohort["is_tricyclic"] = (cohort["exposure_class"] == "Tricyclic").astype(int)
    return cohort


def first_exposure_times(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """Earliest systemic order of the admission's assigned exposure class."""
    return con.execute(
        """
        select c.hadm_id, min(p.starttime) as first_exposure_time
        from mimiciv_analysis.cohort c
        join mimiciv_analysis.antidepressant_prescriptions p
          on p.hadm_id = c.hadm_id and p.is_systemic and p.drug_class = c.exposure_class
        where c.exposure_class in ('Tricyclic', 'SSRI', 'SNRI') and p.starttime is not null
        group by c.hadm_id
        """
    ).fetchdf()


def post_exposure_outcomes(con: duckdb.DuckDBPyConnection, itemids: dict[str, list[int]]) -> pd.DataFrame:
    """Long format, one row per (hadm_id, analyte): the extreme value (and
    its x ULN) over draws at/after first exposure, plus whether any draw
    *before* first exposure already crossed the threshold."""
    con.register("first_exposure", first_exposure_times(con))
    frames = []
    for analyte in PRIMARY_ANALYTES:
        ids = ",".join(str(i) for i in itemids[analyte])
        df = con.execute(
            f"""
            select le.hadm_id, le.charttime, le.valuenum, le.ref_range_upper,
                   le.charttime >= f.first_exposure_time as post
            from mimiciv_hosp.labevents le
            join first_exposure f using (hadm_id)
            where le.itemid in ({ids}) and le.valuenum is not null
            """
        ).fetchdf()
        spec = CLINICAL_THRESHOLDS[analyte]
        df["extreme_x_uln"] = df["valuenum"] / df["ref_range_upper"].replace(0, np.nan)
        df["extreme_value"] = df["valuenum"]
        df["crosses"] = df.apply(lambda r: _crosses_threshold(r, spec), axis=1)

        pre = df[~df["post"]].groupby("hadm_id")["crosses"].max().rename("crossed_before_exposure")
        post = df[df["post"]].copy()
        if ANALYTE_EXTREME_DIRECTION[analyte] == "max":
            idx = post.groupby("hadm_id")["valuenum"].idxmax()
        else:
            idx = post.groupby("hadm_id")["valuenum"].idxmin()
        ext = post.loc[idx, ["hadm_id", "extreme_value", "extreme_x_uln", "crosses"]]
        ext = ext.merge(pre, on="hadm_id", how="left")
        ext["crossed_before_exposure"] = ext["crossed_before_exposure"].astype("boolean").fillna(False).astype(bool)
        ext["analyte"] = analyte
        frames.append(ext)
    return pd.concat(frames, ignore_index=True)


def extreme_before_exposure(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """Share of whole-admission extremes (task 4) charted before the first
    order of the assigned class -- the size of the timing problem."""
    df = con.execute(
        "select o.hadm_id, o.analyte, o.time_to_extreme_hours, c.exposure_class, c.admittime "
        "from mimiciv_analysis.outcomes_long o join mimiciv_analysis.cohort c using (hadm_id) "
        "where c.exposure_class in ('Tricyclic', 'SSRI', 'SNRI')"
    ).fetchdf().merge(first_exposure_times(con), on="hadm_id")
    exposure_h = (df["first_exposure_time"] - df["admittime"]).dt.total_seconds() / 3600
    df["extreme_before_exposure"] = df["time_to_extreme_hours"] < exposure_h
    df["group"] = np.where(df["exposure_class"] == "Tricyclic", "Tricyclic", "SSRI+SNRI")
    return (
        df.groupby(["analyte", "group"])["extreme_before_exposure"]
        .agg(n="size", proportion="mean")
        .reset_index()
    )


def _fit(sub: pd.DataFrame, analysis: str, analyte: str) -> dict:
    sub = sub.dropna(subset=["crosses"]).copy()
    sub["crosses"] = sub["crosses"].astype(int)
    fit = adjusted_odds_ratio(sub, "crosses", "is_tricyclic", ["charlson_comorbidity_index"])
    fit["analysis"], fit["analyte"] = analysis, analyte
    return fit


def _matched_or(sub: pd.DataFrame, analysis: str, analyte: str) -> dict:
    """Unadjusted logistic OR within the matched cohort (balance is by design)."""
    sub = sub.dropna(subset=["crosses"]).copy()
    sub["crosses"] = sub["crosses"].astype(int)
    model = smf.logit("crosses ~ is_tricyclic", data=sub).fit(disp=0)
    ci = model.conf_int().loc["is_tricyclic"]
    return {
        "analysis": analysis,
        "analyte": analyte,
        "n_treated": int(sub["is_tricyclic"].sum()),
        "n_untreated": int((sub["is_tricyclic"] == 0).sum()),
        "n_events_treated": int(sub.loc[sub["is_tricyclic"] == 1, "crosses"].sum()),
        "n_events_untreated": int(sub.loc[sub["is_tricyclic"] == 0, "crosses"].sum()),
        "adjusted_or": float(np.exp(model.params["is_tricyclic"])),
        "adjusted_or_ci_low": float(np.exp(ci[0])),
        "adjusted_or_ci_high": float(np.exp(ci[1])),
        "adjusted_p_value": float(model.pvalues["is_tricyclic"]),
        "status": "ok",
    }


def run_sensitivity(con: duckdb.DuckDBPyConnection, verbose: bool = True) -> pd.DataFrame:
    cohort = _base_cohort(con)
    whole = con.execute("select * from mimiciv_analysis.outcomes_long").fetchdf().merge(cohort, on="hadm_id")
    post = post_exposure_outcomes(con, resolve_lab_itemids(con, verbose=False)).merge(cohort, on="hadm_id")
    pairs = pd.read_csv(OUTPUT_DIR / "propensity_matched_pairs.csv")
    matched_ids = set(pairs["treated_hadm_id"]) | set(pairs["control_hadm_id"])

    rows = []
    for analyte in PRIMARY_ANALYTES:
        spec = CLINICAL_THRESHOLDS[analyte]
        w = whole[whole["analyte"] == analyte].copy()
        w["crosses"] = w.apply(lambda r: _crosses_threshold(r, spec), axis=1)
        p = post[post["analyte"] == analyte]

        rows.append(_fit(w, "whole_admission (original)", analyte))
        rows.append(_fit(p, "post_exposure", analyte))
        rows.append(_fit(p[~p["crossed_before_exposure"]], "post_exposure_incident", analyte))
        no_low = ~((w["exposure_class"] == "Tricyclic") & (w["dose_stratum"] == "low_dose"))
        rows.append(_fit(w[no_low], "excl_low_dose_tricyclic", analyte))
        p_no_low = ~((p["exposure_class"] == "Tricyclic") & (p["dose_stratum"] == "low_dose"))
        rows.append(_fit(p[p_no_low], "post_exposure_excl_low_dose_tricyclic", analyte))
        rows.append(_matched_or(w[w["hadm_id"].isin(matched_ids)], "propensity_matched", analyte))

    cols = ["analysis", "analyte", "n_treated", "n_untreated", "n_events_treated", "n_events_untreated",
            "adjusted_or", "adjusted_or_ci_low", "adjusted_or_ci_high", "adjusted_p_value", "status"]
    df = pd.DataFrame(rows)[cols]
    if verbose:
        print("=== Reviewer sensitivity analyses: OR(Tricyclic vs SSRI+SNRI) for threshold crossing ===")
        print(df.to_string(index=False))
    return df


def baseline_sodium(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    return con.execute(
        """
        select case when c.exposure_class = 'Tricyclic' then 'Tricyclic' else 'SSRI+SNRI' end as grp,
               count(o.baseline_value) as n, avg(o.baseline_value) as mean_baseline,
               avg(o.extreme_value) as mean_extreme, avg(o.delta_from_baseline) as mean_delta
        from mimiciv_analysis.outcomes_long o join mimiciv_analysis.cohort c using (hadm_id)
        where o.analyte = 'Sodium' and c.exposure_class in ('Tricyclic', 'SSRI', 'SNRI')
          and o.baseline_value is not null
        group by 1
        """
    ).fetchdf()


def main() -> int:
    con = duckdb.connect(str(DB_PATH), read_only=True)
    try:
        timing = extreme_before_exposure(con)
        print("=== Whole-admission extreme charted before first exposure ===")
        print(timing.to_string(index=False))
        timing.to_csv(TIMING_CSV, index=False)

        base = baseline_sodium(con)
        print("\n=== Sodium baseline / extreme / delta by group ===")
        print(base.to_string(index=False))
        base.to_csv(BASELINE_CSV, index=False)

        df = run_sensitivity(con)
        df.to_csv(SENSITIVITY_CSV, index=False)
        print(f"\nWrote {len(df)} rows -> {SENSITIVITY_CSV}")
    finally:
        con.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
