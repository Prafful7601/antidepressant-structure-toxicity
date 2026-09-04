"""Propensity-score-matched sensitivity analysis for the primary Tricyclic
vs (SSRI+SNRI) comparison.

Why this exists: the severity-adjusted comparison (task 7,
:mod:`mimic4_pipeline.analysis`) adjusts for one covariate (Charlson index).
Confounding by indication -- why a clinician chose a tricyclic over an
SSRI/SNRI in the first place -- is not resolved by that alone, and the
analysis module's own docstring says so. This module matches
Tricyclic-exposed admissions to SSRI/SNRI-exposed admissions on a full
propensity model (age, sex, emergency-admission status, and every
individual Charlson comorbidity component -- not just the summary index),
reports the standard balance diagnostic (standardized mean difference,
before and after matching), and re-runs the primary comparison on the
matched cohort so the unadjusted, Charlson-adjusted, and matched results
can all be read side by side.

Run standalone:

    python3 -m mimic4_pipeline.propensity
"""

from __future__ import annotations

import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mimic4_pipeline.analysis import CLINICAL_THRESHOLDS, PRIMARY_ANALYTES, _crosses_threshold, safe_two_group_comparison
from mimic4_pipeline.config import DB_PATH, OUTPUT_DIR, SQL_DIR

COVARIATE_COLS: list[str] = [
    "age",
    "is_male",
    "is_emergency_admission",
    "myocardial_infarct",
    "congestive_heart_failure",
    "peripheral_vascular_disease",
    "cerebrovascular_disease",
    "dementia",
    "chronic_pulmonary_disease",
    "rheumatic_disease",
    "peptic_ulcer_disease",
    "mild_liver_disease",
    "diabetes_without_cc",
    "diabetes_with_cc",
    "paraplegia",
    "renal_disease",
    "malignant_cancer",
    "severe_liver_disease",
    "metastatic_solid_tumor",
    "aids",
]

MATCHED_CSV = OUTPUT_DIR / "propensity_matched_pairs.csv"
BALANCE_CSV = OUTPUT_DIR / "propensity_balance.csv"
MATCHED_COMPARISON_CSV = OUTPUT_DIR / "analysis_matched_comparison.csv"
MATCHED_THRESHOLD_CSV = OUTPUT_DIR / "analysis_matched_threshold_crossing.csv"


def _read_sql(name: str) -> str:
    return (SQL_DIR / "propensity" / name).read_text()


def build_covariates(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """One row per admission: propensity covariates joined to the primary
    exposure_class (restricted to Tricyclic/SSRI/SNRI, mirroring
    :func:`mimic4_pipeline.analysis.run_primary_comparison`)."""
    covariates = con.execute(_read_sql("covariates.sql")).fetchdf()
    cohort = con.execute(
        "select hadm_id, exposure_class from mimiciv_analysis.cohort "
        "where exposure_class in ('Tricyclic', 'SSRI', 'SNRI')"
    ).fetchdf()
    df = cohort.merge(covariates, on="hadm_id", how="left")
    df["is_tricyclic"] = (df["exposure_class"] == "Tricyclic").astype(int)
    return df


def fit_propensity_scores(df: pd.DataFrame, treatment_col: str, covariate_cols: list[str]) -> pd.DataFrame:
    """Fit ``treatment_col ~ covariate_cols`` via logistic regression and
    return ``df`` with a ``propensity_score`` column added. Rows with any
    missing covariate are dropped (reported by the caller via row counts,
    not silently)."""
    d = df.dropna(subset=[treatment_col, *covariate_cols]).copy()
    formula = f"{treatment_col} ~ " + " + ".join(covariate_cols)
    model = smf.logit(formula, data=d).fit(disp=0, maxiter=200)
    d["propensity_score"] = model.predict(d)
    return d


def standardized_mean_diff(df: pd.DataFrame, treatment_col: str, covariate_cols: list[str]) -> pd.DataFrame:
    """Standard PSM balance diagnostic: for each covariate, the difference
    in means between treated and control, scaled by the pooled standard
    deviation. |SMD| < 0.1 is the conventional "well balanced" threshold;
    this is exactly what should shrink after matching if the match worked.
    """
    rows = []
    treated = df[df[treatment_col] == 1]
    control = df[df[treatment_col] == 0]
    for col in covariate_cols:
        t, c = treated[col].astype(float), control[col].astype(float)
        pooled_sd = np.sqrt((t.var(ddof=1) + c.var(ddof=1)) / 2)
        smd = (t.mean() - c.mean()) / pooled_sd if pooled_sd > 0 else 0.0
        rows.append({"covariate": col, "mean_treated": t.mean(), "mean_control": c.mean(), "smd": smd})
    return pd.DataFrame(rows)


def match_on_propensity(
    df: pd.DataFrame,
    treatment_col: str,
    ps_col: str,
    id_col: str = "hadm_id",
    caliper: float | None = None,
    batch_size: int = 25,
    seed: int = 0,
) -> pd.DataFrame:
    """1:1 nearest-neighbor matching without replacement, on the logit of
    the propensity score (the standard scale to match on -- raw
    probabilities compress extreme propensities together).

    Processed in random-order batches, rebuilding a KD-tree over the
    still-available controls after each batch: an efficient approximation
    to exact sequential greedy matching (rebuilding after every single
    match is also correct, just much slower at this scale -- ~7k treated
    vs ~80k controls) with a small enough batch size that within-batch
    collisions (two treated units both nearest to the same control) are
    rare, not eliminated by construction.

    ``caliper`` defaults to 0.2 x the SD of the logit propensity score
    (Austin 2011's commonly cited convention); a treated unit with no
    control within the caliper goes unmatched, not force-matched to a
    poor control.
    """
    eps = 1e-6
    d = df[[id_col, treatment_col, ps_col]].dropna().copy()
    p = d[ps_col].clip(eps, 1 - eps)
    d["logit_ps"] = np.log(p / (1 - p))

    if caliper is None:
        caliper = 0.2 * d["logit_ps"].std()

    treated = d[d[treatment_col] == 1].sample(frac=1, random_state=seed)
    control = d[d[treatment_col] == 0]

    treated_ids = treated[id_col].to_numpy()
    treated_lps = treated["logit_ps"].to_numpy()
    control_ids = control[id_col].to_numpy()
    control_lps = control["logit_ps"].to_numpy()
    used = np.zeros(len(control_ids), dtype=bool)

    pairs = []
    for start in range(0, len(treated_ids), batch_size):
        avail_idx = np.flatnonzero(~used)
        if len(avail_idx) == 0:
            break
        tree = cKDTree(control_lps[avail_idx].reshape(-1, 1))
        batch_ids = treated_ids[start : start + batch_size]
        batch_lps = treated_lps[start : start + batch_size]
        dist, idx = tree.query(batch_lps.reshape(-1, 1), k=1)
        for t_id, dd, ii in zip(batch_ids, np.atleast_1d(dist), np.atleast_1d(idx)):
            if dd > caliper:
                continue
            c_pos = avail_idx[ii]
            if used[c_pos]:
                continue  # claimed by another unit earlier in this same batch
            used[c_pos] = True
            pairs.append({"treated_hadm_id": t_id, "control_hadm_id": control_ids[c_pos], "logit_ps_distance": dd})

    return pd.DataFrame(pairs, columns=["treated_hadm_id", "control_hadm_id", "logit_ps_distance"])


def _matched_hadm_ids(matched_pairs: pd.DataFrame) -> tuple[set, set]:
    return set(matched_pairs["treated_hadm_id"]), set(matched_pairs["control_hadm_id"])


def matched_primary_comparison(
    con: duckdb.DuckDBPyConnection, matched_pairs: pd.DataFrame, verbose: bool = True
) -> pd.DataFrame:
    """Re-runs the primary delta-from-baseline comparison
    (:func:`mimic4_pipeline.analysis.safe_two_group_comparison`) restricted
    to the matched cohort."""
    treated_ids, control_ids = _matched_hadm_ids(matched_pairs)
    outcomes = con.execute("select * from mimiciv_analysis.outcomes_long").fetchdf()

    rows = []
    if verbose:
        print(f"=== Matched comparison: Tricyclic vs (SSRI+SNRI), {len(treated_ids)} matched pairs ===")
    for analyte in PRIMARY_ANALYTES:
        sub = outcomes[outcomes["analyte"] == analyte]
        a_vals = sub.loc[sub["hadm_id"].isin(treated_ids), "delta_from_baseline"].tolist()
        b_vals = sub.loc[sub["hadm_id"].isin(control_ids), "delta_from_baseline"].tolist()
        result = safe_two_group_comparison(a_vals, b_vals, analyte)
        rows.append(result)
        if verbose:
            if result["status"] == "ok":
                print(f"  {analyte:16s} n_a={result['n_group_a']:4d} n_b={result['n_group_b']:4d} "
                      f"mean_diff={result['mean_diff_a_minus_b']:+.2f} "
                      f"[{result['mean_diff_ci_low']:+.2f}, {result['mean_diff_ci_high']:+.2f}] "
                      f"p(secondary)={result['p_value_mannwhitney_secondary']:.4f}")
            else:
                print(f"  {analyte:16s} SKIPPED -- {result['status']}")
    return pd.DataFrame(rows)


def matched_threshold_crossing(
    con: duckdb.DuckDBPyConnection, matched_pairs: pd.DataFrame, verbose: bool = True
) -> pd.DataFrame:
    """Re-runs the clinical-threshold-crossing proportions on the matched
    cohort (same pre-specified thresholds as task 7, not re-decided here)."""
    from mimic4_pipeline.analysis import _wilson_ci

    treated_ids, control_ids = _matched_hadm_ids(matched_pairs)
    outcomes = con.execute("select * from mimiciv_analysis.outcomes_long").fetchdf()

    rows = []
    for analyte, spec in CLINICAL_THRESHOLDS.items():
        sub = outcomes[outcomes["analyte"] == analyte].copy()
        sub["crosses"] = sub.apply(lambda r: _crosses_threshold(r, spec), axis=1)
        for label, ids in [("Tricyclic", treated_ids), ("SSRI+SNRI", control_ids)]:
            g = sub[sub["hadm_id"].isin(ids)].dropna(subset=["crosses"])
            n, count = len(g), int(g["crosses"].sum())
            lo, hi = _wilson_ci(count, n)
            rows.append(
                {
                    "analyte": analyte,
                    "threshold": spec["description"],
                    "group": label,
                    "n_evaluable": n,
                    "n_crossing": count,
                    "proportion_crossing": (count / n) if n else None,
                    "proportion_ci_low": lo if n else None,
                    "proportion_ci_high": hi if n else None,
                }
            )
    df = pd.DataFrame(rows)
    if verbose:
        print("\n=== Matched clinical threshold crossing ===")
        print(df.to_string(index=False))
    return df


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(DB_PATH), read_only=True)
    try:
        covariates = build_covariates(con)
        n_before = len(covariates)
        fitted = fit_propensity_scores(covariates, "is_tricyclic", COVARIATE_COLS)
        n_dropped = n_before - len(fitted)
        print(f"Propensity model: {len(fitted)} admissions ({n_dropped} dropped for missing covariates), "
              f"{int(fitted['is_tricyclic'].sum())} Tricyclic / {int((fitted['is_tricyclic']==0).sum())} SSRI+SNRI")

        balance_before = standardized_mean_diff(fitted, "is_tricyclic", COVARIATE_COLS)

        matched_pairs = match_on_propensity(fitted, "is_tricyclic", "propensity_score")
        matched_pairs.to_csv(MATCHED_CSV, index=False)
        treated_ids, control_ids = _matched_hadm_ids(matched_pairs)
        n_treated_total = int(fitted["is_tricyclic"].sum())
        print(f"\nMatched {len(matched_pairs)} of {n_treated_total} Tricyclic admissions "
              f"({len(matched_pairs) / n_treated_total:.1%}) within caliper.")

        matched_fitted = fitted[fitted["hadm_id"].isin(treated_ids | control_ids)]
        balance_after = standardized_mean_diff(matched_fitted, "is_tricyclic", COVARIATE_COLS)

        balance = balance_before.merge(balance_after, on="covariate", suffixes=("_before", "_after"))
        balance.to_csv(BALANCE_CSV, index=False)
        print("\n=== Covariate balance (standardized mean difference; |SMD|<0.1 is conventionally well-balanced) ===")
        print(balance[["covariate", "smd_before", "smd_after"]].round(3).to_string(index=False))
        n_imbalanced_before = int((balance["smd_before"].abs() >= 0.1).sum())
        n_imbalanced_after = int((balance["smd_after"].abs() >= 0.1).sum())
        print(f"\nCovariates with |SMD|>=0.1: {n_imbalanced_before} before matching -> {n_imbalanced_after} after matching")

        comparison_df = matched_primary_comparison(con, matched_pairs)
        comparison_df.to_csv(MATCHED_COMPARISON_CSV, index=False)
        print(f"\nWrote {len(comparison_df)} rows -> {MATCHED_COMPARISON_CSV}")

        threshold_df = matched_threshold_crossing(con, matched_pairs)
        threshold_df.to_csv(MATCHED_THRESHOLD_CSV, index=False)
        print(f"Wrote {len(threshold_df)} rows -> {MATCHED_THRESHOLD_CSV}")
    finally:
        con.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
