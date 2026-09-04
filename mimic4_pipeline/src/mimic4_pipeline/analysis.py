"""Task 7: analysis skeleton.

Primary comparison: tricyclic-exposed vs (SSRI + SNRI)-exposed admissions,
adjusted for admission severity, on each of the 4 outcome analytes'
delta-from-baseline (task 4) -- absolute values are dominated by the
admission's primary reason, not exposure, per the brief. Uses
``mimiciv_analysis.cohort.exposure_class``, the dose-reclassified primary
assignment (task 3), not the raw multi-class-exclusion cohort.

This is written to run end-to-end on the demo **even though one comparison
group is empty there** (task 3's Tricyclic arm stays at zero even after
dose-based reclassification -- the demo's 100 patients contain no
standard-dose amitriptyline order at all, single- or multi-class; see
``mimic4_pipeline.cohort``'s module docstring) -- every statistical step is
guarded and reports *why* it did not run rather than crashing or emitting a
spurious p-value. See :func:`safe_two_group_comparison`.

Effect sizes with confidence intervals are the PRIMARY result; p-values are
secondary (2026-08-23 revision). The feasibility calculation
(:func:`feasibility_projection`) found a minimum detectable Cohen's d of
about 0.034 at full-database scale -- small enough that statistical
significance is essentially guaranteed there, which means **p-values alone
would not indicate clinical relevance at that scale, and the binding
constraint on this study is confounding control, not statistical power**.
Two things follow, both implemented here:
  1. Report the effect size itself (mean difference in delta-from-baseline,
     with a 95% CI) as primary, not just whether p<0.05.
  2. Pre-specify a clinically meaningful effect threshold per analyte
     BEFORE looking at any comparison result (:data:`CLINICAL_THRESHOLDS`),
     so "clinically meaningful" cannot be dredged post hoc from whatever
     the data happen to show. ALT/AST/bilirubin use the standard
     Hy's-Law-style multiple-of-ULN hepatic injury signal (task 4 already
     carries each measurement's own reference range); sodium uses the
     standard absolute threshold for moderate hyponatremia.
"""

from __future__ import annotations

import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from scipy import stats as sstats
import statsmodels.formula.api as smf
from statsmodels.stats.power import TTestIndPower
from statsmodels.stats.proportion import proportion_confint

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mimic4_pipeline.config import DB_PATH, OUTPUT_DIR

# PhysioNet, physionet.org/content/mimiciv/ (MIMIC-IV v3.1, released
# 2024-10-11), accessed 2026-08-23: "546,028 hospitalizations". This is an
# external, cited assumption for the feasibility projection below, not a
# number derived from the demo -- re-run this projection against the real
# admission count once full access is granted and the actual database is
# built (config.py's DB_PATH is the only thing that needs to change to do
# that; this constant should then be replaced with the true count, e.g.
# `select count(*) from mimiciv_hosp.admissions`).
FULL_MIMIC_IV_ADMISSIONS = 546_028

PRIMARY_ANALYTES = ["ALT", "AST", "Total Bilirubin", "Sodium"]
MIN_GROUP_SIZE_FOR_TEST = 2  # below this, neither an effect-size CI nor a p-value is meaningful

# Pre-specified clinically meaningful effect thresholds, decided before
# looking at any comparison group's outcomes (see module docstring).
# ALT/AST/bilirubin: standard Hy's-Law-style multiple-of-ULN hepatic injury
# criteria (Hy's Law itself pairs >=3x ULN transaminases with >=2x ULN
# bilirubin as the drug-induced-liver-injury signal threshold used in FDA
# guidance). Sodium: standard absolute threshold for moderate,
# clinically-actionable hyponatremia -- there is no "multiple of ULN"
# convention for sodium.
CLINICAL_THRESHOLDS: dict[str, dict] = {
    "ALT": {"metric": "extreme_x_uln", "op": ">=", "value": 3.0, "description": ">=3x ULN (Hy's Law hepatocellular injury threshold)"},
    "AST": {"metric": "extreme_x_uln", "op": ">=", "value": 3.0, "description": ">=3x ULN (Hy's Law hepatocellular injury threshold)"},
    "Total Bilirubin": {"metric": "extreme_x_uln", "op": ">=", "value": 2.0, "description": ">=2x ULN (Hy's Law hyperbilirubinemia threshold)"},
    "Sodium": {"metric": "extreme_value", "op": "<=", "value": 130.0, "description": "<=130 mEq/L (moderate hyponatremia)"},
}

ANALYSIS_CSV = OUTPUT_DIR / "analysis_primary_comparison.csv"
FEASIBILITY_CSV = OUTPUT_DIR / "feasibility_projection.csv"


def _clean(values) -> np.ndarray:
    return np.asarray(
        [x for x in values if x is not None and not (isinstance(x, float) and x != x)], dtype=float
    )


def safe_two_group_comparison(a, b, label: str) -> dict:
    """Effect size (mean difference, 95% CI) as the PRIMARY result; a
    Mann-Whitney U p-value as a secondary, distribution-robust check.
    Guarded against every degenerate case instead of crashing or emitting a
    spurious result:
      - either group empty
      - either group below MIN_GROUP_SIZE_FOR_TEST
      - all values identical in both groups (zero variance)
    """
    a = _clean(a)
    b = _clean(b)

    result = {
        "outcome": label,
        "n_group_a": len(a),
        "n_group_b": len(b),
        "mean_a": float(np.mean(a)) if len(a) else None,
        "mean_b": float(np.mean(b)) if len(b) else None,
        # Reported alongside the mean every time, not just when something
        # looks off: a mean-difference CI can be dominated by a handful of
        # extreme-outlier admissions in distributions this skewed (delta
        # values in the thousands of U/L do occur), while the median (and
        # the secondary rank-based test) are robust to exactly that. A
        # mean/median divergence is itself diagnostic, not a footnote.
        "median_a": float(np.median(a)) if len(a) else None,
        "median_b": float(np.median(b)) if len(b) else None,
        "mean_diff_a_minus_b": None,
        "mean_diff_ci_low": None,
        "mean_diff_ci_high": None,
        "p_value_mannwhitney_secondary": None,
        "status": None,
    }

    if len(a) == 0 or len(b) == 0:
        result["status"] = f"insufficient_data: group {'a' if len(a) == 0 else 'b'} is empty (n=0)"
        return result
    if len(a) < MIN_GROUP_SIZE_FOR_TEST or len(b) < MIN_GROUP_SIZE_FOR_TEST:
        result["status"] = (
            f"insufficient_data: need >={MIN_GROUP_SIZE_FOR_TEST} per group for a meaningful "
            f"effect-size CI or test (have {len(a)} and {len(b)})"
        )
        return result
    if np.all(a == a[0]) and np.all(b == b[0]) and a[0] == b[0]:
        result["status"] = "degenerate: both groups constant and equal -- no variation to test"
        return result

    # Primary: Welch's t effect size (mean difference) with its 95% CI --
    # does not assume equal variances, appropriate since group sizes and
    # variances differ substantially between exposure classes here.
    t_result = sstats.ttest_ind(a, b, equal_var=False)
    ci = t_result.confidence_interval(confidence_level=0.95)
    result["mean_diff_a_minus_b"] = float(np.mean(a) - np.mean(b))
    result["mean_diff_ci_low"] = float(ci.low)
    result["mean_diff_ci_high"] = float(ci.high)

    # Secondary: Mann-Whitney U p-value (rank-based, robust to the skewed
    # distributions typical of liver enzymes) -- reported as a robustness
    # check on the effect size above, not as the primary result.
    try:
        _, p = sstats.mannwhitneyu(a, b, alternative="two-sided")
        result["p_value_mannwhitney_secondary"] = float(p)
    except ValueError as exc:
        result["status"] = f"effect_size_ok_but_secondary_test_failed: {exc}"
        return result

    result["status"] = "ok"
    return result


def _wilson_ci(count: int, n: int) -> tuple[float, float]:
    if n == 0:
        return (0.0, 0.0)
    lo, hi = proportion_confint(count, n, alpha=0.05, method="wilson")
    return float(lo), float(hi)


def _crosses_threshold(row: pd.Series, spec: dict) -> bool | None:
    value = row.get(spec["metric"])
    if value is None or (isinstance(value, float) and value != value):
        return None
    return value >= spec["value"] if spec["op"] == ">=" else value <= spec["value"]


def clinical_threshold_crossing(merged: pd.DataFrame, verbose: bool = True) -> pd.DataFrame:
    """Proportion of each group crossing its analyte's PRE-SPECIFIED
    clinically meaningful threshold (:data:`CLINICAL_THRESHOLDS`), with a
    Wilson CI per group -- the clinical-relevance counterpart to the raw
    mean-difference effect size in :func:`run_primary_comparison`.
    """
    rows = []
    for analyte, spec in CLINICAL_THRESHOLDS.items():
        sub = merged[merged["analyte"] == analyte].copy()
        sub["crosses"] = sub.apply(lambda r: _crosses_threshold(r, spec), axis=1)
        for group_label, cls_filter in [("Tricyclic", ["Tricyclic"]), ("SSRI+SNRI", ["SSRI", "SNRI"])]:
            g = sub[sub["exposure_class"].isin(cls_filter)].dropna(subset=["crosses"])
            n = len(g)
            count = int(g["crosses"].sum())
            lo, hi = _wilson_ci(count, n)
            rows.append(
                {
                    "analyte": analyte,
                    "threshold": spec["description"],
                    "group": group_label,
                    "n_evaluable": n,
                    "n_crossing": count,
                    "proportion_crossing": (count / n) if n else None,
                    "proportion_ci_low": lo if n else None,
                    "proportion_ci_high": hi if n else None,
                }
            )
    df = pd.DataFrame(rows)
    if verbose:
        print("\n=== Clinical threshold crossing (pre-specified, not data-dredged) ===")
        for analyte, spec in CLINICAL_THRESHOLDS.items():
            print(f"  {analyte}: {spec['description']}")
        print(df.to_string(index=False))
    return df


def run_primary_comparison(
    con: duckdb.DuckDBPyConnection, exclude_overdose: bool = False, verbose: bool = True
) -> pd.DataFrame:
    """Tricyclic vs (SSRI+SNRI), effect size (mean difference + 95% CI,
    primary) and Mann-Whitney p-value (secondary), per analyte, on
    delta_from_baseline.

    ``exclude_overdose`` runs the brief's task-3 "both with and without"
    overdose cohorts (`is_overdose_admission`, task 3) as this comparison,
    not just as an available column nobody reads: TCA/SSRI overdose
    admissions are exactly the kind of case that can produce an extreme,
    admission-reason-driven (not exposure-driven) enzyme spike -- see the
    mean/median divergence this surfaced on the first full-database run,
    documented in the README.

    (Severity-adjusted regression is the natural next step once a
    comparison group has enough rows to fit one -- it does now.)
    """
    cohort = con.execute(
        "select hadm_id, exposure_class, is_overdose_admission from mimiciv_analysis.cohort "
        "where exposure_class is not null"
    ).fetchdf()
    if exclude_overdose:
        cohort = cohort[~cohort["is_overdose_admission"]]
    outcomes = con.execute("select * from mimiciv_analysis.outcomes_long").fetchdf()

    merged = outcomes.merge(cohort, on="hadm_id", how="inner")
    group_a = merged[merged["exposure_class"] == "Tricyclic"]
    group_b = merged[merged["exposure_class"].isin(["SSRI", "SNRI"])]

    if verbose:
        label = "excluding overdose admissions" if exclude_overdose else "including overdose admissions"
        print(f"=== Primary comparison ({label}): Tricyclic vs (SSRI+SNRI), delta from baseline ===")
        print("  Effect size (mean difference, 95% CI) is PRIMARY; the Mann-Whitney p-value is secondary.")
        print(f"  Group A (Tricyclic) admissions available: {group_a['hadm_id'].nunique()}")
        print(f"  Group B (SSRI+SNRI) admissions available: {group_b['hadm_id'].nunique()}")

    rows = []
    for analyte in PRIMARY_ANALYTES:
        a_vals = group_a.loc[group_a["analyte"] == analyte, "delta_from_baseline"].tolist()
        b_vals = group_b.loc[group_b["analyte"] == analyte, "delta_from_baseline"].tolist()
        result = safe_two_group_comparison(a_vals, b_vals, analyte)
        rows.append(result)
        if verbose:
            if result["status"] == "ok":
                print(f"  {analyte:16s} n_a={result['n_group_a']:3d} n_b={result['n_group_b']:3d} "
                      f"mean_diff={result['mean_diff_a_minus_b']:+.2f} "
                      f"[{result['mean_diff_ci_low']:+.2f}, {result['mean_diff_ci_high']:+.2f}] "
                      f"median_a={result['median_a']:+.2f} median_b={result['median_b']:+.2f} "
                      f"p(secondary)={result['p_value_mannwhitney_secondary']:.4f}")
            else:
                print(f"  {analyte:16s} SKIPPED -- {result['status']}")

    threshold_df = clinical_threshold_crossing(merged, verbose=verbose)
    suffix = "_excl_overdose" if exclude_overdose else ""
    threshold_path = OUTPUT_DIR / f"analysis_clinical_threshold_crossing{suffix}.csv"
    threshold_df.to_csv(threshold_path, index=False)
    if verbose:
        print(f"\nWrote {len(threshold_df)} rows -> {threshold_path}")

    return pd.DataFrame(rows)


def adjusted_odds_ratio(
    df: pd.DataFrame, outcome_col: str, treatment_col: str, covariate_cols: list[str]
) -> dict:
    """Core, DB-independent logistic-regression fit: ``outcome_col`` (must
    be 0/1) on ``treatment_col`` (must be 0/1) plus ``covariate_cols``.
    Returns the treatment's adjusted odds ratio with a 95% CI, guarded
    against the same class of degenerate inputs as
    :func:`safe_two_group_comparison` -- separated by design so this can be
    unit-tested without a database (see ``tests/test_analysis.py``).
    """
    result = {
        "n_treated": int(df[treatment_col].sum()) if len(df) else 0,
        "n_untreated": int((df[treatment_col] == 0).sum()) if len(df) else 0,
        "n_events_treated": None,
        "n_events_untreated": None,
        "adjusted_or": None,
        "adjusted_or_ci_low": None,
        "adjusted_or_ci_high": None,
        "adjusted_p_value": None,
        "status": None,
    }
    df = df.dropna(subset=[outcome_col, treatment_col, *covariate_cols])
    n_treated = int(df[treatment_col].sum())
    n_untreated = int((df[treatment_col] == 0).sum())
    result["n_treated"] = n_treated
    result["n_untreated"] = n_untreated

    if n_treated < MIN_GROUP_SIZE_FOR_TEST or n_untreated < MIN_GROUP_SIZE_FOR_TEST:
        result["status"] = f"insufficient_data: n_treated={n_treated}, n_untreated={n_untreated}"
        return result

    result["n_events_treated"] = int(df.loc[df[treatment_col] == 1, outcome_col].sum())
    result["n_events_untreated"] = int(df.loc[df[treatment_col] == 0, outcome_col].sum())

    if df[outcome_col].nunique() < 2:
        result["status"] = "degenerate: no variation in the outcome (all events or none)"
        return result
    # Perfect or quasi-perfect separation (e.g. every event falls in one
    # group) makes the MLE diverge -- check for it explicitly rather than
    # let statsmodels emit an unreliable coefficient with a huge SE.
    for grp in (0, 1):
        sub = df[df[treatment_col] == grp][outcome_col]
        if len(sub) and sub.nunique() < 2:
            result["status"] = f"degenerate: no outcome variation within {treatment_col}={grp}"
            return result

    formula = f"{outcome_col} ~ {treatment_col} + " + " + ".join(covariate_cols)
    try:
        model = smf.logit(formula, data=df).fit(disp=0, maxiter=100)
    except Exception as exc:  # noqa: BLE001 -- report, don't crash on a convergence/separation failure
        result["status"] = f"model_failed: {exc}"
        return result

    if not model.mle_retvals.get("converged", True):
        result["status"] = "did_not_converge"
        return result

    coef = model.params.get(treatment_col)
    if coef is None or not np.isfinite(coef) or abs(coef) > 20:
        # An extreme (near-infinite-odds-ratio) coefficient is the signature
        # of separation statsmodels' convergence check alone doesn't always
        # catch -- treat it as a failure rather than report exp(20)-ish as
        # if it were a real effect size.
        result["status"] = "unstable_fit: coefficient diverged (likely separation)"
        return result

    ci = model.conf_int().loc[treatment_col]
    result["adjusted_or"] = float(np.exp(coef))
    result["adjusted_or_ci_low"] = float(np.exp(ci[0]))
    result["adjusted_or_ci_high"] = float(np.exp(ci[1]))
    result["adjusted_p_value"] = float(model.pvalues[treatment_col])
    result["status"] = "ok"
    return result


def severity_adjusted_comparison(
    con: duckdb.DuckDBPyConnection, exclude_overdose: bool = False, verbose: bool = True
) -> pd.DataFrame:
    """Severity-adjusted Tricyclic-vs-(SSRI+SNRI) comparison (brief, task 7:
    "adjusted for admission severity") -- fit on the pre-specified
    clinical-threshold-crossing indicator (:data:`CLINICAL_THRESHOLDS`), not
    the raw delta, because the raw delta's *mean* is dominated by rare
    extreme-outlier admissions (see :func:`run_primary_comparison`'s
    mean/median divergence, documented in the README) that a linear
    regression would inherit; logistic regression on a binary,
    already-robust-by-construction outcome does not have that problem.

    Charlson comorbidity index is the covariate (per
    :mod:`mimic4_pipeline.severity`'s recommendation: it covers the whole
    cohort, unlike the ICU-only first-day severity scores).
    """
    cohort = con.execute(
        "select hadm_id, exposure_class, is_overdose_admission from mimiciv_analysis.cohort "
        "where exposure_class in ('Tricyclic', 'SSRI', 'SNRI')"
    ).fetchdf()
    if exclude_overdose:
        cohort = cohort[~cohort["is_overdose_admission"]]
    cohort["is_tricyclic"] = (cohort["exposure_class"] == "Tricyclic").astype(int)

    severity = con.execute(
        "select hadm_id, charlson_comorbidity_index from mimiciv_analysis.severity"
    ).fetchdf()
    outcomes = con.execute("select * from mimiciv_analysis.outcomes_long").fetchdf()
    merged = outcomes.merge(cohort, on="hadm_id", how="inner").merge(severity, on="hadm_id", how="left")

    if verbose:
        label = "excluding overdose admissions" if exclude_overdose else "including overdose admissions"
        print(f"=== Severity-adjusted comparison ({label}): logistic regression, "
              f"crosses_threshold ~ is_tricyclic + charlson_comorbidity_index ===")

    rows = []
    for analyte in PRIMARY_ANALYTES:
        spec = CLINICAL_THRESHOLDS[analyte]
        sub = merged[merged["analyte"] == analyte].copy()
        sub["crosses"] = sub.apply(lambda r: _crosses_threshold(r, spec), axis=1)
        sub = sub.dropna(subset=["crosses"])
        sub["crosses"] = sub["crosses"].astype(int)

        fit = adjusted_odds_ratio(sub, "crosses", "is_tricyclic", ["charlson_comorbidity_index"])
        fit["analyte"] = analyte
        rows.append(fit)
        if verbose:
            if fit["status"] == "ok":
                print(f"  {analyte:16s} adjusted OR(Tricyclic)={fit['adjusted_or']:.3f} "
                      f"[{fit['adjusted_or_ci_low']:.3f}, {fit['adjusted_or_ci_high']:.3f}] "
                      f"p={fit['adjusted_p_value']:.4f}  (n={fit['n_treated']}+{fit['n_untreated']})")
            else:
                print(f"  {analyte:16s} SKIPPED -- {fit['status']}")

    cols = ["analyte", "n_treated", "n_untreated", "n_events_treated", "n_events_untreated",
            "adjusted_or", "adjusted_or_ci_low", "adjusted_or_ci_high", "adjusted_p_value", "status"]
    return pd.DataFrame(rows)[cols]


def feasibility_projection(con: duckdb.DuckDBPyConnection, verbose: bool = True) -> pd.DataFrame:
    """Given the demo's observed per-class exposure counts, project onto
    the full MIMIC-IV database and compute the minimum detectable effect
    size (Cohen's d) for the primary Tricyclic-vs-(SSRI+SNRI) comparison at
    alpha=0.05, power=0.80.

    **Reframed 2026-08-23**: the minimum detectable d here (~0.034) is far
    below any clinically meaningful effect (see :data:`CLINICAL_THRESHOLDS`)
    -- at full database scale, statistical significance is essentially
    guaranteed regardless of whether the effect is clinically meaningful.
    **The binding constraint on this study is therefore confounding
    control (severity adjustment, dose-ambiguous exposure classification,
    co-medications), not statistical power.** This projection is kept
    because it demonstrates that directly -- not because power is what
    still needs solving.

    Classes with zero observed admissions (Tricyclic, in the demo) cannot
    be projected by simple scaling -- 0 x any scale factor is still 0, an
    artifact of small-sample sparsity, not evidence the true rate is zero.
    Those use the Wilson-interval upper bound on a 0/n proportion instead,
    reported explicitly as an upper-bound-only projection.
    """
    n_total_admissions_demo = con.execute("select count(*) from mimiciv_hosp.admissions").fetchone()[0]
    class_counts = con.execute(
        "select exposure_class, count(*) as n from mimiciv_analysis.cohort "
        "where exposure_class is not null group by exposure_class"
    ).fetchdf()
    counts_by_class = dict(zip(class_counts["exposure_class"], class_counts["n"]))

    rows = []
    for cls in ["Tricyclic", "SSRI", "SNRI", "RIMA", "Atypical"]:
        n_demo = int(counts_by_class.get(cls, 0))
        lo_p, hi_p = _wilson_ci(n_demo, n_total_admissions_demo)
        point_p = n_demo / n_total_admissions_demo
        projected_point = round(point_p * FULL_MIMIC_IV_ADMISSIONS)
        projected_lo = round(lo_p * FULL_MIMIC_IV_ADMISSIONS)
        projected_hi = round(hi_p * FULL_MIMIC_IV_ADMISSIONS)
        rows.append(
            {
                "exposure_class": cls,
                "n_demo": n_demo,
                "n_demo_total_admissions": n_total_admissions_demo,
                "projected_n_full_db_point": projected_point if n_demo > 0 else None,
                "projected_n_full_db_low_95ci": projected_lo,
                "projected_n_full_db_high_95ci": projected_hi,
                "note": "point estimate is 0/n -- Wilson upper bound only, not a real estimate" if n_demo == 0 else "",
            }
        )
    projection_df = pd.DataFrame(rows)

    if verbose:
        print(f"\n=== Feasibility projection (scaled from {n_total_admissions_demo} demo admissions to "
              f"{FULL_MIMIC_IV_ADMISSIONS:,} full-MIMIC-IV admissions -- PhysioNet v3.1, cited in module docstring) ===")
        print(projection_df.to_string(index=False))

    # Minimum detectable effect size for the primary comparison, using the
    # projected group sizes (Tricyclic uses its Wilson-upper-bound
    # projection -- a best case, since the point estimate is degenerate)
    # and the observed outcome SD from the demo's SSRI+SNRI group as the
    # best available estimate of population variance.
    outcomes = con.execute("select * from mimiciv_analysis.outcomes_long").fetchdf()
    cohort = con.execute(
        "select hadm_id, exposure_class from mimiciv_analysis.cohort where exposure_class is not null"
    ).fetchdf()
    merged = outcomes.merge(cohort, on="hadm_id", how="inner")

    power_rows = []
    n_tricyclic_proj = projection_df.loc[projection_df["exposure_class"] == "Tricyclic", "projected_n_full_db_high_95ci"].iloc[0]
    n_ssri = projection_df.loc[projection_df["exposure_class"] == "SSRI", "projected_n_full_db_point"].iloc[0] or 0
    n_snri = projection_df.loc[projection_df["exposure_class"] == "SNRI", "projected_n_full_db_point"].iloc[0] or 0
    n_comparator_proj = n_ssri + n_snri

    for analyte in PRIMARY_ANALYTES:
        comparator_vals = merged.loc[
            (merged["exposure_class"].isin(["SSRI", "SNRI"])) & (merged["analyte"] == analyte),
            "delta_from_baseline",
        ].dropna()
        if len(comparator_vals) < 2 or n_tricyclic_proj < 2 or n_comparator_proj < 2:
            power_rows.append(
                {
                    "analyte": analyte,
                    "observed_sd_comparator": float(comparator_vals.std()) if len(comparator_vals) else None,
                    "n_tricyclic_projected": int(n_tricyclic_proj),
                    "n_comparator_projected": int(n_comparator_proj),
                    "min_detectable_cohens_d": None,
                    "min_detectable_raw_units": None,
                    "status": "insufficient_data: fewer than 2 observed comparator values or a projected group size < 2",
                }
            )
            continue

        sd = float(comparator_vals.std())
        ratio = n_comparator_proj / n_tricyclic_proj
        try:
            d = TTestIndPower().solve_power(
                effect_size=None, nobs1=n_tricyclic_proj, ratio=ratio, alpha=0.05, power=0.80, alternative="two-sided"
            )
            power_rows.append(
                {
                    "analyte": analyte,
                    "observed_sd_comparator": sd,
                    "n_tricyclic_projected": int(n_tricyclic_proj),
                    "n_comparator_projected": int(n_comparator_proj),
                    "min_detectable_cohens_d": float(d),
                    "min_detectable_raw_units": float(d) * sd,
                    "status": "ok",
                }
            )
        except Exception as exc:  # noqa: BLE001 -- power solver can fail to converge; report, don't crash
            power_rows.append(
                {
                    "analyte": analyte,
                    "observed_sd_comparator": sd,
                    "n_tricyclic_projected": int(n_tricyclic_proj),
                    "n_comparator_projected": int(n_comparator_proj),
                    "min_detectable_cohens_d": None,
                    "min_detectable_raw_units": None,
                    "status": f"power_solver_failed: {exc}",
                }
            )

    power_df = pd.DataFrame(power_rows)
    if verbose:
        print(f"\n=== Minimum detectable effect size at full-DB projected scale (alpha=0.05, power=0.80) ===")
        print("  Tricyclic group size uses its Wilson-upper-bound projection (best case, since the demo's")
        print("  observed count is 0). Outcome SD is the demo's own observed SSRI+SNRI variability.")
        print(power_df.to_string(index=False))
        min_d = power_df["min_detectable_cohens_d"].dropna()
        if len(min_d):
            print(f"\n  NOTE: min detectable d ~= {min_d.min():.3f} is far below a typical 'small' effect")
            print("  (Cohen's convention: d=0.2). At full scale, statistical significance is essentially")
            print("  guaranteed -- p<0.05 alone will NOT indicate clinical relevance. The pre-specified")
            print("  clinical thresholds (CLINICAL_THRESHOLDS) and the effect-size CIs above are what")
            print("  should drive interpretation, not p-values. The binding constraint on this study is")
            print("  confounding control, not statistical power.")

    return pd.concat([projection_df.assign(_section="class_projection"), power_df.assign(_section="power")], ignore_index=True)


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(DB_PATH), read_only=True)
    try:
        # Brief, task 3: "produce the cohort both with and without" overdose
        # admissions. Both variants of the primary comparison are run here,
        # not just the column left sitting unused.
        comparison_df = run_primary_comparison(con, exclude_overdose=False)
        comparison_df.to_csv(ANALYSIS_CSV, index=False)
        print(f"\nWrote {len(comparison_df)} rows -> {ANALYSIS_CSV}")

        print()
        comparison_excl_df = run_primary_comparison(con, exclude_overdose=True)
        excl_path = OUTPUT_DIR / "analysis_primary_comparison_excl_overdose.csv"
        comparison_excl_df.to_csv(excl_path, index=False)
        print(f"\nWrote {len(comparison_excl_df)} rows -> {excl_path}")

        print()
        adjusted_df = severity_adjusted_comparison(con, exclude_overdose=False)
        adjusted_path = OUTPUT_DIR / "analysis_severity_adjusted.csv"
        adjusted_df.to_csv(adjusted_path, index=False)
        print(f"\nWrote {len(adjusted_df)} rows -> {adjusted_path}")

        print()
        feasibility_df = feasibility_projection(con)
        feasibility_df.to_csv(FEASIBILITY_CSV, index=False)
        print(f"Wrote {len(feasibility_df)} rows -> {FEASIBILITY_CSV}")
    finally:
        con.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
