"""Task 6: severity adjustment.

Uses mimic-code's own concept scripts rather than reimplementing severity
scores (brief, task 6). All of SOFA, SAPS-II, APS-III, SIRS, OASIS, LODS,
and Charlson are already materialized into ``mimiciv_derived`` by
mimic-code's ``concepts_duckdb/duckdb.sql``, run as part of
``scripts/build_database.sh`` (task 1) -- this module ports the *specific*
scores needed onto the cohort, joins them to `hadm_id`, and documents which
concept files were used (see ``sql/severity/severity_scores.sql``).

Coverage note that matters for task 7: SOFA/SAPS-II/APS-III/SIRS/OASIS/LODS
are all first-ICU-day scores -- an admission with no ICU stay has no value
for any of them, and in a hospital-wide (not ICU-only) cohort that is a
large fraction of admissions (147 of 275 in the demo). Charlson, by
contrast, is derived from diagnosis codes and covers every admission
regardless of ICU status. Recommendation for task 7: use **Charlson as the
primary severity/comorbidity covariate** for the full cohort, with
SOFA/APS-III as a documented ICU-only sensitivity analysis on the subset
that has one -- not the reverse, which would silently restrict the whole
study to ICU stays.
"""

from __future__ import annotations

import sys
from pathlib import Path

import duckdb
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mimic4_pipeline.config import DB_PATH, OUTPUT_DIR, SQL_DIR

SEVERITY_CSV = OUTPUT_DIR / "severity.csv"


def _read_sql(name: str) -> str:
    return (SQL_DIR / "severity" / name).read_text()


def build_severity(con: duckdb.DuckDBPyConnection, verbose: bool = True) -> pd.DataFrame:
    icu_scores = con.execute(_read_sql("severity_scores.sql")).fetchdf()

    charlson = con.execute(
        "select hadm_id, charlson_comorbidity_index from mimiciv_derived.charlson"
    ).fetchdf()

    try:
        cohort_hadm_ids = con.execute("select hadm_id from mimiciv_analysis.cohort").fetchdf()["hadm_id"]
    except duckdb.Error as exc:
        raise RuntimeError(
            "mimiciv_analysis.cohort not found -- run `python3 -m mimic4_pipeline.cohort` first (task 3)."
        ) from exc

    severity = pd.DataFrame({"hadm_id": cohort_hadm_ids})
    severity = severity.merge(charlson, on="hadm_id", how="left")
    severity = severity.merge(icu_scores.drop(columns=["stay_id"]), on="hadm_id", how="left")
    severity["has_icu_stay"] = severity["first_day_sofa"].notna()

    if verbose:
        n = len(severity)
        n_icu = int(severity["has_icu_stay"].sum())
        print(f"=== Severity coverage over {n} cohort admissions ===")
        print(f"  Charlson comorbidity index: {severity['charlson_comorbidity_index'].notna().sum()} / {n} (covers the full cohort)")
        print(f"  Admissions with a first-ICU-day severity score (SOFA/SAPS-II/APS-III/SIRS/OASIS/LODS): {n_icu} / {n}")
        if n:
            print(f"    ({n_icu / n:.0%} -- these scores are ICU-only; see module docstring for the task-7 recommendation)")

    return severity


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(DB_PATH), read_only=False)
    try:
        severity = build_severity(con)
        severity.to_csv(SEVERITY_CSV, index=False)
        print(f"\nWrote {len(severity)} rows -> {SEVERITY_CSV}")

        con.execute("create schema if not exists mimiciv_analysis")
        con.execute("create or replace table mimiciv_analysis.severity as select * from severity")
        print("Persisted -> mimiciv_analysis.severity")
    finally:
        con.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
