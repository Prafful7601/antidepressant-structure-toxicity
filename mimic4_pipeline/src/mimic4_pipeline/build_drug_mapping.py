"""Build the antidepressant drug-exposure mapping from actual
``mimiciv_hosp.prescriptions`` data and write the reviewable audit CSV.

Run standalone:
    python3 -m mimic4_pipeline.build_drug_mapping

What this does, in order (mirrors the brief, section 2):
1. Pull every distinct ``drug`` value from prescriptions and run the
   word-boundary matcher (:mod:`mimic4_pipeline.drug_mapping`) against it.
2. Pull full prescription rows for every matched drug-text value, apply the
   systemic-route filter, and flag combination products.
3. Print, prominently and early, the per-agent prescription and admission
   counts -- this is the binding constraint on the whole study (whether
   enough exposed admissions exist), so it comes before the detailed
   review material below rather than after it -- and how many admissions
   have more than one antidepressant class (which determines whether the
   cohort step's multi-class exclusion rule actually matters).
4. Print every distinct drug value that matched, and every rejected
   near-miss / psychiatric-review candidate, for manual review.
5. Cross-check via formulary_drug_cd/gsn/ndc: look for prescriptions rows
   sharing a code with a confirmed match but a different (unmatched) drug
   text -- these are printed as "recovered via code" candidates.
6. Write ``output/drug_mapping_audit.csv`` and persist the same table into
   the DuckDB database as ``mimiciv_analysis.antidepressant_prescriptions``
   so downstream SQL (cohort, outcomes) can join against it directly.
"""

from __future__ import annotations

import sys
from pathlib import Path

import duckdb
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mimic4_pipeline.config import DB_PATH, OUTPUT_DIR, SQL_DIR
from mimic4_pipeline.drug_mapping import (
    is_systemic,
    match_drug_text,
    near_miss_report,
    psychiatric_review_candidates,
)

AUDIT_CSV = OUTPUT_DIR / "drug_mapping_audit.csv"


def _read_sql(name: str) -> str:
    return (SQL_DIR / "exposure" / name).read_text()


def _sql_list(values) -> str:
    return ",".join("'" + str(v).replace("'", "''") + "'" for v in values)


def build_mapping(con: duckdb.DuckDBPyConnection, verbose: bool = True) -> pd.DataFrame:
    # One row per (drug, formulary_drug_cd, gsn, ndc, route, form_rx) combo,
    # not per distinct drug string -- the same drug text can carry several
    # different codes/routes, and all of them need to reach the code-based
    # cross-check below.
    distinct_df = con.execute(_read_sql("distinct_prescription_drug_values.sql")).fetchdf()
    distinct_values = sorted(distinct_df["drug"].dropna().unique().tolist())

    # --- Step 1: text match ---------------------------------------------
    matched_agent = {}
    for v in distinct_values:
        m = match_drug_text(v)
        if m is not None:
            matched_agent[v] = m

    if verbose:
        print(f"Distinct (drug, code, route, form) combinations in prescriptions: {len(distinct_df)}")
        print(f"Distinct drug text values among them: {len(distinct_values)}")
        print(f"Distinct drug text values matched to a curated agent: {len(matched_agent)}")

    # --- Step 2: pull full prescription rows, apply route/form filter --
    all_matched_drug_values = list(matched_agent.keys())
    rx_sql = _read_sql("prescriptions_for_matched_drugs.sql").format(
        matched_drug_values=_sql_list(all_matched_drug_values) or "''"
    )
    rx = con.execute(rx_sql).fetchdf()

    rx["agent"] = rx["drug"].map(lambda d: matched_agent[d].agent)
    rx["drug_class"] = rx["drug"].map(lambda d: matched_agent[d].drug_class)
    rx["alias_matched"] = rx["drug"].map(lambda d: matched_agent[d].alias_matched)
    rx["alias_type"] = rx["drug"].map(lambda d: matched_agent[d].alias_type)
    rx["is_combination_product"] = rx["drug"].map(lambda d: matched_agent[d].is_combination_product)
    rx["is_systemic"] = rx.apply(lambda r: is_systemic(r["route"], r["form_rx"]), axis=1)

    # --- Step 3: the binding-constraint numbers, printed early and ------
    # prominently, before any of the detailed review material below.
    if verbose:
        n_nonsystemic = (~rx["is_systemic"]).sum()
        n_combo = rx["is_combination_product"].sum()
        counts = (
            rx.groupby(["drug_class", "agent"])
            .agg(n_prescriptions=("hadm_id", "size"), n_admissions=("hadm_id", "nunique"), n_patients=("subject_id", "nunique"))
            .sort_values(["drug_class", "agent"])
        )
        print("\n" + "=" * 78)
        print("PER-AGENT PRESCRIPTION / ADMISSION COUNTS -- the binding constraint on")
        print("whether this study can detect anything (feasibility, task 7)")
        print("=" * 78)
        print(counts.to_string() if not counts.empty else "  (no agent matched any prescription)")
        matched_classes = sorted(rx["drug_class"].unique().tolist())
        for cls in ["Tricyclic", "SSRI", "SNRI", "RIMA", "Atypical"]:
            if cls not in matched_classes:
                print(f"  NOTE: zero agents matched in class '{cls}' -- no exposed admissions for this class")
        print(f"\n  Total matched prescription rows: {len(rx)}")
        print(f"  -> Non-systemic route/form (excluded from exposure): {n_nonsystemic}")
        print(f"  -> Combination-product rows (flagged, not excluded here): {n_combo}")

        # Multi-class admissions: directly determines whether task 3's
        # "exclude admissions with >1 exposure class" rule actually bites.
        systemic_rx = rx[rx["is_systemic"]]
        classes_per_hadm = systemic_rx.groupby("hadm_id")["drug_class"].nunique()
        agents_per_hadm = systemic_rx.groupby("hadm_id")["agent"].nunique()
        n_multi_class = int((classes_per_hadm > 1).sum())
        # Same class, but more than one distinct agent within it (e.g. two
        # different SSRIs in one admission) -- a separate concern from the
        # multi-class exclusion rule, worth knowing about even for
        # admissions the multi-class rule would keep.
        n_multi_agent_same_class = int(((classes_per_hadm == 1) & (agents_per_hadm > 1)).sum())
        print(f"\n  Admissions with exactly 1 antidepressant class: {int((classes_per_hadm == 1).sum())}")
        print(f"  Admissions with >1 antidepressant class (would be excluded by the task-3 rule): {n_multi_class}")
        if n_multi_class:
            multi_hadm = classes_per_hadm[classes_per_hadm > 1].index.tolist()
            detail = (
                systemic_rx[systemic_rx["hadm_id"].isin(multi_hadm)]
                .groupby("hadm_id")["drug_class"]
                .unique()
                .apply(sorted)
            )
            for hadm_id, classes in detail.items():
                print(f"    hadm_id {hadm_id}: {classes}")
        print(
            f"  Admissions exposed to >1 agent WITHIN the same class (e.g. two SSRIs): "
            f"{n_multi_agent_same_class}"
        )
        print("=" * 78)

    # --- Step 4: detailed review material -------------------------------
    if verbose:
        print(f"\n=== ACCEPTED text matches ({len(matched_agent)} distinct drug values) ===")
        for v, m in sorted(matched_agent.items()):
            combo = " [COMBINATION PRODUCT]" if m.is_combination_product else ""
            print(f"  {v!r:45s} -> {m.agent} ({m.drug_class}, via {m.alias_type} alias {m.alias_matched!r}){combo}")

        rejected = near_miss_report(distinct_values)
        print(f"\n=== REJECTED near-misses ({len(rejected)}) ===")
        if not rejected:
            print("  (none -- no drug value contained a curated agent/brand name as an unbounded substring)")
        for v, reason in rejected:
            print(f"  {v!r:45s} -> {reason}")

        review = psychiatric_review_candidates(distinct_values)
        print(f"\n=== Psychiatric-medication review candidates ({len(review)}) ===")
        print("  (Superset of the near-misses above, plus brand-name hits and other known")
        print("  antidepressants not in the curated 28 and a generic-name-suffix heuristic.")
        print("  None of these are auto-included -- review by eye, especially once the full")
        print("  database's much larger drug-value list is run through this same pass.)")
        if not review:
            print("  (none)")
        for v, reason in review:
            print(f"  {v!r:45s} -> {reason}")

    # --- Step 5: formulary_drug_cd / gsn / ndc recovery cross-check ----
    matched_rows = distinct_df[distinct_df["drug"].isin(matched_agent)]
    formulary_prefixes = sorted({str(c)[:4] for c in matched_rows["formulary_drug_cd"].dropna()})
    gsn_values = sorted({str(g) for g in matched_rows["gsn"].dropna()})
    recovered = pd.DataFrame(columns=["drug", "formulary_drug_cd", "gsn", "ndc"])
    if formulary_prefixes or gsn_values:
        prefix_expr = " or ".join(f"formulary_drug_cd like '{p}%'" for p in formulary_prefixes) or "false"
        sql = _read_sql("prescriptions_by_code.sql").format(
            formulary_prefixes=prefix_expr,
            gsn_list=_sql_list(gsn_values) or "''",
            matched_drug_values=_sql_list(matched_agent.keys()) or "''",
        )
        recovered = con.execute(sql).fetchdf()

    if verbose:
        print(f"\n=== Code-based recovery cross-check (formulary_drug_cd/gsn shared with a confirmed match) ===")
        if recovered.empty:
            print("  (none -- no additional drug-text spelling shares a code with a confirmed match)")
        else:
            print("  WARNING: the following drug text values were NOT matched by name, but share a")
            print("  formulary_drug_cd or gsn with a confirmed antidepressant match. Review before including:")
            print(recovered.to_string(index=False))
    # This pipeline does not auto-include code-recovered rows without a
    # human confirming them (a shared code is strong but not certain
    # evidence, e.g. GSN reuse across reformulations) -- report only.

    return rx


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(DB_PATH), read_only=False)
    try:
        rx = build_mapping(con)
        rx.to_csv(AUDIT_CSV, index=False)
        print(f"\nWrote {len(rx)} rows -> {AUDIT_CSV}")

        con.execute("create schema if not exists mimiciv_analysis")
        con.execute("create or replace table mimiciv_analysis.antidepressant_prescriptions as select * from rx")
        print("Persisted -> mimiciv_analysis.antidepressant_prescriptions (in the DuckDB database)")
    finally:
        con.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
