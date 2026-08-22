#!/usr/bin/env python3
"""Day-one feasibility check.

Run this FIRST, the moment full MIMIC-IV access is granted and
``config.py``'s ``_DEFAULT_DB_PATH`` (or ``MIMIC4_DB_PATH``) points at the
full database -- before building anything else on top of it. The whole
study depends on the Tricyclic exposure arm existing at scale; this answers
that question in one read-only pass, printed to stdout, in under a minute:

    python3 scripts/day_one_counts.py

Nothing is written to the database and no upstream script needs to have
been run first. It reuses, rather than reimplements, the exact logic used
everywhere else in this pipeline:
  - task 2's drug-text matcher (``build_drug_mapping.build_mapping``)
  - task 3's dose-based exposure classification
    (``cohort._exposure_per_admission``) -- the PRIMARY rule -- alongside
    the pure-exclusion SENSITIVITY rule for direct comparison
  - task 4's lab itemid resolution (``outcomes.resolve_lab_itemids``)

Ends with a plain-language verdict on whether the Tricyclic arm is large
enough to support the planned comparison, and a fallback if not.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import duckdb
import pandas as pd
from statsmodels.stats.power import TTestIndPower

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mimic4_pipeline.build_drug_mapping import build_mapping
from mimic4_pipeline.cohort import _exposure_per_admission
from mimic4_pipeline.config import DB_PATH
from mimic4_pipeline.drug_mapping import AGENT_CLASS
from mimic4_pipeline.outcomes import resolve_lab_itemids

TRICYCLIC_AGENTS = sorted(a for a, c in AGENT_CLASS.items() if c == "Tricyclic")
DOSE_STRATA = ["low_dose", "standard_dose", "unknown_dose", "not_applicable"]
EXPOSURE_CLASSES = ["Tricyclic", "SSRI", "SNRI", "RIMA", "Atypical"]

# Cohen's convention for a "small" effect -- used here purely as a
# transparent, standard go/no-go convention for the plain-language verdict
# below (a function of sample size alone, no outcome data needed, which is
# what keeps this script fast). This is deliberately simpler than, and
# separate from, the pre-specified clinical thresholds in analysis.py
# (task 7), which are for interpreting an actual result, not gating
# whether the comparison is worth running at all.
VERDICT_EFFECT_SIZE = 0.2
VERDICT_ALPHA = 0.05
VERDICT_TARGET_POWER = 0.80


def _print_header(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def main() -> int:
    t0 = time.time()
    print(f"Database: {DB_PATH}")
    if not DB_PATH.exists():
        print(f"ERROR: no database file at {DB_PATH}.", file=sys.stderr)
        print("Build one first (scripts/build_database.sh), or point config.py /", file=sys.stderr)
        print("MIMIC4_DB_PATH at the full MIMIC-IV build.", file=sys.stderr)
        return 1

    con = duckdb.connect(str(DB_PATH), read_only=True)
    try:
        n_patients = con.execute("select count(*) from mimiciv_hosp.patients").fetchone()[0]
        n_total_admissions = con.execute("select count(*) from mimiciv_hosp.admissions").fetchone()[0]
        if n_patients == 100:
            _print_header("WARNING: this looks like the 100-patient MIMIC-IV DEMO, not the full DB")
            print("Everything below is a demo-scale sanity check, not the real answer --")
            print("re-point config.py / MIMIC4_DB_PATH at the full build and rerun.")

        rx = build_mapping(con, verbose=False)
        exposure = _exposure_per_admission(rx)
        n_exposed = len(exposure)

        resolved_itemids = resolve_lab_itemids(con, verbose=False)
        exposed_ids = exposure[["hadm_id"]]
        con.register("exposed_hadm", exposed_ids)
        lab_availability = {}
        for analyte, itemids in resolved_itemids.items():
            if not itemids:
                lab_availability[analyte] = 0
                continue
            n = con.execute(
                f"""
                select count(distinct le.hadm_id)
                from mimiciv_hosp.labevents le
                join exposed_hadm eh using (hadm_id)
                where le.itemid in ({','.join(str(i) for i in itemids)})
                  and le.valuenum is not null
                """
            ).fetchone()[0]
            lab_availability[analyte] = n
    finally:
        con.close()

    _print_header("ADMISSIONS")
    print(f"  Total admissions in this database:                {n_total_admissions:>8d}")
    print(f"  Exposed admissions (>=1 systemic curated rx):      {n_exposed:>8d}")

    _print_header("PRIMARY RULE: exposure class after dose-based reclassification")
    primary_counts = exposure["exposure_class"].value_counts(dropna=True)
    for cls in EXPOSURE_CLASSES:
        marker = "  <<< THE NUMBER THIS STUDY DEPENDS ON" if cls == "Tricyclic" else ""
        print(f"  {cls:12s} {int(primary_counts.get(cls, 0)):>8d}{marker}")
    method_counts = exposure["exposure_assignment_method"].value_counts()
    for method in ["excluded_multiclass", "excluded_ambiguous_low_dose_only"]:
        print(f"  [excluded]  {method:35s} {int(method_counts.get(method, 0)):>6d}")

    tricyclic_hadm_ids = set(exposure.loc[exposure["exposure_class"] == "Tricyclic", "hadm_id"])
    n_tricyclic = len(tricyclic_hadm_ids)

    _print_header(f"TRICYCLIC ARM DETAIL: {n_tricyclic} admissions (primary rule)")
    if n_tricyclic == 0:
        print("  (no admissions -- nothing to break down)")
    else:
        tricyclic_rx = rx[
            rx["is_systemic"] & rx["hadm_id"].isin(tricyclic_hadm_ids) & (rx["drug_class"] == "Tricyclic")
        ]
        print("  By agent (admissions containing each agent; can exceed the total above if")
        print("  more than one TCA was given in the same admission):")
        for agent in TRICYCLIC_AGENTS:
            n = tricyclic_rx.loc[tricyclic_rx["agent"] == agent, "hadm_id"].nunique()
            print(f"    {agent:16s} {n:>8d}")

        print("\n  By dose stratum (only amitriptyline/doxepin have dose ambiguity; every other")
        print("  TCA agent is 'not_applicable' -- any dose counts as genuine antidepressant use):")
        tricyclic_exposure = exposure[exposure["hadm_id"].isin(tricyclic_hadm_ids)]
        dose_counts = tricyclic_exposure["dose_stratum"].value_counts()
        for stratum in DOSE_STRATA:
            print(f"    {stratum:16s} {int(dose_counts.get(stratum, 0)):>8d}")

    _print_header("SENSITIVITY RULE (comparison only): pure exclusion, n_classes_raw == 1")
    sensitivity = exposure[exposure["n_classes_raw"] == 1]
    sensitivity_classes = sensitivity["classes_present"]
    sensitivity_counts_by_class = {}
    for cls in EXPOSURE_CLASSES:
        n = int((sensitivity_classes == cls).sum())
        sensitivity_counts_by_class[cls] = n
        print(f"  {cls:12s} {n:>8d}")

    _print_header("LAB AVAILABILITY (exposed admissions with >=1 non-null value)")
    for analyte in ["ALT", "AST", "Total Bilirubin", "Sodium"]:
        n = lab_availability.get(analyte, 0)
        pct = f"({n / n_exposed:.0%} of exposed)" if n_exposed else ""
        print(f"  {analyte:16s} {n:>8d}  {pct}")

    # --- verdict -----------------------------------------------------
    n_comparator = int(primary_counts.get("SSRI", 0)) + int(primary_counts.get("SNRI", 0))
    _print_header("VERDICT")
    print(f"  Tricyclic (primary rule): {n_tricyclic}   |   Tricyclic (sensitivity rule): "
          f"{sensitivity_counts_by_class['Tricyclic']}")
    print(f"  SSRI+SNRI comparator (primary rule): {n_comparator}")

    fallback_ranked = sorted(
        (c for c in EXPOSURE_CLASSES if c != "Tricyclic"),
        key=lambda c: primary_counts.get(c, 0),
        reverse=True,
    )
    fallback_a, fallback_b = fallback_ranked[0], fallback_ranked[1]
    fallback_text = (
        f"FALLBACK: compare {fallback_a} (n={int(primary_counts.get(fallback_a, 0))}) vs "
        f"{fallback_b} (n={int(primary_counts.get(fallback_b, 0))}) as the primary comparison "
        f"instead; report Tricyclic descriptively (case count, agent, dose-stratum breakdown "
        f"above) rather than as a formal statistical arm."
    )

    if n_tricyclic == 0:
        print("\n  NOT LARGE ENOUGH: zero Tricyclic admissions. The primary Tricyclic vs")
        print("  SSRI+SNRI comparison cannot run at all.")
        print(f"  {fallback_text}")
    elif n_tricyclic < 2 or n_comparator < 2:
        print(f"\n  NOT LARGE ENOUGH: too few admissions to assess (Tricyclic n={n_tricyclic},")
        print(f"  comparator n={n_comparator}; need >=2 in each to even estimate an effect).")
        print(f"  {fallback_text}")
    else:
        achieved_power = TTestIndPower().power(
            effect_size=VERDICT_EFFECT_SIZE,
            nobs1=n_tricyclic,
            ratio=n_comparator / n_tricyclic,
            alpha=VERDICT_ALPHA,
            alternative="two-sided",
        )
        print(f"\n  Achieved power to detect a conventional 'small' effect (Cohen's d="
              f"{VERDICT_EFFECT_SIZE}, alpha={VERDICT_ALPHA}) at these observed sample sizes: "
              f"{achieved_power:.2f}")
        if achieved_power >= VERDICT_TARGET_POWER:
            print(f"\n  LARGE ENOUGH: achieved power ({achieved_power:.2f}) meets the "
                  f"{VERDICT_TARGET_POWER:.0%} convention. The primary Tricyclic vs SSRI+SNRI")
            print("  comparison is worth running as planned.")
        elif achieved_power >= 0.5:
            print(f"\n  BORDERLINE: achieved power ({achieved_power:.2f}) is below the "
                  f"{VERDICT_TARGET_POWER:.0%} convention but not negligible.")
            print("  The primary comparison can run, but only a moderate-or-larger effect would")
            print(f"  be detectable. Consider running it AND the fallback side by side. {fallback_text}")
        else:
            print(f"\n  NOT LARGE ENOUGH: achieved power is only {achieved_power:.2f} for even a")
            print("  conventional small effect.")
            print(f"  {fallback_text}")

    print(f"\n(finished in {time.time() - t0:.1f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
