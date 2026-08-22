"""Task 3: cohort construction.

Builds one row per ``hadm_id`` with an antidepressant exposure class
assigned, through an explicit, logged filtering sequence (a CONSORT-style
flow -- every step below prints/records how many admissions it removed or
annotated). Run standalone:

    python3 -m mimic4_pipeline.cohort

Requires ``mimic4_pipeline.build_drug_mapping`` to have already been run
(it reads ``mimiciv_analysis.antidepressant_prescriptions``, task 2's
output table) and the database to have mimic-code's derived concepts
(``mimiciv_derived.charlson``), which ``scripts/build_database.sh`` builds
by default.

Documented rules (brief, task 3; revised 2026-08-23)
------------------------------------------------------
- **Multiple exposure classes in one admission are resolved by dose-based
  reclassification, PRIMARILY -- not by exclusion.** The original rule
  (exclude any admission with >1 antidepressant class) turned out to
  systematically delete the Tricyclic exposure group: every one of the
  demo's 4 amitriptyline admissions co-occurs with an SSRI, and low-dose
  amitriptyline (like low-dose trazodone) alongside an SSRI is a normal
  inpatient prescribing pattern for pain or sleep, not genuine dual-class
  antidepressant therapy. So: a class present in an admission is dropped
  from consideration -- not the admission -- when *every* agent
  contributing to that class is one of the three dose-ambiguous agents
  (amitriptyline/doxepin/trazodone) at low dose (see
  :func:`mimic4_pipeline.dose.class_is_dose_ambiguous_low`). If exactly one
  class remains after dropping such classes, the admission is assigned to
  it (``exposure_assignment_method = "dose_reclassified"``). If more than
  one class remains (a genuine multi-class exposure -- e.g. two classes
  each backed by a non-dose-ambiguous agent, or a dose-ambiguous agent at
  standard/unknown dose), the admission is excluded, same as before.
  **The demo empirically justifies the rule, with an important caveat
  worth stating precisely rather than glossing over**: of the 26 raw
  multi-class admissions, dose-based reclassification recovers 14 as
  single-class SSRI exposure (10 Atypical+SSRI, 4 SSRI+Tricyclic -- all 4
  of the demo's amitriptyline orders are 10mg or 50mg, both below the
  75mg antidepressant-dose threshold, so every SSRI+Tricyclic admission
  resolves to SSRI, not Tricyclic). **This does NOT create a non-zero
  Tricyclic group in the demo** -- by the rule's own logic, low-dose
  amitriptyline is presumptively non-psychiatric, so those 4 admissions'
  genuine antidepressant exposure IS the SSRI. The Tricyclic arm stays at
  zero for a separate, genuine reason: the demo's 100 patients contain no
  standard-dose (>=75mg) amitriptyline order at all, single-class or
  multi-class. What the rule fixes is different and still real: it stops
  15% of raw exposed admissions (14/103) from being discarded as spurious
  "dual-class" exposure when they are not, growing the SSRI/Atypical
  comparator groups (SSRI 29->43) rather than the Tricyclic group. See
  :func:`_exposure_per_admission`.
- **Pure exclusion (the original rule) is retained as a documented
  sensitivity analysis**, not discarded -- every admission's raw class
  count (`n_classes_raw`) and dose-reclassified assignment method are both
  kept in the output table, so the exclusion-only cohort is reproducible
  by filtering on `n_classes_raw == 1` alone.
- **Low-dose use that does NOT resolve a multi-class admission (i.e. a
  single-class admission's only class is low-dose-ambiguous) is still
  flagged as a separate stratum (`dose_stratum`), never silently
  dropped** -- see :mod:`mimic4_pipeline.dose` for the thresholds and
  rationale.
- **Overdose admissions are flagged (`is_overdose_admission`), not
  removed** -- both cohorts (with and without) are one boolean-column
  filter apart on the same table, per the brief's request for both.
- **Pre-existing liver disease** is pulled from mimic-code's own Charlson
  concept (``mimiciv_derived.charlson``) rather than re-derived by hand.
"""

from __future__ import annotations

import sys
from pathlib import Path

import duckdb
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mimic4_pipeline.config import DB_PATH, OUTPUT_DIR, SQL_DIR
from mimic4_pipeline.dose import (
    DOSE_THRESHOLDS_MG_PER_DAY,
    class_is_dose_ambiguous_low,
    classify_dose_stratum,
    estimate_daily_dose_mg,
)

COHORT_CSV = OUTPUT_DIR / "cohort.csv"
CONSORT_LOG_PATH = OUTPUT_DIR / "cohort_consort_flow.txt"


def _read_sql(name: str) -> str:
    return (SQL_DIR / "cohort" / name).read_text()


class ConsortLog:
    """Accumulates the CONSORT-style flow so it can be printed and saved."""

    def __init__(self):
        self.lines: list[str] = []

    def step(self, label: str, n: int, note: str = ""):
        line = f"{n:>6}  {label}"
        if note:
            line += f"  -- {note}"
        self.lines.append(line)
        print(line)

    def text(self) -> str:
        return "\n".join(self.lines) + "\n"


def _overall_dose_stratum(agents: list[str], doses_by_agent: dict[str, list]) -> tuple[str, dict[str, str]]:
    """Worst-case (most "standard-dose-like") stratum across every
    dose-ambiguous agent present, plus the per-agent stratum dict needed by
    :func:`class_is_dose_ambiguous_low`."""
    per_agent = {
        agent: classify_dose_stratum(agent, doses_by_agent[agent])
        for agent in DOSE_THRESHOLDS_MG_PER_DAY
        if agent in agents
    }
    strata = list(per_agent.values())
    if not strata:
        overall = "not_applicable"
    elif "standard_dose" in strata:
        overall = "standard_dose"
    elif "unknown_dose" in strata and "low_dose" not in strata:
        overall = "unknown_dose"
    else:
        overall = "low_dose"
    return overall, per_agent


def _exposure_per_admission(rx: pd.DataFrame) -> pd.DataFrame:
    """One row per hadm_id: which class(es)/agent(s) it was exposed to
    (systemic-route rows only), the raw class count, the dose-reclassified
    primary exposure_class (see module docstring), and per-agent dose
    strata."""
    systemic = rx[rx["is_systemic"]].copy()
    systemic["est_daily_dose_mg"] = systemic.apply(
        lambda r: estimate_daily_dose_mg(r["dose_val_rx"], r["doses_per_24_hrs"]), axis=1
    )

    rows = []
    for hadm_id, g in systemic.groupby("hadm_id"):
        classes = sorted(g["drug_class"].unique().tolist())
        agents = sorted(g["agent"].unique().tolist())
        doses_by_agent = {a: g.loc[g["agent"] == a, "est_daily_dose_mg"].tolist() for a in agents}

        dose_stratum, dose_stratum_by_agent = _overall_dose_stratum(agents, doses_by_agent)

        agents_by_class = {
            cls: sorted(g.loc[g["drug_class"] == cls, "agent"].unique().tolist()) for cls in classes
        }
        ambiguous_low_classes = [
            cls for cls in classes
            if class_is_dose_ambiguous_low(agents_by_class[cls], dose_stratum_by_agent)
        ]

        if len(classes) == 1:
            exposure_class = classes[0]
            method = "single_class"
            reclassified_away_classes: list[str] = []
        else:
            resolved_classes = [c for c in classes if c not in ambiguous_low_classes]
            if len(resolved_classes) == 1:
                exposure_class = resolved_classes[0]
                method = "dose_reclassified"
                reclassified_away_classes = ambiguous_low_classes
            elif len(resolved_classes) == 0:
                # every class present is fully explained by low-dose,
                # dose-ambiguous agents (e.g. low-dose amitriptyline AND
                # low-dose trazodone, no other class) -- nothing to assign.
                exposure_class = None
                method = "excluded_ambiguous_low_dose_only"
                reclassified_away_classes = ambiguous_low_classes
            else:
                # genuine multi-class exposure remains even after dropping
                # low-dose-ambiguous classes -- same outcome as the
                # original pure-exclusion rule.
                exposure_class = None
                method = "excluded_multiclass"
                reclassified_away_classes = ambiguous_low_classes

        # dose_stratum reported for the FINAL assigned class's own agents
        # (post-reclassification), not the raw admission-wide worst case --
        # a reclassified SSRI assignment is "not_applicable" for dose, since
        # SSRIs have no dose ambiguity; the dropped class's dose story lives
        # in reclassified_away_classes/agents instead.
        if exposure_class is not None:
            final_class_agents = agents_by_class.get(exposure_class, [])
            final_strata = [dose_stratum_by_agent[a] for a in final_class_agents if a in dose_stratum_by_agent]
            if not final_strata:
                final_dose_stratum = "not_applicable"
            elif "standard_dose" in final_strata:
                final_dose_stratum = "standard_dose"
            elif "unknown_dose" in final_strata and "low_dose" not in final_strata:
                final_dose_stratum = "unknown_dose"
            else:
                final_dose_stratum = "low_dose"
        else:
            final_dose_stratum = dose_stratum

        rows.append(
            {
                "hadm_id": hadm_id,
                "n_classes_raw": len(classes),
                "n_agents": len(agents),
                "classes_present": "|".join(classes),
                "agents_present": "|".join(agents),
                "exposure_class": exposure_class,
                "exposure_assignment_method": method,
                "reclassified_away_classes": "|".join(reclassified_away_classes),
                "reclassified_away_agents": "|".join(
                    sorted({a for c in reclassified_away_classes for a in agents_by_class[c]})
                ),
                "dose_stratum": final_dose_stratum,
            }
        )
    return pd.DataFrame(rows)


def build_cohort(con: duckdb.DuckDBPyConnection, verbose: bool = True) -> tuple[pd.DataFrame, ConsortLog]:
    log = ConsortLog()

    admissions = con.execute(_read_sql("base_admissions.sql")).fetchdf()
    log.step("All hospital admissions", len(admissions))

    try:
        rx = con.execute(
            "select * from mimiciv_analysis.antidepressant_prescriptions"
        ).fetchdf()
    except duckdb.Error as exc:
        raise RuntimeError(
            "mimiciv_analysis.antidepressant_prescriptions not found -- run "
            "`python3 -m mimic4_pipeline.build_drug_mapping` first (task 2)."
        ) from exc

    exposure = _exposure_per_admission(rx)
    cohort = admissions.merge(exposure, on="hadm_id", how="left")
    n_exposed = cohort["n_classes_raw"].notna().sum()
    log.step("Admissions with >=1 systemic curated-antidepressant prescription", int(n_exposed))

    n_raw_multi_class = int((cohort["n_classes_raw"] > 1).sum())
    log.step(
        "  of which have >1 antidepressant class present (raw, before dose reclassification)",
        n_raw_multi_class,
    )

    n_dose_reclassified = int((cohort["exposure_assignment_method"] == "dose_reclassified").sum())
    log.step(
        "    -> resolved by dose-based reclassification (PRIMARY rule, see module docstring)",
        n_dose_reclassified,
        "e.g. low-dose amitriptyline/trazodone + SSRI -> assigned to SSRI, not excluded",
    )
    n_excluded_multiclass = int((cohort["exposure_assignment_method"] == "excluded_multiclass").sum())
    n_excluded_ambig = int((cohort["exposure_assignment_method"] == "excluded_ambiguous_low_dose_only").sum())
    log.step(
        "    -> still excluded: genuine multi-class exposure remains after reclassification",
        n_excluded_multiclass,
    )
    log.step(
        "    -> still excluded: every class present was low-dose-ambiguous-only (nothing to assign)",
        n_excluded_ambig,
    )
    log.step(
        "  [sensitivity] admissions the ORIGINAL pure-exclusion rule would drop (n_classes_raw > 1)",
        n_raw_multi_class,
        "reproduce via: cohort[cohort.n_classes_raw == 1] -- ignores dose reclassification entirely",
    )

    n_excluded_total = n_excluded_multiclass + n_excluded_ambig
    cohort = cohort[cohort["n_classes_raw"].isna() | cohort["exposure_class"].notna()].copy()
    n_single_class = int(cohort["exposure_class"].notna().sum())
    log.step("Primary exposed population (single-class after dose reclassification)", n_single_class)

    for cls in ["Tricyclic", "SSRI", "SNRI", "RIMA", "Atypical"]:
        n = int((cohort["exposure_class"] == cls).sum())
        log.step(f"  of which {cls}", n)

    dose_counts = cohort.loc[cohort["exposure_class"].notna(), "dose_stratum"].value_counts()
    for stratum, n in dose_counts.items():
        if stratum != "not_applicable":
            log.step(f"  dose_stratum = {stratum} (amitriptyline/doxepin/trazodone only)", int(n))

    # Overdose flag -- annotated, not dropped, and computed over ALL
    # admissions (not just the exposed ones) so it can be joined in
    # regardless of cohort membership.
    overdose_df = con.execute(_read_sql("overdose_admissions.sql")).fetchdf()
    overdose_hadm_ids = set(overdose_df["hadm_id"])
    cohort["is_overdose_admission"] = cohort["hadm_id"].isin(overdose_hadm_ids)
    n_overdose_in_cohort = int((cohort["is_overdose_admission"] & cohort["exposure_class"].notna()).sum())
    log.step(
        "  of the single-class-exposed, flagged as antidepressant-overdose admission",
        n_overdose_in_cohort,
        "flagged not removed -- both with/without-overdose cohorts available via this column",
    )

    # Liver disease covariate -- annotated, not used for filtering.
    liver = con.execute(_read_sql("liver_disease_covariate.sql")).fetchdf()
    cohort = cohort.merge(liver, on="hadm_id", how="left")
    n_liver = int((cohort["has_liver_disease"] == True).sum())  # noqa: E712
    log.step(
        "  of all admissions in the table, with pre-existing liver disease (Charlson)",
        n_liver,
        "covariate only, not a filter",
    )

    cohort["n_classes_raw"] = cohort["n_classes_raw"].fillna(0).astype(int)
    cohort["n_agents"] = cohort["n_agents"].fillna(0).astype(int)
    cohort["dose_stratum"] = cohort["dose_stratum"].fillna("not_exposed")
    cohort["exposure_assignment_method"] = cohort["exposure_assignment_method"].fillna("not_exposed")
    for col in ["classes_present", "agents_present", "reclassified_away_classes", "reclassified_away_agents"]:
        cohort[col] = cohort[col].fillna("")

    if verbose:
        print(f"\nFinal cohort table: {len(cohort)} rows (one per admission), "
              f"{n_single_class} single-class-exposed ({n_dose_reclassified} of them via dose "
              f"reclassification), {n_excluded_total} excluded as unresolvable multi-class.")

    return cohort, log


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(DB_PATH), read_only=False)
    try:
        cohort, log = build_cohort(con)
        cohort.to_csv(COHORT_CSV, index=False)
        CONSORT_LOG_PATH.write_text(log.text())
        print(f"\nWrote {len(cohort)} rows -> {COHORT_CSV}")
        print(f"Wrote CONSORT flow -> {CONSORT_LOG_PATH}")

        con.execute("create schema if not exists mimiciv_analysis")
        con.execute("create or replace table mimiciv_analysis.cohort as select * from cohort")
        print("Persisted -> mimiciv_analysis.cohort")
    finally:
        con.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
