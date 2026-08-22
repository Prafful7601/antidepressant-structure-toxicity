# MIMIC-IV antidepressant exposure & organ-toxicity pipeline

## What this is

A retrospective analysis pipeline that tests whether antidepressant
organ-system toxicity signals -- previously identified from FDA FAERS
spontaneous adverse-event reports, in the parent repository's
FAERS/cheminformatics work -- are reproducible in *measured* clinical data,
using the MIMIC-IV critical-care database. Spontaneous-report signals
(FAERS) are subject to reporting bias, confounding by indication, and no
denominator; this pipeline builds the corresponding cohort, exposure, and
outcome extraction needed to check those same signals against actual
admissions, prescriptions, and lab trajectories.

Given 28 curated antidepressants across 5 pharmacological classes
(Tricyclic, SSRI, SNRI, RIMA, Atypical), it identifies which admissions
were exposed to which class, extracts hepatic (ALT/AST/bilirubin) and
sodium (hyponatremia) outcomes, adjusts for admission severity, and runs a
primary Tricyclic-vs-SSRI/SNRI comparison with a pre-registered feasibility
calculation for when the full database is available.

## Why a demo-only build

**Status:** built and validated against the open **MIMIC-IV Clinical
Database Demo v2.2** (100 patients, no credentialing required). Full,
credentialed MIMIC-IV access is pending with PhysioNet. The demo will not
produce a meaningful clinical finding -- 100 patients yield only a handful
of antidepressant-exposed admissions -- but every step of the pipeline runs
end-to-end against it, and the demo run itself surfaced real, reportable
findings (a curation trap in the 28-agent list, a labeling error in
MIMIC-IV's own lab dictionary, a cohort-definition rule that was deleting
real exposure -- see the task sections below). When full access is
granted, change one line (see [Moving to full
MIMIC-IV](#moving-to-full-mimic-iv-once-credentialed) below), rerun, and
run `scripts/day_one_counts.py` first to check whether the Tricyclic arm
exists at scale before building anything further on top of it.

## No patient data in this repository

**This repository contains no patient data of any kind.** Only code (SQL
and Python) and documentation are committed; `data/`, `output/`, and
`external/` (all gitignored, see `.gitignore`) hold the downloaded/built
database and every derived artifact locally, and are never pushed. The
open MIMIC-IV Clinical Database Demo (ODbL-licensed, used throughout this
README) requires no credentialing. **The full MIMIC-IV database is
separate, requires its own PhysioNet credentialing and signed data use
agreement, and is not distributed here or by this project in any form** --
running this pipeline against it requires independently completing that
process at [physionet.org](https://physionet.org/) and downloading the
data yourself (see [Moving to full
MIMIC-IV](#moving-to-full-mimic-iv-once-credentialed)).

## Status / citation

This code accompanies a forthcoming paper and is released as the
reproducible analysis pipeline behind it. It is a companion to the parent
repository's FAERS/cheminformatics structure-toxicity pipeline (same
author, same 28-agent curated list); see that repository's own README for
the FAERS-side methodology and citation.

## Setup

```bash
cd mimic4_pipeline
python3 -m venv ../.venv        # or reuse the parent repo's venv
source ../.venv/bin/activate
pip install --upgrade pip        # a stock venv's default pip (e.g. 21.2.4) cannot
                                  # do a pyproject.toml-only editable install; a
                                  # confusing "setup.py not found" error otherwise
pip install -e ".[dev]"          # adds pytest
```

## Build the database

```bash
./scripts/build_database.sh
```

This script:
1. Installs/uses the DuckDB CLI (pinned to 1.4.5, matching the `duckdb`
   Python package).
2. Clones [MIT-LCP/mimic-code](https://github.com/MIT-LCP/mimic-code)
   (pinned to commit `8bcbd190ca75670cd5281f9ead3611ae1cefb73e`) into
   `external/mimic-code/` -- **not vendored into this repo**, so it never
   gets committed; re-cloned by anyone who runs this script.
3. Downloads the MIMIC-IV Clinical Database Demo v2.2 from PhysioNet
   (open, ODbL-licensed, no login) into `data/raw/`.
4. Runs mimic-code's own `mimic-iv/buildmimic/duckdb/build_mimic.sh`
   **rather than a hand-written loader** (per the brief: use mimic-code's
   own DuckDB build script when one exists). That script creates the
   schema, loads `hosp/` and `icu/`, and derives mimic-code's "concepts"
   (severity scores, comorbidities, first-day summaries -- see task 6)
   into a `mimiciv_derived` schema, all inside one DuckDB file.
5. Runs an independent, strict row-count validation
   (`scripts/validate_row_counts.py`) against the row counts published for
   the demo (mimic-code's own validator only warns; this one exits
   non-zero on any mismatch, so a bad build cannot pass silently).

Output: `data/mimic4_demo.duckdb` (gitignored). Confirmed on 2026-08-23: all
28 published demo tables (hosp + icu) match exactly.

Re-running the script is safe/idempotent: the demo download is skipped if
already present, and mimic-code's build script resumes from wherever it
left off (it tracks completed steps inside the database itself).

## Reproduce everything

The full pipeline, in order, from a clean checkout:

```bash
cd mimic4_pipeline
python3 -m venv ../.venv && source ../.venv/bin/activate
pip install --upgrade pip && pip install -e ".[dev]"

./scripts/build_database.sh                       # task 1
python3 -m pytest tests/ -v                        # all unit + integration tests
python3 -m mimic4_pipeline.build_drug_mapping       # task 2
python3 -m mimic4_pipeline.cohort                   # task 3
python3 -m mimic4_pipeline.outcomes --shape both    # task 4
# task 5: investigated, deferred -- see that section below, nothing to run
python3 -m mimic4_pipeline.severity                 # task 6
python3 -m mimic4_pipeline.analysis                 # task 7
```

Each step reads from tables the previous step wrote into the DuckDB
database (`mimiciv_analysis.*`) and prints its own summary as it runs; run
them in this order on a first pass. Every step is also independently
re-runnable once its inputs exist (each `main()` fails with a clear error
naming the missing upstream table/command if you skip one).

## Configuration -- the one line that changes for full MIMIC-IV

**`src/mimic4_pipeline/config.py`** is the single source of truth for the
database path. Every module imports `DB_PATH` (or calls `get_db_path()`)
from there instead of hardcoding a path. See
[Moving to full MIMIC-IV](#moving-to-full-mimic-iv-once-credentialed).

## Repository layout

```
mimic4_pipeline/
  README.md
  LICENSE
  pyproject.toml
  scripts/
    build_database.sh        # task 1: download demo + build DuckDB (see above)
    validate_row_counts.py   # strict row-count gate, called by build_database.sh
    day_one_counts.py        # run FIRST on full access: Tricyclic arm by agent + dose stratum + verdict
  sql/
    exposure/                 # task 2: drug-mapping queries
    cohort/                   # task 3: cohort/poisoning/liver-disease queries
    outcomes/                 # task 4: lab itemid resolution + extraction
    cardiac/                  # task 5: investigation only (deferred, see below)
    severity/                 # task 6: ports mimic-code's severity concepts onto the cohort
  src/mimic4_pipeline/
    config.py                 # DB path + demo row-count constants
    drug_mapping.py           # task 2: the free-text matcher (see below)
    build_drug_mapping.py     # task 2: runs the matcher against the DB, writes the audit CSV
    dose.py                   # task 3: free-text dose parsing + low-dose stratification rule
    cohort.py                 # task 3: CONSORT-style cohort construction
    outcomes.py               # task 4: lab itemid resolution + baseline/trajectory extraction
    severity.py               # task 6: severity score coverage over the cohort
    analysis.py                # task 7: guarded comparison + feasibility/power projection
  tests/
    test_drug_mapping.py      # unit tests for the matcher, incl. deliberate traps
    test_dose.py              # unit tests for dose parsing + stratification
    test_cohort.py            # unit tests for dose-based multi-class reclassification
    test_outcomes.py          # integration tests against the built demo DB
    test_analysis.py          # unit tests for the guarded comparison + clinical thresholds
  data/     (gitignored)       # raw demo download + built .duckdb file
  output/   (gitignored)       # audit CSVs, figures, results
  external/ (gitignored)       # pinned mimic-code clone
```

Data and everything derived from it are gitignored (`data/`, `output/`,
`external/`, `*.duckdb`, `*.csv`) -- **only code is committed**, per the
PhysioNet DUA and the brief's own constraint. `output/drug_mapping_audit.csv`
mentioned below is a local, reviewable file; it is not committed.

## Task 2: antidepressant exposure mapping

```bash
python3 -m mimic4_pipeline.build_drug_mapping
```

Maps the 28 curated antidepressants to the actual free-text `drug` values in
`mimiciv_hosp.prescriptions`. This was treated as its own deliverable (see
the module docstring in `src/mimic4_pipeline/drug_mapping.py` for the full
design rationale):

### Curation finding worth reporting in the paper's methods

**Matching is word-boundary regex, not raw substring matching, because four
of the 28 curated agents contain another curated agent's name as a raw
substring**: `imipramine`⊂`trimipramine`, `citalopram`⊂`escitalopram`,
`venlafaxine`⊂`desvenlafaxine`, `milnacipran`⊂`levomilnacipran`. This was
discovered by exhaustively checking all 28×27 ordered pairs, not assumed.
A naive `"citalopram" in drug_text.lower()` check would have silently
misclassified every escitalopram order as citalopram -- a real,
demonstrated risk in this specific curated list, not a hypothetical.
Word-boundary matching (`\b...\b`) resolves all four correctly; this is
proven by `tests/test_drug_mapping.py::test_curated_agent_names_do_not_collide_undetected`
and `test_longer_agent_name_not_misclassified_as_embedded_shorter_agent`,
not just asserted. **This is a curation step worth its own line in the
paper's methods.**

### Brand names are review hints, never matcher input

Only generic names drive automatic matching. The compiled brand-name list
(`BRAND_NAMES` in `drug_mapping.py`) is hand-curated and not verified
against a formulary, so a brand-name hit **never** by itself classifies a
row as exposure -- it only surfaces the drug-text value in the
"psychiatric-medication review candidates" section of the printed audit
(`psychiatric_review_candidates()`), alongside two other hint sources: a
short list of real antidepressants outside the curated 28
(`OTHER_KNOWN_ANTIDEPRESSANTS_FOR_REVIEW`), and a generic-name-suffix
heuristic (`ANTIDEPRESSANT_SUFFIX_STEMS`, e.g. `-oxetine`, `-triptyline`).
All three are deliberately over-inclusive nets for a human to look at, not
a classifier. **When full MIMIC-IV access arrives, rerun this script and
read the review-candidates section over the full (much larger) drug-value
list by eye** rather than trusting any compiled list to auto-include
anything.

### Other design points

- Combination products (e.g. an amitriptyline-chlordiazepoxide product) are
  flagged (`is_combination_product`), not silently included or dropped, so
  the cohort step can make an explicit decision.
- Non-systemic routes/forms (topical, ophthalmic, otic, nasal-local,
  cream/ointment/patch) are excluded from exposure via `is_systemic()`.
  None of the 28 curated agents has an approved topical use, so any such
  route indicates miscoding or a compounded product.
- `formulary_drug_cd` / `gsn` / `ndc` are used as an independent recovery
  cross-check: after text matching, any *other* prescriptions row sharing
  one of those codes with a confirmed match is printed for manual review
  (not auto-included).

### Demo result (2026-08-23)

16 distinct drug-text values matched, all by generic name, spanning **10 of
the 28 agents** (amitriptyline, bupropion, citalopram, duloxetine,
escitalopram, fluoxetine, mirtazapine, paroxetine, sertraline, trazodone)
across **172 prescription-order rows**. Zero near-misses were rejected and
zero psychiatric-review candidates or code-recovery candidates were found
-- i.e. the text-based match was already complete and clean in the demo.
All 172 rows were on a systemic route (`PO` or `PO/NG`); none were flagged
as combination products. This 10/28 agent coverage is expected and is
exactly the small-sample limitation the brief anticipated for the demo --
not a bug in the matcher.

**Per-agent prescription/admission counts** (the binding constraint on
whether this study can detect anything -- see task 7's feasibility
calculation):

| class | agent | n_prescriptions | n_admissions | n_patients |
|---|---|---:|---:|---:|
| Atypical | bupropion | 18 | 12 | 4 |
| Atypical | mirtazapine | 12 | 9 | 8 |
| Atypical | trazodone | 64 | 45 | 29 |
| SNRI | duloxetine | 8 | 6 | 2 |
| SSRI | citalopram | 9 | 8 | 6 |
| SSRI | escitalopram | 16 | 14 | 5 |
| SSRI | fluoxetine | 20 | 15 | 4 |
| SSRI | paroxetine | 9 | 9 | 2 |
| SSRI | sertraline | 12 | 9 | 5 |
| Tricyclic | amitriptyline | 4 | 4 | 2 |

Zero agents matched in the RIMA class (moclobemide, clorgyline -- expected;
moclobemide is not US-marketed and clorgyline was never marketed at all).

**Multi-class admissions: 26 of 103 exposed admissions (~25%) have more than
one antidepressant class present**, overwhelmingly trazodone (Atypical)
co-prescribed alongside an SSRI (22 of the 26), plus 3 SSRI+Tricyclic and 1
more Atypical+SSRI. Only 2 admissions have >1 agent *within* the same class.
This is a large enough fraction that **the task-3 multi-class exclusion rule
is not a formality here -- it removes about a quarter of otherwise-exposed
admissions**, which is exactly why the cohort step documents it as an
explicit, reviewable decision rather than an incidental filter. (Trazodone's
frequent co-prescription for sleep alongside a scheduled antidepressant is
also the reason task 3 stratifies low-dose trazodone/amitriptyline/doxepin
separately rather than dropping it.)

Full row-level detail (including which prescription row matched via which
alias) is in `output/drug_mapping_audit.csv` after running the script above
(not committed -- see gitignore note above).

### Tests

Unit tests (`tests/test_drug_mapping.py`, 27 tests) cover: the 28-agent/class
composition, all 4 real substring collisions found among the curated
agents, case-insensitivity, salt-form and release-modifier suffixes,
deliberate synthetic substring traps (an agent name embedded in a longer,
unrelated fictional compound with no word boundary on either side),
combination-product flagging, confirmation that brand names alone do NOT
drive a match (only surface as review candidates), the brand-name substring
trap, the other-known-antidepressant and suffix-stem review heuristics,
non-systemic route/form exclusion, and the near-miss reporter.

```bash
python3 -m pytest tests/ -v
```

## Task 3: cohort construction

```bash
python3 -m mimic4_pipeline.cohort
```

One row per `hadm_id`, built through an explicit, logged CONSORT-style flow
(`src/mimic4_pipeline/cohort.py`; SQL in `sql/cohort/`). Documented rules:

### Multi-class admissions: dose-based reclassification is now PRIMARY (revised 2026-08-23)

The original rule (exclude any admission with >1 antidepressant class) was
replaced after review: it was throwing away real single-class exposure, not
just genuine ambiguity. **A class present in an admission is now dropped
from consideration -- not the whole admission -- when every agent behind
that class is one of the three dose-ambiguous agents
(amitriptyline/doxepin/trazodone) at low dose** (see
`mimic4_pipeline.dose.class_is_dose_ambiguous_low`). If exactly one class
remains, the admission is assigned to it
(`exposure_assignment_method = "dose_reclassified"`). If more than one
class remains (a genuine multi-class exposure, or a dose-ambiguous agent at
standard/unknown dose), the admission is still excluded -- same outcome as
the original rule, now reached only when it's actually warranted.

**Pure exclusion is retained as a sensitivity analysis**, not discarded:
every admission's pre-reclassification class count (`n_classes_raw`) is
kept in the output table, so the original cohort is reproduced by filtering
on `n_classes_raw == 1` alone, ignoring `exposure_class` and
`exposure_assignment_method` entirely.

**Demo result, and an important caveat worth stating precisely rather than
glossing over:** of 26 raw multi-class admissions, 14 are recovered by
dose-based reclassification (10 Atypical+SSRI -> SSRI, 4 SSRI+Tricyclic ->
SSRI), growing the primary exposed population from 77 to **91** and the
SSRI group from 29 to **43**. The remaining 12 stay excluded as genuine
multi-class exposure. **This does NOT create a non-zero Tricyclic group** --
it can't, by the rule's own logic: all 4 of the demo's amitriptyline orders
are 10mg or 50mg, both below the 75mg antidepressant-dose threshold, so
every SSRI+Tricyclic admission's genuine antidepressant exposure *is* the
SSRI, and correctly resolves there. **The Tricyclic arm staying at zero is
a separate, real finding: the demo's 100 patients contain no
standard-dose amitriptyline order at all**, single-class or multi-class --
not an artifact of the exclusion rule (which is why fixing the rule didn't
fix the count). What the rule does fix is real: 14 of 103 exposed
admissions (comparator-group SSRI/Atypical, not Tricyclic) were being
discarded as spurious "dual-class" exposure when they were not.

### Other rules

- **Low-dose use that does NOT resolve a multi-class admission (a
  single-class admission whose only class is low-dose-ambiguous) is still
  flagged as a separate stratum (`dose_stratum`), never dropped.** Daily
  dose is estimated from free-text `dose_val_rx` × `doses_per_24_hrs`
  (missing frequency treated as once-daily -- the conservative floor);
  thresholds (amitriptyline/doxepin <75 mg/day, trazodone <150 mg/day) and
  the dose-string parser (handles ranges like `"25-50"` by taking the upper
  bound) live in `src/mimic4_pipeline/dose.py`, unit-tested in
  `tests/test_dose.py` and `tests/test_cohort.py`. (One real bug caught by
  writing those tests: `classify_dose_stratum` filtered out `None` but not
  float `NaN` -- and pandas silently turns `None` into `NaN` inside a float
  column after a groupby, so an unparseable dose was slipping through as if
  it were a real, always-low value. Fixed and regression-tested.)
- **Overdose admissions are flagged (`is_overdose_admission`), not
  removed** -- both "with" and "without" cohorts are one boolean filter
  apart on the same table. Identified via ICD-10 T43.0/.1/.2 **"Poisoning
  by..."** subcodes (tricyclic/tetracyclic, MAOI, other/unspecified incl.
  SSRI/SNRI) plus ICD-9 969.0x, resolved from `d_icd_diagnoses` at runtime
  rather than hardcoded (`sql/cohort/poisoning_diagnosis_codes.sql`).
  Deliberately narrower than the full T43.x range: excludes T43's "Adverse
  effect of..." / "Underdosing of..." subcodes, which share the same
  numeric prefix but are not overdoses -- confirmed against the demo, which
  has exactly one T43-prefixed diagnosis (`T43505S`, an antipsychotic
  adverse-effect code) that a naive `T43%` filter would have wrongly
  counted as an antidepressant overdose.
- **Pre-existing liver disease** is pulled from mimic-code's own Charlson
  concept (`mimiciv_derived.charlson`, built by `build_database.sh`) rather
  than re-derived from ICD codes by hand (brief's task-6 instruction,
  applied here since task 3 needs it as a covariate).

**Full demo result (2026-08-23):** of 275 admissions, 12 excluded as
genuine unresolvable multi-class exposure, leaving 263 in the table; **91
single-class-exposed** (SSRI 43, Atypical 42, SNRI 6, Tricyclic **0**, RIMA
0), 14 of them via dose reclassification. 25 admissions in the low-dose
stratum, 1 standard-dose (all trazodone -- no standard-dose TCA exposure
exists in the demo, see above). Zero overdose admissions (as expected for
100 patients). 33 of 263 admissions carry pre-existing liver disease per
Charlson. Full flow: `output/cohort_consort_flow.txt` after running the
script; table: `output/cohort.csv` (gitignored, not committed).

## Task 4: outcome extraction

```bash
python3 -m mimic4_pipeline.outcomes --shape both   # or: baseline | trajectory
```

ALT, AST, total bilirubin, and serum sodium from `hosp.labevents` joined to
`hosp.d_labitems` (`src/mimic4_pipeline/outcomes.py`; SQL in
`sql/outcomes/`). Itemids are **resolved at runtime**, never hardcoded
(`resolve_lab_itemids()`), and every candidate considered is printed:

- All 4 analytes had **more than one itemid in the dictionary** (e.g. ALT:
  50861 *and* 53084, both `Blood`/`Chemistry`, same label) -- a lab-system
  coding migration where only one member of each pair actually has data in
  the demo. The resolver unions same-fluid-same-category itemids that have
  data rather than picking one arbitrarily, so a full-database period using
  the other itemid isn't silently dropped.
- **Sodium has a genuine specimen/method duplicate**, exactly as the brief
  anticipated: blood-chemistry sodium (itemid 50983, 3007 rows) vs
  whole-blood/blood-gas sodium (50824, 141 rows). The resolver picks the
  higher-volume blood-chemistry measurement and reports the excluded
  blood-gas alternative rather than silently merging them.
- The AST label is **misspelled "Asparate Aminotransferase"** in MIMIC-IV's
  own dictionary (alongside a correctly-spelled, unpopulated duplicate) --
  found by querying the real dictionary, not assumed from memory.
- `valuenum` (MIMIC's own numeric-parsed version of the free-text `value`
  field) is filtered to non-null, which is exactly the brief's "handle
  non-numeric strings in value" instruction -- no extra parsing needed.

Two shapes, both derived from one long-format table
(`mimiciv_analysis.outcomes_long`) so trajectory's delta can reference
baseline directly:
- **Baseline**: first value within 24h of admission
  (`output/outcomes_baseline.csv`).
- **Trajectory**: the *clinically relevant* extreme value over the
  admission and its time-to-extreme, plus delta from baseline
  (`output/outcomes_trajectory.csv`). For ALT/AST/bilirubin that's the
  maximum (hepatotoxicity); **for sodium it's the minimum**, since
  SSRI/SNRI-associated hyponatremia (SIADH) is a drop, not a rise --
  "peak" is interpreted as peak severity in the direction that matters for
  each analyte, documented in `ANALYTE_EXTREME_DIRECTION`. Each row also
  carries `extreme_x_uln`: the extreme value expressed as a multiple of its
  own draw's upper limit of normal (`ref_range_upper`, confirmed
  100%-populated and consistent for all 4 resolved itemids in the demo) --
  the standard Hy's-Law-style way of expressing a hepatic signal, used by
  task 7's pre-specified clinical thresholds instead of an arbitrary
  raw-unit cutoff.

**Demo result (2026-08-23):** of 263 cohort admissions, ALT/AST/bilirubin
each resolve for ~139 admissions (~108 with a baseline value in the first
24h); sodium resolves for 238 admissions (219 with baseline) since it's
drawn far more routinely than a liver panel.

Tests: `tests/test_outcomes.py` (integration, against the built demo DB;
skips cleanly if the DB hasn't been built) checks the itemid resolution,
the 24h baseline window, delta-from-baseline sign/nullness, the
extreme-value-over-ULN computation, and the one-row-per-`hadm_id` pivot
shapes.

## Task 5: cardiac conduction (QTc) -- investigation, not yet built

**Recommendation: not feasible from core MIMIC-IV (`hosp`+`icu`) alone. A
real QTc endpoint requires MIMIC-IV-ECG -- and one premise in the brief
about that module is factually wrong (see below). This is a decision point
for you before anything downstream gets built on it.**

### What's in core MIMIC-IV

Searched `icu.d_items` for any QT/QTc/QRS/PR-interval/conduction label
(`sql/cardiac/search_d_items_conduction.sql`). Exactly one hit: itemid
`224359`, label `"QTc"`, category `"Routine Vital Signs"`, unit `sec`. No
QRS duration, no PR interval, anywhere in `d_items`. (`mimiciv_derived.rhythm`,
mimic-code's own concept, gives categorical heart rhythm/ectopy type --
useful for an arrhythmia endpoint, but not QTc.)

Checked how populated it is (`sql/cardiac/qtc_chartevents_sparsity.sql`).
**In the demo: 14 rows total, across only 4 patients / 4 admissions**, out
of 668,862 `chartevents` rows and 140 icustays. Values (0.39-0.50 sec, i.e.
390-500 ms) are physiologically plausible, so this is genuine QTc data, not
a miscoded field -- it's just a manually-charted nursing vital-sign entry,
not a systematic per-ECG measurement, and far too sparse to build an
endpoint on. (Whether it stays this sparse as a *fraction* of admissions at
full scale is a reasonable extrapolation -- same charting system, same
practice pattern -- but not something confirmable without the full
database; the row *count* would obviously be larger.)

### Correction to the brief

The brief states MIMIC-IV-ECG "is a separate credentialed module I have
not applied for." **This is not correct as of the current release.**
Checked directly against PhysioNet (2026-08-23):

| Version | Released | Access |
|---|---|---|
| v0.1 | 2022-12-23 | Credentialed |
| v0.2 | 2023-02-08 | Credentialed |
| v0.3 | 2023-07-21 | Credentialed |
| **v1.0** | **2023-09-15** | **Open (ODbL v1.0)** |

v1.0 dropped the free-text cardiologist-report field that had required
credentialing and became fully open -- *"Anyone can access the files, as
long as they conform to the terms of the specified license"* -- the same
open-access class as the MIMIC-IV demo this pipeline already uses, not the
credentialed DUA that gates full MIMIC-IV. No separate PhysioNet
application is needed. **When this gets built, pin explicitly to v1.0**
(`https://physionet.org/content/mimic-iv-ecg/1.0/`), not an unversioned
URL, since earlier versions are a different (credentialed, report-bearing)
dataset.

### What MIMIC-IV-ECG actually offers, and its real limitations

`machine_measurements.csv` (~174 MB, ~800,000 ECGs from ~160,000 subjects)
does **not** provide a precomputed QTc column. Its data dictionary lists
per-ECG wave onset/end timestamps in ms: `rr_interval`, `p_onset`, `p_end`,
`qrs_onset`, `qrs_end`, `t_end` (plus axes and machine-generated text
reports). QT interval, PR interval, QRS duration, and QTc would all have to
be **derived** here, not read off a column:
- `QT = t_end - qrs_onset`
- `PR = qrs_onset - p_onset`
- `QRS duration = qrs_end - qrs_onset`
- `QTc` = QT corrected by `rr_interval` via a chosen formula (Bazett is
  simplest/most common; Fridericia is a common, more accurate alternative
  at heart-rate extremes) -- **this choice would need to be a documented
  methodological decision**, not a silent default.

Two real limitations beyond that:
1. **Links to MIMIC-IV by `subject_id` only -- no `hadm_id`.** Matching an
   ECG to a specific hospital admission means joining `subject_id` +
   `ecg_time` falling within `[admittime, dischtime]`, which is doable but
   imperfect and needs its own documented rule and reported match rate.
2. Per PhysioNet's own description, **only ~55% of these ECGs overlap with
   a hospital admission at all** (~25% overlap an ED visit instead) --
   so even before any antidepressant-exposure filtering, a large fraction
   of the ~800k ECGs are simply not usable for an admission-level analysis.

### Recommendation and decision

Given the above, a QTc endpoint is buildable, but it is a materially larger
piece of work than tasks 1-4: a second dataset to download, a derivation
step with a methodological choice to defend, and an imperfect
subject+time-window match to `hadm_id` with its own documented loss rate.

**Decision (2026-08-23): deferred.** Proceeding with tasks 6-7 without a
cardiac arm; returning to MIMIC-IV-ECG once the rest of the pipeline is
validated end-to-end. Requirements for when that happens, so this isn't
re-litigated from scratch:

1. **Don't default to Bazett alone.** It over-corrects at high heart rates,
   and ICU patients are frequently tachycardic. Compute Fridericia as well
   and report both; make the primary correction a config flag (same
   one-flag-changes-it pattern as `config.py`'s `DB_PATH`), not a hardcoded
   choice.
2. **Check whether `machine_measurements.csv` already includes a QTc field
   before deriving one from the onset/end timestamps.** The columns
   surfaced during this investigation (`rr_interval`, `p_onset`, `p_end`,
   `qrs_onset`, `qrs_end`, `t_end`, axes) came from a webpage summary of the
   data dictionary, not a full read of the actual
   `machine_measurements_data_dictionary.csv` -- confirm directly against
   that file once downloaded. If a machine-computed QTc exists, prefer it
   and compare against the derived value rather than assuming which is
   better.
3. **Report the ECG-to-admission match rate explicitly**, plus how many
   *exposed* admissions have at least one matched ECG. If that number is
   small, the cardiac arm dies right there and that needs to surface early,
   not after the rest of the analysis is built on top of it.
4. **QT-prolonging co-medications are the dominant confounder and must be
   covariates, not an afterthought**: antipsychotics, ondansetron,
   methadone, macrolides, and fluoroquinolones are all common in an ICU
   population, plus potassium/magnesium/calcium derangement. Without these,
   a cardiac result is uninterpretable regardless of how carefully QTc
   itself is derived.

## Task 6: severity adjustment

```bash
python3 -m mimic4_pipeline.severity
```

Uses mimic-code's own concept scripts rather than reimplementing severity
scores (brief). Everything needed (SOFA, SAPS-II, APS-III, SIRS, OASIS,
LODS, Charlson) is **already materialized** into `mimiciv_derived` by
mimic-code's `concepts_duckdb/duckdb.sql`, run as part of
`scripts/build_database.sh` in task 1 -- task 6's job was porting the
specific scores onto the cohort and citing which concept files were used
(`sql/severity/severity_scores.sql` cites each one: `firstday/first_day_sofa.sql`,
`score/sapsii.sql`, `score/apsiii.sql`, `score/sirs.sql`, `score/oasis.sql`,
`score/lods.sql`, `comorbidity/charlson.sql`).

**Coverage finding that matters for task 7:** SOFA/SAPS-II/APS-III/SIRS/
OASIS/LODS are all first-**ICU-day** scores -- an admission with no ICU
stay has none of them. In a hospital-wide (not ICU-only) cohort that's a
large fraction: **only 124 of 263 cohort admissions (47%) have an ICU stay
at all** in the demo. Charlson, derived from diagnosis codes rather than
ICU-day physiology, covers the full cohort (263/263). 9 of the demo's 275
admissions have more than one ICU stay; severity is joined via the
*first* stay per admission (`sql/severity/severity_scores.sql`).
**Recommendation for task 7: Charlson is the primary severity/comorbidity
covariate for the full cohort; SOFA/APS-III are a documented ICU-only
sensitivity analysis on the subset that has one** -- not the reverse, which
would silently restrict the whole study to ICU stays.

## Task 7: analysis skeleton

```bash
python3 -m mimic4_pipeline.analysis
```

Primary comparison: Tricyclic-exposed vs (SSRI+SNRI)-exposed admissions
(using `mimiciv_analysis.cohort.exposure_class`, the dose-reclassified
primary assignment from task 3, not the raw exclusion-only cohort), on each
outcome analyte's delta-from-baseline (task 4). **Every statistical step is
guarded** (`safe_two_group_comparison()`, unit-tested in
`tests/test_analysis.py`) against empty groups, single-observation groups,
NaN/None values, and zero-variance groups -- each returns an explicit
`status` string and nothing else populated, rather than crashing or
fabricating a result. This is exercised for real on the demo, not just in
unit tests: **the Tricyclic group is still empty even after task 3's dose
reclassification** (see task 3's section above for why that's a genuine
finding, not a bug), so the primary comparison's group A is empty for every
analyte, and the script's own output shows exactly that -- `SKIPPED --
insufficient_data: group a is empty (n=0)` for all 4 analytes -- runs
cleanly to completion regardless.

### Effect sizes with CIs are primary; p-values are secondary (revised 2026-08-23)

The feasibility calculation below found a minimum detectable effect of
Cohen's d ≈ 0.034 at full database scale -- roughly *six times smaller*
than Cohen's own convention for a "small" effect (d=0.2). **That means
statistical significance is essentially guaranteed at full scale regardless
of whether the effect is clinically meaningful, so p-values alone would not
be an interpretable result there.** Two changes follow:

1. **`safe_two_group_comparison()` now reports the effect size itself --
   mean difference in delta-from-baseline with a 95% CI (Welch's t,
   unequal-variance) -- as the primary result.** A Mann-Whitney U p-value
   is still computed and reported, but labeled and treated as secondary: a
   distribution-robustness check on the effect size, not the headline
   number.
2. **A clinically meaningful effect threshold is pre-specified per analyte,
   before looking at any comparison result**, so "clinically meaningful"
   can't be dredged post hoc from whatever the full database happens to
   show (`CLINICAL_THRESHOLDS` in `analysis.py`, tested in
   `tests/test_analysis.py`):

   | Analyte | Threshold | Basis |
   |---|---|---|
   | ALT | ≥3× ULN | Hy's Law hepatocellular injury criterion |
   | AST | ≥3× ULN | Hy's Law hepatocellular injury criterion |
   | Total Bilirubin | ≥2× ULN | Hy's Law hyperbilirubinemia criterion |
   | Sodium | ≤130 mEq/L | Standard moderate-hyponatremia threshold (no ULN-multiple convention applies to sodium) |

   The ULN multiples use each lab draw's *own* `ref_range_upper` (task 4's
   `extreme_x_uln`), not a fixed external constant -- confirmed
   100%-populated in the demo. `clinical_threshold_crossing()` reports, per
   exposure group, the proportion of admissions crossing this threshold
   with a Wilson 95% CI (`output/analysis_clinical_threshold_crossing.csv`)
   -- e.g. in the demo, 7.7% of SSRI+SNRI-exposed admissions with an
   evaluable ALT cross the ≥3×ULN threshold (95% CI 2.1-24.1%; Tricyclic is
   empty, same reason as above).

### Feasibility/power calculation

Given the primary comparison can't run on the demo, this projection is
still the most substantive quantitative output of this stage -- but its
conclusion is now the opposite of "we need more power":

- Per-class exposure counts are projected from the demo onto the full
  MIMIC-IV database's published admission count (**546,028**, PhysioNet
  v3.1, cited and dated in `analysis.py`'s docstring -- an external,
  cited assumption, not derived from the demo; re-run this projection
  against the real count once the full database is built).
- **Classes with zero observed admissions (Tricyclic, RIMA) cannot be
  projected by simple scaling** -- 0 x any scale factor is still 0, an
  artifact of small-sample sparsity, not evidence of a true zero rate.
  Those use a Wilson-interval upper bound on the 0/275 proportion instead
  (**≤7,522** admissions, 95% CI), reported explicitly as an upper-bound
  projection, not a point estimate.
- **Even using that conservative upper-bound floor for the Tricyclic
  group, the projected full-database sample (7,522 vs ~97,000 SSRI+SNRI)
  would detect a minimum Cohen's d ≈ 0.034** at alpha=0.05/power=0.80 --
  translating to about 2.2 U/L for ALT, 1.3 U/L for AST, 0.016 mg/dL for
  bilirubin, and 0.25 mEq/L for sodium (using the demo's own observed
  SSRI+SNRI outcome variance as the best available estimate).
- **Conclusion: the binding constraint on this study is confounding
  control, not statistical power.** A d of 0.034 is trivially detectable --
  the full database will find *some* statistically significant difference
  on almost any outcome almost regardless of whether antidepressant class
  is actually driving it. What determines whether that difference means
  anything is (a) whether it clears the pre-specified clinical thresholds
  above, and (b) how well severity, dose-ambiguous exposure, and
  confounding co-medications are controlled -- not sample size. This is
  exactly why the effect-size-with-CI and clinical-threshold-crossing
  outputs above are the ones to read, not a bare p<0.05.

Outputs: `output/analysis_primary_comparison.csv`,
`output/analysis_clinical_threshold_crossing.csv`,
`output/feasibility_projection.csv`.

## Moving to full MIMIC-IV (once credentialed)

1. Download the full MIMIC-IV `hosp/` and `icu/` data from PhysioNet (this
   repo's `build_database.sh` only fetches the open demo).
2. Build it into DuckDB with the same mimic-code pipeline this project
   already uses:
   ```bash
   cd external/mimic-code/mimic-iv/buildmimic/duckdb
   ./build_mimic.sh /path/to/full/mimic-iv/data /path/to/mimic4_full.duckdb
   ```
3. In `src/mimic4_pipeline/config.py`, change `_DEFAULT_DB_PATH` to point at
   that file. **That is the only line in this codebase that should need to
   change.** (Or set the `MIMIC4_DB_PATH` environment variable instead of
   editing the file, for a one-off run.)
4. Rerun `scripts/validate_row_counts.py` -- it will need updating first,
   since `DEMO_EXPECTED_ROW_COUNTS` in `config.py` is specific to the
   100-patient demo; the full database has no equivalent fixed row-count
   oracle (its size is documented as approximate/changing across
   PhysioNet releases), so that gate is demo-only and should be replaced
   with a sanity check (e.g. patient count in the expected order of
   magnitude) for the full build.
5. **Before building anything further, run the day-one feasibility check:**
   ```bash
   python3 scripts/day_one_counts.py
   ```
   Answers the one question the whole study depends on -- does the
   Tricyclic exposure arm exist at scale -- in one read-only pass, printed
   to stdout, in under a minute (0.9s on the demo), with nothing else
   needing to have been run first. Reuses task 2's matcher, task 3's
   dose-based classification, and task 4's lab-itemid resolution exactly
   (no reimplementation, no risk of drift). Prints: total and exposed
   admission counts; per-class counts under the primary (dose-reclassified)
   rule; the Tricyclic arm broken down by individual agent and by dose
   stratum, not just a single headline count; the same counts under the
   pure-exclusion sensitivity rule for comparison; how many exposed
   admissions have at least one ALT/AST/bilirubin/sodium value; and a
   plain-language verdict (using achieved statistical power at the observed
   sample sizes for a conventional small effect, Cohen's d=0.2) on whether
   the Tricyclic arm is large enough to support the primary comparison,
   naming a data-driven fallback pairwise comparison if not. Demo run
   (2026-08-23, included for reference -- the script itself warns when it's
   looking at the demo rather than the full DB): 0 Tricyclic admissions
   either way, exactly matching `mimic4_pipeline.cohort`'s own numbers (see
   task 3 above) -- a consistency check on the script itself, not a new
   finding; verdict on the demo is "NOT LARGE ENOUGH" with fallback SSRI
   (n=43) vs Atypical (n=42).

## Provenance

- mimic-code: [MIT-LCP/mimic-code](https://github.com/MIT-LCP/mimic-code),
  commit `8bcbd190ca75670cd5281f9ead3611ae1cefb73e`. Used for its
  `mimic-iv/buildmimic/duckdb` build pipeline (task 1) and its
  `mimic-iv/concepts_duckdb` severity/comorbidity concept views (task 6).
  MIT-LCP/mimic-code is MIT-licensed.
- MIMIC-IV Clinical Database Demo v2.2:
  [physionet.org/content/mimic-iv-demo/2.2](https://physionet.org/content/mimic-iv-demo/2.2/),
  ODbL-licensed, no credentialing required.

## License

This project's code is MIT-licensed (see `LICENSE`). It does not cover, and
does not grant any rights to, MIMIC-IV data itself -- see [No patient
data](#no-patient-data-in-this-repository) above.
