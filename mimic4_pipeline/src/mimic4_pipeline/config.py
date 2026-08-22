"""Single source of truth for the DuckDB database path.

Every other module imports ``DB_PATH`` (or calls :func:`get_db_path`) from
here instead of hardcoding a path. To move from the MIMIC-IV Clinical
Database Demo to the full, credentialed MIMIC-IV database once PhysioNet
access is granted:

1. Build the full database with mimic-code's ``buildmimic/duckdb`` pipeline
   (see ``mimic4_pipeline/scripts/build_database.sh``), pointed at the full
   data download instead of the demo.
2. Change ``_DEFAULT_DB_PATH`` below to that file's path. That is the only
   line in this codebase that should ever need to change.

An environment variable override (``MIMIC4_DB_PATH``) is also supported for
one-off runs (e.g. CI) without editing this file.
"""

from __future__ import annotations

import os
from pathlib import Path

# Directory containing this file's package, i.e. mimic4_pipeline/src/mimic4_pipeline
_PACKAGE_DIR = Path(__file__).resolve().parent
# mimic4_pipeline/ (the project root for this sub-pipeline)
PROJECT_ROOT = _PACKAGE_DIR.parent.parent

# --------------------------------------------------------------------------
# THE ONE LINE THAT CHANGES when full MIMIC-IV access is granted:
_DEFAULT_DB_PATH = PROJECT_ROOT / "data" / "mimic4_demo.duckdb"
# --------------------------------------------------------------------------

SQL_DIR = PROJECT_ROOT / "sql"
OUTPUT_DIR = PROJECT_ROOT / "output"
DATA_DIR = PROJECT_ROOT / "data"

# Row counts published for the MIMIC-IV Clinical Database Demo v2.2
# (https://physionet.org/content/mimic-iv-demo/2.2/), used to fail loudly if
# the load produced something other than the demo. Sourced from mimic-code's
# mimic-iv/buildmimic/postgres/validate_demo.sql (MIT-LCP/mimic-code, commit
# 8bcbd190ca75670cd5281f9ead3611ae1cefb73e), which mirrors the counts in the
# demo's own documentation.
DEMO_EXPECTED_ROW_COUNTS: dict[str, int] = {
    "mimiciv_hosp.admissions": 275,
    "mimiciv_hosp.d_hcpcs": 89200,
    "mimiciv_hosp.d_icd_diagnoses": 109775,
    "mimiciv_hosp.d_icd_procedures": 85257,
    "mimiciv_hosp.d_labitems": 1622,
    "mimiciv_hosp.diagnoses_icd": 4506,
    "mimiciv_hosp.drgcodes": 454,
    "mimiciv_hosp.emar": 35835,
    "mimiciv_hosp.emar_detail": 72018,
    "mimiciv_hosp.hcpcsevents": 61,
    "mimiciv_hosp.labevents": 107727,
    "mimiciv_hosp.microbiologyevents": 2899,
    "mimiciv_hosp.omr": 2964,
    "mimiciv_hosp.patients": 100,
    "mimiciv_hosp.pharmacy": 15306,
    "mimiciv_hosp.poe": 45154,
    "mimiciv_hosp.poe_detail": 3795,
    "mimiciv_hosp.prescriptions": 18087,
    "mimiciv_hosp.procedures_icd": 722,
    "mimiciv_hosp.services": 319,
    "mimiciv_hosp.transfers": 1190,
    "mimiciv_icu.icustays": 140,
    "mimiciv_icu.d_items": 4014,
    "mimiciv_icu.chartevents": 668862,
    "mimiciv_icu.datetimeevents": 15280,
    "mimiciv_icu.inputevents": 20404,
    "mimiciv_icu.outputevents": 9362,
    "mimiciv_icu.procedureevents": 1468,
}


def get_db_path() -> Path:
    """Resolve the DuckDB database path.

    Resolution order:
      1. ``MIMIC4_DB_PATH`` environment variable, if set.
      2. ``_DEFAULT_DB_PATH`` above.
    """
    env_override = os.environ.get("MIMIC4_DB_PATH")
    if env_override:
        return Path(env_override).expanduser().resolve()
    return _DEFAULT_DB_PATH.resolve()


DB_PATH = get_db_path()
