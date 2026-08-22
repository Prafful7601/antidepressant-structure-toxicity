#!/usr/bin/env python3
"""Fail loudly if the built database's row counts don't match the published
MIMIC-IV Clinical Database Demo v2.2 counts.

mimic-code's own build_mimic.sh runs a row-count check too, but only warns
(non-zero exit is not guaranteed) on mismatch. This script is the strict
gate: any missing table or wrong count is a hard failure (exit 1), so a
partial or corrupted load cannot silently pass as a successful build.

Run standalone at any time:
    python3 scripts/validate_row_counts.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import duckdb

from mimic4_pipeline.config import DEMO_EXPECTED_ROW_COUNTS, get_db_path


def main() -> int:
    db_path = get_db_path()
    if not db_path.exists():
        print(f"ERROR: database file not found at {db_path}", file=sys.stderr)
        print("  Run scripts/build_database.sh first.", file=sys.stderr)
        return 1

    con = duckdb.connect(str(db_path), read_only=True)
    failures = []
    print(f"Validating {db_path} against MIMIC-IV demo v2.2 published row counts:\n")
    print(f"{'table':<32} {'expected':>10} {'observed':>10}  status")
    print("-" * 68)
    for table, expected in DEMO_EXPECTED_ROW_COUNTS.items():
        try:
            observed = con.execute(f"SELECT count(*) FROM {table}").fetchone()[0]
        except duckdb.Error as exc:
            observed = None
            status = "FAILED (table missing/unreadable)"
            failures.append((table, expected, observed, str(exc)))
        else:
            if observed == expected:
                status = "OK"
            else:
                status = "FAILED"
                failures.append((table, expected, observed, None))
        print(f"{table:<32} {expected:>10} {str(observed):>10}  {status}")
    con.close()

    print()
    if failures:
        print(f"ROW COUNT VALIDATION FAILED: {len(failures)} table(s) mismatched.", file=sys.stderr)
        for table, expected, observed, err in failures:
            detail = f" ({err})" if err else ""
            print(f"  - {table}: expected {expected}, got {observed}{detail}", file=sys.stderr)
        print(
            "\nThis means the database was NOT built as the MIMIC-IV demo v2.2 "
            "(wrong data version, partial load, or a schema/loader bug). Do not "
            "proceed with analysis until this is resolved.",
            file=sys.stderr,
        )
        return 1

    print(f"All {len(DEMO_EXPECTED_ROW_COUNTS)} tables match published demo row counts.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
