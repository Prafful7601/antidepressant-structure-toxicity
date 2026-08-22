"""Task 4: outcome extraction from ``hosp.labevents``.

Four analytes: ALT, AST, total bilirubin (hepatotoxicity), and serum sodium
(hyponatremia -- the SIADH-mediated signal associated with SSRIs/SNRIs in
particular). Itemids are resolved at runtime against ``d_labitems``, never
hardcoded (brief, task 4) -- see :func:`resolve_lab_itemids`, which prints
every candidate itemid it considered and which it picked, for review.

Two output shapes (brief: "controlled by a flag"):
  - ``baseline``: first value within 24h of admission, per analyte.
  - ``trajectory``: the clinically relevant extreme value over the whole
    admission (see :data:`ANALYTE_EXTREME_DIRECTION`) and its
    time-to-extreme, plus the delta from baseline -- because, per the
    brief, absolute enzyme values in an ICU cohort are dominated by the
    admission's primary reason, not by antidepressant exposure, so the
    *change* from the patient's own baseline is the more meaningful signal.

Both shapes are computed from the same underlying long-format table
(:func:`extract_outcomes`); ``shape`` only controls which wide CSV gets
written by :func:`main`.
"""

from __future__ import annotations

import sys
from pathlib import Path

import duckdb
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mimic4_pipeline.config import DB_PATH, OUTPUT_DIR, SQL_DIR

# For ALT/AST/total bilirubin, an elevated (max) value is the
# hepatotoxicity-relevant extreme. For sodium, a *depressed* (min) value is
# the antidepressant-relevant extreme (SIADH-mediated hyponatremia,
# classically associated with SSRIs/SNRIs) -- "peak" is interpreted as
# "peak severity in the direction that matters for this analyte", not
# literally the maximum in both cases.
ANALYTE_EXTREME_DIRECTION: dict[str, str] = {
    "ALT": "max",
    "AST": "max",
    "Total Bilirubin": "max",
    "Sodium": "min",
}

BASELINE_WINDOW_HOURS = 24


def _read_sql(subdir: str, name: str) -> str:
    return (SQL_DIR / subdir / name).read_text()


def resolve_lab_itemids(con: duckdb.DuckDBPyConnection, verbose: bool = True) -> dict[str, list[int]]:
    """Resolve ALT/AST/Total-Bilirubin/Sodium itemids from d_labitems at
    runtime, cross-checked against which candidates actually have rows in
    labevents. Prints every candidate considered and the resolution
    reasoning -- review this output (brief: "print what you resolve and
    let me confirm").
    """
    candidates = con.execute(_read_sql("outcomes", "lab_itemid_candidates.sql")).fetchdf()
    resolved: dict[str, list[int]] = {}

    if verbose:
        print("=== Lab itemid resolution (task 4) ===")

    for analyte in candidates["analyte"].unique():
        g = candidates[candidates["analyte"] == analyte].copy()
        g["n_labevents_rows"] = [
            con.execute(f"select count(*) from mimiciv_hosp.labevents where itemid = {i}").fetchone()[0]
            for i in g["itemid"]
        ]
        if verbose:
            print(f"\n{analyte}:")
            print(g[["itemid", "label", "fluid", "category", "n_labevents_rows"]].to_string(index=False))

        # Within one fluid+category combination, more than one itemid with
        # data indicates a lab-system coding migration for the *same*
        # measurement (observed for ALT/AST/bilirubin/sodium in this demo:
        # a newer 53xxx/52xxx itemid exists in the dictionary alongside the
        # populated one, with zero rows in the demo) -- union all itemids
        # with data in the highest-volume fluid+category rather than
        # picking one arbitrarily, so a future period using the other
        # itemid isn't silently dropped. A *different* fluid/category (e.g.
        # whole-blood/blood-gas sodium vs blood-chemistry sodium) is a
        # different specimen/method and is reported but excluded by
        # default.
        g["fluid_category"] = g["fluid"].fillna("?") + " / " + g["category"].fillna("?")
        totals = g.groupby("fluid_category")["n_labevents_rows"].sum().sort_values(ascending=False)
        if totals.empty or totals.iloc[0] == 0:
            resolved[analyte] = []
            if verbose:
                print(f"  -> WARNING: no candidate itemid has any data for {analyte}.")
            continue

        chosen_fc = totals.index[0]
        chosen = g[(g["fluid_category"] == chosen_fc) & (g["n_labevents_rows"] > 0)]
        resolved[analyte] = sorted(int(i) for i in chosen["itemid"])
        if verbose:
            print(f"  -> RESOLVED: itemid(s) {resolved[analyte]}  ({chosen_fc})")
            for fc in totals.index:
                if fc != chosen_fc and totals[fc] > 0:
                    print(f"     excluded alternate specimen/method: {fc} ({int(totals[fc])} rows) -- review if needed")

    return resolved


def extract_outcomes(
    con: duckdb.DuckDBPyConnection,
    resolved_itemids: dict[str, list[int]],
    hadm_ids: list[int] | None = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """Long-format table: one row per (hadm_id, analyte), with baseline
    (first value within 24h of admission) and trajectory (extreme value,
    time-to-extreme, delta from baseline) columns together.
    """
    frames = []
    for analyte, itemids in resolved_itemids.items():
        if not itemids:
            continue
        sql = _read_sql("outcomes", "labevents_for_itemids.sql").format(
            itemid_list=",".join(str(i) for i in itemids)
        )
        le = con.execute(sql).fetchdf()
        if hadm_ids is not None:
            le = le[le["hadm_id"].isin(hadm_ids)]
        if le.empty:
            if verbose:
                print(f"  {analyte}: 0 labevents rows in scope")
            continue

        le["hours_from_admit"] = (le["charttime"] - le["admittime"]).dt.total_seconds() / 3600.0
        direction = ANALYTE_EXTREME_DIRECTION[analyte]

        rows = []
        for hadm_id, g in le.groupby("hadm_id"):
            g = g.sort_values("charttime")
            baseline_window = g[(g["hours_from_admit"] >= 0) & (g["hours_from_admit"] <= BASELINE_WINDOW_HOURS)]
            baseline_row = baseline_window.iloc[0] if not baseline_window.empty else None

            extreme_idx = g["valuenum"].idxmax() if direction == "max" else g["valuenum"].idxmin()
            extreme_row = g.loc[extreme_idx]

            baseline_value = baseline_row["valuenum"] if baseline_row is not None else None
            delta = (extreme_row["valuenum"] - baseline_value) if baseline_value is not None else None

            # Extreme value expressed as a multiple of its own draw's ULN
            # (ref_range_upper) -- the standard Hy's-Law-style hepatic
            # signal metric, used downstream (task 7) to define a
            # pre-specified clinically meaningful effect rather than an
            # arbitrary raw-unit threshold. Sodium has no meaningful
            # "multiple of ULN" reading (hyponatremia is assessed on the
            # absolute scale), but the field is computed uniformly here and
            # simply not used for sodium's threshold downstream.
            extreme_uln = extreme_row.get("ref_range_upper")
            extreme_x_uln = (
                extreme_row["valuenum"] / extreme_uln if extreme_uln not in (None, 0) and not pd.isna(extreme_uln) else None
            )

            rows.append(
                {
                    "hadm_id": hadm_id,
                    "analyte": analyte,
                    "extreme_direction": direction,
                    "n_measurements": len(g),
                    "baseline_value": baseline_value,
                    "baseline_hours_from_admit": baseline_row["hours_from_admit"] if baseline_row is not None else None,
                    "extreme_value": extreme_row["valuenum"],
                    "extreme_ref_range_upper": extreme_uln,
                    "extreme_x_uln": extreme_x_uln,
                    "time_to_extreme_hours": extreme_row["hours_from_admit"],
                    "delta_from_baseline": delta,
                }
            )
        frames.append(pd.DataFrame(rows))
        if verbose:
            n_hadm = len(frames[-1])
            n_baseline = frames[-1]["baseline_value"].notna().sum()
            print(f"  {analyte}: {len(le)} labevents rows -> {n_hadm} admissions, {n_baseline} with a baseline value in the first {BASELINE_WINDOW_HOURS}h")

    if not frames:
        return pd.DataFrame(
            columns=[
                "hadm_id", "analyte", "extreme_direction", "n_measurements",
                "baseline_value", "baseline_hours_from_admit", "extreme_value",
                "extreme_ref_range_upper", "extreme_x_uln",
                "time_to_extreme_hours", "delta_from_baseline",
            ]
        )
    return pd.concat(frames, ignore_index=True)


def pivot_baseline(long_df: pd.DataFrame) -> pd.DataFrame:
    """Wide shape: one row per hadm_id, one baseline-value column per analyte."""
    if long_df.empty:
        return pd.DataFrame(columns=["hadm_id"])
    wide = long_df.pivot(index="hadm_id", columns="analyte", values="baseline_value")
    wide.columns = [f"{c.lower().replace(' ', '_')}_baseline" for c in wide.columns]
    return wide.reset_index()


def pivot_trajectory(long_df: pd.DataFrame) -> pd.DataFrame:
    """Wide shape: one row per hadm_id, extreme/time-to-extreme/delta columns per analyte."""
    if long_df.empty:
        return pd.DataFrame(columns=["hadm_id"])
    pieces = []
    for value_col, suffix in [
        ("extreme_value", "extreme"),
        ("extreme_x_uln", "extreme_x_uln"),
        ("time_to_extreme_hours", "time_to_extreme_hours"),
        ("delta_from_baseline", "delta_from_baseline"),
        ("n_measurements", "n_measurements"),
    ]:
        w = long_df.pivot(index="hadm_id", columns="analyte", values=value_col)
        w.columns = [f"{c.lower().replace(' ', '_')}_{suffix}" for c in w.columns]
        pieces.append(w)
    wide = pd.concat(pieces, axis=1)
    return wide.reset_index()


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shape", choices=["baseline", "trajectory", "both"], default="both")
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(DB_PATH), read_only=False)
    try:
        resolved = resolve_lab_itemids(con)

        try:
            cohort_hadm_ids = con.execute("select hadm_id from mimiciv_analysis.cohort").fetchdf()["hadm_id"].tolist()
        except duckdb.Error as exc:
            raise RuntimeError(
                "mimiciv_analysis.cohort not found -- run `python3 -m mimic4_pipeline.cohort` first (task 3)."
            ) from exc

        print(f"\n=== Extracting outcomes for {len(cohort_hadm_ids)} cohort admissions ===")
        long_df = extract_outcomes(con, resolved, hadm_ids=cohort_hadm_ids)

        con.execute("create schema if not exists mimiciv_analysis")
        con.execute("create or replace table mimiciv_analysis.outcomes_long as select * from long_df")
        print("Persisted -> mimiciv_analysis.outcomes_long")

        if args.shape in ("baseline", "both"):
            baseline_df = pivot_baseline(long_df)
            path = OUTPUT_DIR / "outcomes_baseline.csv"
            baseline_df.to_csv(path, index=False)
            print(f"Wrote {len(baseline_df)} rows -> {path}")

        if args.shape in ("trajectory", "both"):
            trajectory_df = pivot_trajectory(long_df)
            path = OUTPUT_DIR / "outcomes_trajectory.csv"
            trajectory_df.to_csv(path, index=False)
            print(f"Wrote {len(trajectory_df)} rows -> {path}")
    finally:
        con.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
