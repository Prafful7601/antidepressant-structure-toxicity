"""
Top up faers_reports.csv for drugs that were truncated by transient openFDA
network errors during the main fetch. Re-fetches ONLY the listed drugs (now with
retry logic in fetch_drug), replaces their rows in faers_reports.csv, and saves.

Usage:
    python topup_faers.py duloxetine venlafaxine ...
    (no args -> uses the DEFAULT_DEFICIENT list below)
"""

import sys
import pandas as pd
from faers_disproportionality import DRUG_SYNONYMS
from fetch_faers_openfda import fetch_drug

# Drugs truncated in the 2500/drug run (venlafaxine failed at 0).
DEFAULT_DEFICIENT = ["venlafaxine", "amitriptyline", "reboxetine",
                     "desipramine", "paroxetine"]

# canonical (lowercase) -> synonyms, mirroring fetch_faers_openfda.DRUGS
SYNS_BY_CANON = {syns[0]: syns for syns in DRUG_SYNONYMS.values()}


def main():
    targets = [d.lower() for d in (sys.argv[1:] or DEFAULT_DEFICIENT)]
    df = pd.read_csv("faers_reports.csv", dtype=str)
    print(f"Existing: {len(df)} rows, {df['primaryid'].nunique()} reports.\n")

    new_frames = []
    for canon in targets:
        if canon not in SYNS_BY_CANON:
            print(f"  ? unknown drug '{canon}', skipping")
            continue
        before = df[df["drugname"] == canon]["primaryid"].nunique()
        rows = fetch_drug(canon, SYNS_BY_CANON[canon])
        sub = pd.DataFrame(rows, columns=["primaryid", "drugname", "pt"])
        after = sub["primaryid"].nunique()
        print(f"  {canon:<14} {before:>5} -> {after:>5} reports "
              f"({len(sub):>6} rows)")
        new_frames.append(sub)

    # Drop old rows for the re-fetched drugs, append fresh ones.
    df = df[~df["drugname"].isin(targets)]
    out = pd.concat([df] + new_frames, ignore_index=True).drop_duplicates()
    out.to_csv("faers_reports.csv", index=False)
    print(f"\nWrote faers_reports.csv: {len(out)} rows, "
          f"{out['primaryid'].nunique()} unique reports, "
          f"{out['drugname'].nunique()} drugs.")


if __name__ == "__main__":
    main()
