"""
Uncapped FAERS sensitivity analysis for the 2,500-report-per-drug cap.

fetch_faers_openfda.py pages through openFDA without a sort parameter, so for
the 18 of 27 drugs where the 2,500-report cap binds, the sample is the first
2,500 records in openFDA's index order, not a random sample. This script
re-derives every (drug, organ system) 2x2 table from openFDA's
meta.results.total -- the count of ALL matching reports, uncapped -- using the
same drug search (build_search) and the same MedDRA PT -> organ mapping
(ORGAN_SYSTEM_PT), then compares the uncapped RORs with faers_signals.csv.

Background note: the pooled comparator is the sum of the other drugs' totals,
so a report listing two antidepressants counts in both (an approximation);
and uncapped, the high-volume SSRIs dominate the background, whereas the cap
weighted drugs roughly equally. Differences therefore reflect background
composition as well as sampling.

Outputs:
  faers_uncapped_counts.csv  per-drug total and per-organ report counts
  faers_uncapped_ror.csv     uncapped ROR / lower bound vs. the capped values

Usage:
  python faers_uncapped.py            # query openFDA (~135 requests), then compare
  python faers_uncapped.py --no-fetch # reuse faers_uncapped_counts.csv
"""
import json
import sys
import time
import urllib.error
import urllib.request

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from faers_disproportionality import DRUG_SYNONYMS, ORGAN_SYSTEM_PT
from fetch_faers_openfda import API, API_KEY, build_search

COUNTS_CSV = "faers_uncapped_counts.csv"
ROR_CSV = "faers_uncapped_ror.csv"
ORGANS = list(ORGAN_SYSTEM_PT)
PAUSE = 0.4


def total(search):
    url = f"{API}?search={search}&limit=1"
    if API_KEY:
        url += f"&api_key={API_KEY}"
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=90) as resp:
                return json.load(resp)["meta"]["results"]["total"]
        except urllib.error.HTTPError as e:
            if e.code == 404:          # openFDA's "no matches"
                return 0
            time.sleep(3 * (attempt + 1))
        except Exception:
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"gave up on {url[:200]}")


def fetch_counts():
    rows = []
    for drug, syns in DRUG_SYNONYMS.items():
        if drug == "Clorgyline":       # never marketed, zero reports
            continue
        ds = build_search(syns[0], syns)
        row = {"drug_name": drug, "total": total(ds)}
        time.sleep(PAUSE)
        for organ, pts in ORGAN_SYSTEM_PT.items():
            rx = "(" + "+OR+".join(
                "patient.reaction.reactionmeddrapt.exact:%22" + p.upper().replace(" ", "+") + "%22"
                for p in pts) + ")"
            row[organ] = total(ds + "+AND+" + rx)
            time.sleep(PAUSE)
        print(row, flush=True)
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(COUNTS_CSV, index=False)
    return df


def compare(u):
    rows = []
    for _, r in u.iterrows():
        others = u[u.drug_name != r.drug_name]
        for o in ORGANS:
            a = r[o]
            b = r.total - a
            c = others[o].sum()
            d = others.total.sum() - c
            if a == 0:
                ror = lo = np.nan
            else:
                ror = (a / b) / (c / d)
                lo = np.exp(np.log(ror) - 1.96 * np.sqrt(1 / a + 1 / b + 1 / c + 1 / d))
            rows.append(dict(drug_name=r.drug_name, organ_system=o, n_cases_uncapped=a,
                             ror_uncapped=ror, ror_lower_uncapped=lo))
    orig = pd.read_csv("faers_signals.csv")
    df = pd.DataFrame(rows).merge(
        orig[["drug_name", "organ_system", "ror", "ror_lower", "signal", "signal_ror1"]],
        on=["drug_name", "organ_system"])
    df["signal_uncapped"] = ((df.ror_lower_uncapped > 2) & (df.n_cases_uncapped >= 3)).astype(int)
    df["signal_ror1_uncapped"] = (df.ror_uncapped > 1).astype(int)
    df.to_csv(ROR_CSV, index=False)

    ok = df.dropna(subset=["ror_uncapped", "ror"])
    print(f"\nReports (drug-report pairs), uncapped: {int(u.total.sum()):,}")
    print("Spearman, capped vs uncapped ROR:")
    for o in ORGANS:
        s = ok[ok.organ_system == o]
        rho, p = spearmanr(s.ror, s.ror_uncapped)
        print(f"  {o:10s} rho={rho:+.2f}  p={p:.4f}")
    print("\nPrimary signals (lower bound > 2, >= 3 cases):")
    print("  capped:  ", df[df.signal == 1][["drug_name", "organ_system"]].values.tolist())
    print("  uncapped:", df[df.signal_uncapped == 1][["drug_name", "organ_system"]].values.tolist())
    print(f"ROR>1 signals: capped {int(df.signal_ror1.sum())}, uncapped {int(df.signal_ror1_uncapped.sum())}")

    # class-level TCA / (SSRI+SNRI) contrast on the 16 agents shared with the MIMIC-IV exposure list
    tca = ["Imipramine", "Doxepin", "Nortriptyline", "Amitriptyline", "Clomipramine",
           "Desipramine", "Trimipramine", "Amoxapine"]
    comp = ["Duloxetine", "Venlafaxine", "Sertraline", "Escitalopram", "Paroxetine",
            "Fluvoxamine", "Citalopram", "Fluoxetine"]
    p = df.pivot(index="drug_name", columns="organ_system", values="ror_uncapped")
    print("\nUncapped class ratio, TCA / SSRI+SNRI (16 shared agents):",
          {o: round(p.loc[tca, o].mean() / p.loc[comp, o].mean(), 3) for o in ORGANS})


if __name__ == "__main__":
    counts = pd.read_csv(COUNTS_CSV) if "--no-fetch" in sys.argv else fetch_counts()
    compare(counts)
