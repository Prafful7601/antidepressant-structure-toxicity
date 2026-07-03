"""
Fetch FAERS adverse-event reports for the 26 antidepressants directly from the
openFDA API (https://open.fda.gov/apis/drug/event/), and write them as
faers_reports.csv with the columns faers_disproportionality.py expects:

    primaryid, drugname, pt

One row per (report x drug x reaction). We search openFDA's normalized
generic_name field so brand names/salts/misspellings are already resolved.

Because we only pull reports that mention our 26 drugs, the resulting
disproportionality analysis is a WITHIN-CLASS comparison (each antidepressant
vs. the pooled background of the other antidepressants) - which is exactly what
Objective 2 ("rank drugs by safety within class") calls for.

No API key needed for light use (openFDA allows 1000 requests/day, 240/min).
Uses only the Python standard library - no extra installs.

Usage:
    python fetch_faers_openfda.py
"""

import json
import time
import urllib.request
import urllib.error
import pandas as pd

# Single source of truth for drug names + brand/salt synonyms.
from faers_disproportionality import DRUG_SYNONYMS

API = "https://api.fda.gov/drug/event.json"

# Build (canonical_lowercase, [search_terms]) for each of our 28 drugs.
# search_terms = generic name + every brand/salt synonym, so we can hit BOTH
# the normalized openfda.generic_name field (US-approved drugs) AND the raw
# free-text patient.drug.medicinalproduct field (non-US drugs like moclobemide,
# mianserin, reboxetine, tianeptine, dothiepin that openFDA never normalized).
DRUGS = [(syns[0], syns) for drug, syns in DRUG_SYNONYMS.items()]

import os

PER_REQUEST = 100        # openFDA max per request
MAX_PER_DRUG = 2500      # reports per drug. ~25 requests/drug. Kept modest
                         # because openFDA gets MUCH slower at deep pagination
                         # (high skip values), so a smaller cap is both faster
                         # per-request and lighter overall. ~4x the original 600.
PAUSE = 0.3              # seconds between requests, to be polite to the API

# Optional free API key (https://open.fda.gov/apis/authentication/) raises the
# limit from 1,000 to 120,000 requests/day. Set it in your shell, e.g.:
#   PowerShell:  $env:OPENFDA_API_KEY = "your_key_here"
API_KEY = os.environ.get("OPENFDA_API_KEY", "").strip()


def build_search(canonical, synonyms):
    """OR-query across the normalized generic_name field and the free-text
    medicinalproduct field (one clause per synonym). Built by hand so the
    '+OR+' operators reach openFDA literally rather than percent-encoded."""
    clauses = [f'patient.drug.openfda.generic_name:%22{canonical}%22']
    for term in dict.fromkeys(synonyms):          # de-dupe, keep order
        phrase = term.replace(" ", "+")           # multi-word -> +-joined phrase
        clauses.append(f'patient.drug.medicinalproduct:%22{phrase}%22')
    return "(" + "+OR+".join(clauses) + ")"


MAX_RETRIES = 4          # retry transient network errors before giving up on a page


def _get_page(url):
    """GET one page with retries on transient errors. Returns dict or None.
    Raises urllib.error.HTTPError(404) sentinel to signal 'no more results'."""
    last = None
    for attempt in range(MAX_RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=90) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise                       # caller treats as end-of-results
            last = e
        except Exception as e:              # timeouts, dropped connections, etc.
            last = e
        time.sleep(1.5 * (attempt + 1))     # linear backoff before retry
    raise last if last else RuntimeError("unknown fetch error")


def fetch_drug(canonical, synonyms):
    """Return list of (primaryid, drugname, pt) tuples for one drug."""
    rows = []
    search = build_search(canonical, synonyms)
    for skip in range(0, MAX_PER_DRUG, PER_REQUEST):
        url = f"{API}?search={search}&limit={PER_REQUEST}&skip={skip}"
        if API_KEY:
            url += f"&api_key={API_KEY}"
        try:
            data = _get_page(url)
        except urllib.error.HTTPError as e:
            # 404 = no (more) results for this query; that's normal, stop paging.
            if e.code == 404:
                break
            print(f"  ! HTTP {e.code} for {canonical} at skip={skip} (gave up)")
            break
        except Exception as e:
            print(f"  ! error for {canonical} at skip={skip} after retries: {e}")
            break

        results = data.get("results", [])
        if not results:
            break

        for rpt in results:
            pid = rpt.get("safetyreportid")
            reactions = rpt.get("patient", {}).get("reaction", []) or []
            for rx in reactions:
                pt = rx.get("reactionmeddrapt")
                if pid and pt:
                    rows.append((pid, canonical, pt.lower().strip()))

        # Fewer than a full page means we've reached the end.
        if len(results) < PER_REQUEST:
            break
        time.sleep(PAUSE)
    return rows


def main():
    all_rows = []
    for i, (canonical, synonyms) in enumerate(DRUGS, 1):
        rows = fetch_drug(canonical, synonyms)
        n_reports = len({r[0] for r in rows})
        print(f"[{i:2d}/{len(DRUGS)}] {canonical:<15} {n_reports:>5} reports, "
              f"{len(rows):>6} reaction rows")
        all_rows.extend(rows)
        time.sleep(PAUSE)

    df = pd.DataFrame(all_rows, columns=["primaryid", "drugname", "pt"])
    df = df.drop_duplicates()
    df.to_csv("faers_reports.csv", index=False)
    print(f"\nWrote faers_reports.csv: {len(df)} rows, "
          f"{df['primaryid'].nunique()} unique reports, "
          f"{df['drugname'].nunique()} drugs.")
    print("Now run:  python faers_disproportionality.py")


if __name__ == "__main__":
    main()
