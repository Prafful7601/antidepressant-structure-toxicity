"""
Determine the exact FAERS date range actually captured by fetch_faers_openfda.py
and dump the literal openFDA query string used for every drug.

fetch_faers_openfda.py never applies a receivedate filter and never stores
receivedate, so the effective date range is simply "whatever openFDA served
as of the day the fetch was run" - this script finds out what that range was
without re-downloading all reports: for each drug's exact search query it asks
openFDA for just the earliest and latest `receivedate`
(sort=receivedate:asc/desc, limit=1), two lightweight requests per drug.

Usage:
    python faers_date_range.py
"""

import json
import time
import urllib.request
import urllib.error
import pandas as pd

from faers_disproportionality import DRUG_SYNONYMS
from fetch_faers_openfda import build_search, API, API_KEY

PAUSE = 0.3
MAX_RETRIES = 4


def _get(url):
    last = None
    for attempt in range(MAX_RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None      # no results at all for this query
            last = e
        except Exception as e:
            last = e
        time.sleep(1.0 * (attempt + 1))
    raise last if last else RuntimeError("date-range fetch failed")


def edge_date(search, direction):
    """direction: 'asc' (earliest) or 'desc' (latest) receivedate."""
    url = f"{API}?search={search}&sort=receivedate:{direction}&limit=1"
    if API_KEY:
        url += f"&api_key={API_KEY}"
    data = _get(url)
    if not data:
        return None
    results = data.get("results", [])
    return results[0].get("receivedate") if results else None


def main():
    drugs = [(syns[0], syns) for drug, syns in DRUG_SYNONYMS.items()]
    rows = []
    queries = []
    for canonical, synonyms in drugs:
        search = build_search(canonical, synonyms)
        queries.append((canonical, f"{API}?search={search}"))
        earliest = edge_date(search, "asc")
        time.sleep(PAUSE)
        latest = edge_date(search, "desc")
        time.sleep(PAUSE)
        rows.append({"drug_name": canonical, "earliest_receivedate": earliest,
                     "latest_receivedate": latest})
        print(f"  {canonical:<15} {earliest or 'n/a':>10}  ->  {latest or 'n/a':>10}")

    out = pd.DataFrame(rows)
    out.to_csv("faers_date_range.csv", index=False)

    valid = out.dropna(subset=["earliest_receivedate", "latest_receivedate"])
    overall_min = valid["earliest_receivedate"].min()
    overall_max = valid["latest_receivedate"].max()
    print(f"\nOverall FAERS receive-date range across all queried drugs: "
          f"{overall_min} - {overall_max}")

    with open("faers_query_strings.txt", "w") as f:
        f.write("Exact openFDA query used per drug (drug/event.json?search=...)\n")
        f.write(f"Overall receive-date range retrieved: {overall_min} - {overall_max}\n")
        f.write("Pulled without a date filter (whatever openFDA served as of "
                "run date); this file pins the actual range for reproducibility.\n\n")
        for canonical, q in queries:
            f.write(f"{canonical}:\n  {q}\n\n")
    print("Wrote faers_date_range.csv and faers_query_strings.txt")


if __name__ == "__main__":
    main()
