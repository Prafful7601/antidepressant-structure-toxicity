"""
Independent structural verification of all 28 curated antidepressant SMILES
against PubChem (not just citalopram/escitalopram, the two entries already
known to have been corrected).

For each drug: query PubChem PUG-REST by name (falling back to brand/salt
synonyms from faers_disproportionality.DRUG_SYNONYMS) for its canonical
InChIKey, compute our own SMILES' InChIKey, and compare the first-14-char
connectivity skeleton (salt/stereo-agnostic - the same comparison already
used for the ChEMBL de-duplication and the citalopram fix).

Uses only the free PubChem PUG-REST endpoint, no API key needed.

Usage:
    python pubchem_verify.py
"""

import json
import time
import urllib.request
import urllib.parse
import urllib.error
import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi

from faers_disproportionality import DRUG_SYNONYMS

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PAUSE = 0.25          # PubChem asks for <=5 req/s
MAX_RETRIES = 4


def _get(url):
    last = None
    for attempt in range(MAX_RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise
            last = e
        except Exception as e:
            last = e
        time.sleep(1.0 * (attempt + 1))
    raise last if last else RuntimeError("pubchem fetch failed")


def pubchem_inchikey(name):
    """Return PubChem's canonical InChIKey for a compound name, or None."""
    q = urllib.parse.quote(name)
    url = f"{PUG}/compound/name/{q}/property/InChIKey/JSON"
    try:
        data = _get(url)
    except Exception:
        return None
    props = data.get("PropertyTable", {}).get("Properties", [])
    return props[0]["InChIKey"] if props else None


def main():
    ds = pd.read_csv("antidepressants_dataset.csv")
    rows = []
    for _, row in ds.iterrows():
        drug, smi = row["drug_name"], row["smiles"]
        mol = Chem.MolFromSmiles(smi)
        our_key = inchi.MolToInchiKey(mol) if mol is not None else None
        our_skel = our_key[:14] if our_key else None

        # Try the canonical name first, then every brand/salt synonym.
        candidates = [drug] + DRUG_SYNONYMS.get(drug, [])
        pc_key, matched_name = None, None
        for name in candidates:
            pc_key = pubchem_inchikey(name)
            if pc_key:
                matched_name = name
                break
            time.sleep(PAUSE)
        pc_skel = pc_key[:14] if pc_key else None

        match = (our_skel is not None and pc_skel is not None
                 and our_skel == pc_skel)
        rows.append({
            "drug_name": drug, "smiles": smi,
            "our_inchikey": our_key, "our_skeleton": our_skel,
            "pubchem_query_name": matched_name, "pubchem_inchikey": pc_key,
            "pubchem_skeleton": pc_skel, "skeleton_match": match,
        })
        status = "OK" if match else ("NO PUBCHEM HIT" if pc_key is None else "**MISMATCH**")
        print(f"  {drug:<15} {status}")
        time.sleep(PAUSE)

    out = pd.DataFrame(rows)
    out.to_csv("pubchem_verification.csv", index=False)
    n_match = int(out["skeleton_match"].sum())
    n_nohit = int(out["pubchem_inchikey"].isna().sum())
    n_mismatch = len(out) - n_match - n_nohit
    print(f"\n{n_match}/{len(out)} skeletons match PubChem.")
    if n_nohit:
        print(f"{n_nohit} had no PubChem name hit (need manual CID lookup): "
              f"{out.loc[out['pubchem_inchikey'].isna(), 'drug_name'].tolist()}")
    if n_mismatch:
        print(f"{n_mismatch} MISMATCHES (structure error candidates): "
              f"{out.loc[~out['skeleton_match'] & out['pubchem_inchikey'].notna(), 'drug_name'].tolist()}")
    print("\nWrote pubchem_verification.csv")


if __name__ == "__main__":
    main()
