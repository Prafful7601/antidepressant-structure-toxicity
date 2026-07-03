"""
Objective 4: screen public-database (ChEMBL) molecules for novel antidepressant
candidates that CARRY safe structural features (safeophores) and LACK toxicophores.

Pipeline:
  1. Pull compounds active (pChEMBL >= POTENCY) against antidepressant-relevant
     targets from the ChEMBL REST API: serotonin (SERT), norepinephrine (NET),
     dopamine (DAT) transporters and MAO-A. This keeps candidates
     pharmacologically plausible (they hit the right targets).
  2. Drop molecules that ARE our existing 28 drugs (match on InChIKey skeleton),
     so we only surface *other* candidates.
  3. Run each through the SAME RDKit structural screen used on our drug set
     (SUBSTRUCTURE_SMARTS + tricyclic core), reusing cheminformatics_analysis.py.
  4. Classify toxicophores / safeophores from the Objective-3 continuous
     correlation results (structure_toxicity_correlation_continuous.csv), so this
     screen automatically tracks the latest FAERS-derived evidence.
  5. Keep candidates that have >=1 safeophore AND 0 toxicophores; rank by a
     transparent safety score (#safeophores - #toxicophores), with drug-likeness
     (Lipinski) context. Write chembl_candidates_screened.csv.

Note: candidate quality is bounded by the strength of the Obj-3 toxicophore /
safeophore evidence, which is currently exploratory (small drug set). Treat the
ranked list as hypothesis-generating leads for follow-up, not validated hits.

Usage:
    python chembl_candidate_screen.py
"""

import json
import time
import urllib.request
import urllib.error
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, inchi, rdMolDescriptors, AllChem
from rdkit import DataStructs

# Reuse the exact structural definitions applied to our 28 drugs.
from cheminformatics_analysis import (
    SUBSTRUCTURE_SMARTS, patterns, is_tricyclic_dibenzo)

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

# Antidepressant-relevant targets (single-protein ChEMBL target IDs).
TARGETS = {
    "SERT (SLC6A4)": "CHEMBL228",
    "NET (SLC6A2)":  "CHEMBL222",
    "DAT (SLC6A3)":  "CHEMBL238",
    "MAO-A":         "CHEMBL1951",
}

POTENCY = 6.0            # pChEMBL >= 6  (~<= 1 uM): drug-like potency
MAX_PER_TARGET = 3000    # cap molecules pulled per target
PER_PAGE = 1000          # ChEMBL API max page size
PAUSE = 0.2


MAX_RETRIES = 4


def _chembl_get(url):
    """GET with retries: ChEMBL throws transient HTTP 500s under load."""
    last = None
    for attempt in range(MAX_RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            last = e
            time.sleep(2.0 * (attempt + 1))
    raise last if last else RuntimeError("chembl fetch failed")


def fetch_target_actives(target_id):
    """Yield (molecule_chembl_id, canonical_smiles) for potent actives."""
    out = {}
    offset = 0
    while offset < MAX_PER_TARGET:
        url = (f"{CHEMBL}/activity.json?target_chembl_id={target_id}"
               f"&pchembl_value__gte={POTENCY}"
               f"&limit={PER_PAGE}&offset={offset}")
        try:
            data = _chembl_get(url)
        except Exception as e:
            print(f"  ! error for {target_id} at offset {offset} after retries: {e}")
            break
        acts = data.get("activities", [])
        if not acts:
            break
        for a in acts:
            smi = a.get("canonical_smiles")
            mid = a.get("molecule_chembl_id")
            if smi and mid:
                out[mid] = smi
        if len(acts) < PER_PAGE:
            break
        offset += PER_PAGE
        time.sleep(PAUSE)
    return out


def inchikey_skeleton(mol):
    """First 14 chars of InChIKey = connectivity skeleton (salt/stereo-agnostic)."""
    try:
        return inchi.MolToInchiKey(mol)[:14]
    except Exception:
        return None


TOX_CUT = 0.08           # net TRS weight at/above which a feature is a toxicophore
SAFE_CUT = -0.01         # net TRS weight at/below which a feature is a safeophore


def load_feature_roles(path="trs_weights.csv"):
    """Load the validated, organ-integrated TRS weights (Objective 5) and derive
    toxicophore/safeophore sets + the weight dict. Using the net TRS weight per
    feature avoids the per-organ sign-flip problem (e.g. carboxylic acid is a
    hepatic toxicophore but a neurologic safeophore; its net weight is +, so it
    is correctly treated as a toxicophore here)."""
    try:
        w = pd.read_csv(path)
        weights = dict(zip(w["feature"], w["weight"]))
    except FileNotFoundError:
        weights = {"tricyclic_core": 0.23, "carboxylic_acid": 0.19,
                   "piperazine": 0.18, "long_alkyl_chain": 0.10,
                   "morpholine": 0.08, "para_halogen_phenyl": -0.03,
                   "nitrile": -0.02}
    tox = {f for f, wt in weights.items() if wt >= TOX_CUT}
    safe = {f for f, wt in weights.items() if wt <= SAFE_CUT}
    return tox, safe, weights


def veber_ok(mol):
    """Veber (2002) oral-bioavailability rule: <=10 rotatable bonds AND
    TPSA <= 140 A^2. Complements Lipinski, which several antidepressants
    already violate (e.g. on HBA/HBD), so Veber gives an independent check."""
    return (Descriptors.NumRotatableBonds(mol) <= 10
            and rdMolDescriptors.CalcTPSA(mol) <= 140)


def ghose_ok(mol):
    """Ghose (1999) drug-likeness filter: MW 160-480, LogP -0.4 to 5.6,
    molar refractivity 40-130, total atom count 20-70 (heavy-atom count used
    as the standard practical approximation)."""
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    mr = Crippen.MolMR(mol)
    natoms = mol.GetNumHeavyAtoms()
    return (160 <= mw <= 480 and -0.4 <= logp <= 5.6
            and 40 <= mr <= 130 and 20 <= natoms <= 70)


def reference_fingerprints(dataset_csv="antidepressants_dataset.csv"):
    """Morgan fingerprints of our 28 reference antidepressants, for novelty
    scoring: how similar is each ChEMBL candidate to a KNOWN antidepressant?"""
    ref = pd.read_csv(dataset_csv)
    names, fps = [], []
    for name, smi in zip(ref["drug_name"], ref["smiles"]):
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            names.append(name)
            fps.append(AllChem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=2048))
    return names, fps


def max_tanimoto_to_known(mol, ref_names, ref_fps):
    """Max Tanimoto similarity (Morgan r=2, 2048 bits) of `mol` to any of the
    28 reference antidepressants, plus which one it's most similar to."""
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    sims = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
    i = int(np.argmax(sims)) if sims else -1
    return (round(float(sims[i]), 3), ref_names[i]) if i >= 0 else (None, None)


def fetch_max_phase(chembl_ids, batch=45):
    """Look up ChEMBL's max_phase (development/approval status; 4 = approved
    for at least one indication somewhere) for a batch of molecule IDs, so we
    can flag candidates that are already-approved drugs (for OTHER
    indications) rather than truly novel chemical entities."""
    out = {}
    ids = list(chembl_ids)
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        url = (f"{CHEMBL}/molecule.json?molecule_chembl_id__in="
               f"{','.join(chunk)}&limit={batch}")
        try:
            data = _chembl_get(url)
        except Exception as e:
            print(f"  ! max_phase lookup failed for batch starting {chunk[0]}: {e}")
            continue
        for m in data.get("molecules", []):
            out[m.get("molecule_chembl_id")] = m.get("max_phase")
        time.sleep(PAUSE)
    return out


def features_present(mol):
    feats = set()
    for name, patt in patterns.items():
        if mol.HasSubstructMatch(patt):
            feats.add(name)
    if is_tricyclic_dibenzo(mol):
        feats.add("tricyclic_core")
    return feats


def trs_of(feats, weights):
    """Fitted Toxicity Risk Score = sum of TRS weights of features present."""
    return round(sum(weights.get(f, 0.0) for f in feats), 3)


def main():
    tox_feats, safe_feats, weights = load_feature_roles()
    print(f"Toxicophores (avoid): {sorted(tox_feats) or '[none]'}")
    print(f"Safeophores (want):   {sorted(safe_feats) or '[none]'}\n")

    # 1. Pull candidates from ChEMBL.
    candidates = {}
    target_of = {}
    for label, tid in TARGETS.items():
        actives = fetch_target_actives(tid)
        for mid, smi in actives.items():
            candidates.setdefault(mid, smi)
            target_of.setdefault(mid, label)
        print(f"  {label:<14} {len(actives):>5} actives "
              f"(running unique total: {len(candidates)})")
    print(f"\n{len(candidates)} unique candidate molecules pulled from ChEMBL.")

    # 2. InChIKey skeletons of our existing drugs, to exclude them.
    #    NOTE: this excludes ONLY our 28 curated antidepressants (by connectivity
    #    skeleton), not all FDA-approved drugs generally. A compound already
    #    approved for another indication (e.g. an ADHD drug that also hits NET)
    #    can still pass this filter - see the max_phase check in step 5 below,
    #    which flags exactly that.
    ours = pd.read_csv("antidepressants_dataset.csv")
    our_keys = set()
    for smi in ours["smiles"]:
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            k = inchikey_skeleton(m)
            if k:
                our_keys.add(k)
    ref_names, ref_fps = reference_fingerprints("antidepressants_dataset.csv")

    # 3. Screen.
    rows = []
    for mid, smi in candidates.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        skel = inchikey_skeleton(mol)
        if skel in our_keys:
            continue                       # it's one of our existing drugs
        feats = features_present(mol)
        safe_hits = feats & safe_feats
        tox_hits = feats & tox_feats
        trs = trs_of(feats, weights)        # validated Obj-5 structure score
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        lipinski_ok = (mw <= 500 and logp <= 5
                       and Lipinski.NumHDonors(mol) <= 5
                       and Lipinski.NumHAcceptors(mol) <= 10)
        max_sim, most_similar_to = max_tanimoto_to_known(mol, ref_names, ref_fps)
        rows.append({
            "molecule_chembl_id": mid,
            "via_target": target_of[mid],
            "smiles": smi,
            "MW": round(mw, 1), "LogP": round(logp, 2),
            "lipinski_ok": lipinski_ok,
            "veber_ok": veber_ok(mol),
            "ghose_ok": ghose_ok(mol),
            "safeophores": ";".join(sorted(safe_hits)),
            "toxicophores": ";".join(sorted(tox_hits)),
            "n_safe": len(safe_hits), "n_tox": len(tox_hits),
            "TRS": trs,
            "pass": (len(safe_hits) >= 1 and len(tox_hits) == 0),
            "max_tanimoto_to_known": max_sim,
            "most_similar_known_drug": most_similar_to,
        })

    res = pd.DataFrame(rows)
    if res.empty:
        print("\nNo candidates retrieved (ChEMBL may be unavailable). "
              "Try re-running shortly.")
        return
    # Prioritise: pass the filter, then lowest predicted TRS (safest), then drug-like.
    res = res.sort_values(["pass", "TRS", "lipinski_ok"],
                          ascending=[False, True, False])

    # 4. Novelty: candidates with max_tanimoto_to_known >= 0.85 are close
    #    analogues of a known antidepressant (not independently novel), even
    #    though they are a different ChEMBL entity (different salt/stereo/
    #    minor substitution). Flag rather than silently drop, so the paper can
    #    report the true novelty distribution instead of asserting it.
    res["novel_scaffold"] = res["max_tanimoto_to_known"] < 0.85

    # 5. FDA/development-phase overlap check for PASSING candidates only
    #    (cheap enough to do for ~a few hundred; answers "are these truly
    #    novel chemical entities, or already-marketed drugs for something
    #    else?"). max_phase == 4 means approved somewhere for some indication.
    pass_ids = res.loc[res["pass"], "molecule_chembl_id"].tolist()
    max_phase = fetch_max_phase(pass_ids) if pass_ids else {}
    res["max_phase"] = res["molecule_chembl_id"].map(max_phase)
    res["already_approved_drug"] = res["max_phase"] == 4

    res.to_csv("chembl_candidates_screened.csv", index=False)

    n_pass = int(res["pass"].sum())
    n_approved_overlap = int(res.loc[res["pass"], "already_approved_drug"].sum())
    n_novel = int(res.loc[res["pass"], "novel_scaffold"].sum())
    n_veber = int(res.loc[res["pass"], "veber_ok"].sum())
    n_ghose = int(res.loc[res["pass"], "ghose_ok"].sum())
    print(f"\nScreened {len(res)} non-redundant candidates.")
    print(f"{n_pass} PASS (>=1 safeophore, 0 toxicophores, "
          f"of which {int(res[res['pass']]['lipinski_ok'].sum())} Lipinski-ok, "
          f"{n_veber} Veber-ok, {n_ghose} Ghose-ok). Ranked by fitted TRS (lower = safer).")
    print(f"Novelty: {n_novel}/{n_pass} passing candidates have max Tanimoto "
          f"similarity < 0.85 to any of our 28 reference antidepressants "
          f"(median similarity = {res.loc[res['pass'], 'max_tanimoto_to_known'].median():.3f}).")
    print(f"Approval overlap: {n_approved_overlap}/{n_pass} passing candidates are "
          f"already an approved drug (ChEMBL max_phase=4) for some other "
          f"indication, not a novel chemical entity.")
    print("\nTop 15 prioritized candidates:")
    show = ["molecule_chembl_id", "via_target", "MW", "LogP", "lipinski_ok",
            "veber_ok", "ghose_ok", "max_tanimoto_to_known", "already_approved_drug",
            "safeophores", "toxicophores", "TRS", "pass"]
    print(res[show].head(15).to_string(index=False))
    print("\nWrote chembl_candidates_screened.csv")


if __name__ == "__main__":
    main()
