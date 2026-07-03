"""
Independent sanity check of the full pipeline. Re-derives key numbers from the
raw data (NOT from the pipeline's own outputs) and flags inconsistencies.
"""

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi
from faers_disproportionality import SYN_TO_DRUG, PT_TO_ORGAN, ORGAN_SYSTEM_PT

PASS, FAIL = "PASS", "**FAIL**"
issues = []


def check(name, ok, detail=""):
    print(f"[{PASS if ok else FAIL}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        issues.append(name)


# 1. Raw data integrity ------------------------------------------------------
rep = pd.read_csv("faers_reports.csv", dtype=str)
check("faers_reports.csv has 3 expected columns",
      list(rep.columns) == ["primaryid", "drugname", "pt"], str(list(rep.columns)))
check("no null primaryid/drugname/pt",
      rep[["primaryid", "drugname", "pt"]].notna().all().all())
check("rows de-duplicated", not rep.duplicated().any(),
      f"{rep.duplicated().sum()} dup rows")
n_reports = rep["primaryid"].nunique()
print(f"      {len(rep)} rows, {n_reports} unique reports, "
      f"{rep['drugname'].nunique()} drugs")

# 2. Independently recompute ROR for one drug x organ, compare to pipeline ----
sig = pd.read_csv("faers_signals.csv")
rep["drug"] = rep["drugname"].map(SYN_TO_DRUG)
rep["organ"] = rep["pt"].map(PT_TO_ORGAN)
r = rep.dropna(subset=["drug", "organ"])
total = r["primaryid"].nunique()

# pick the highest-n_cases signal row to test
row = sig.sort_values("n_cases", ascending=False).iloc[0]
drug, organ = row["drug_name"], row["organ_system"]
a = r[(r.drug == drug) & (r.organ == organ)]["primaryid"].nunique()
drug_tot = r[r.drug == drug]["primaryid"].nunique()
organ_tot = r[r.organ == organ]["primaryid"].nunique()
b, c = drug_tot - a, organ_tot - a
d = total - a - b - c
a2, b2, c2, d2 = [x if x else 0.5 for x in (a, b, c, d)]
ror_indep = (a2 / b2) / (c2 / d2)
check(f"ROR re-derived for {drug}/{organ}",
      np.isclose(ror_indep, row["ror"], rtol=1e-3),
      f"independent={ror_indep:.4f} vs pipeline={row['ror']:.4f}")
check(f"n_cases matches for {drug}/{organ}", a == int(row["n_cases"]),
      f"independent a={a} vs pipeline n_cases={row['n_cases']}")

# 3. Signal flag logic -------------------------------------------------------
recomputed = ((sig["ror_lower"] > 2) & (sig["n_cases"] >= 3)).astype(int)
check("primary signal flag = (ror_lower>2 AND n_cases>=3)",
      (recomputed == sig["signal"]).all())
recomputed_legacy = ((sig["ror_lower"] > 1) & (sig["n_cases"] >= 3)).astype(int)
check("legacy signal_ror1 flag = (ror_lower>1 AND n_cases>=3)",
      (recomputed_legacy == sig["signal_ror1"]).all())
check("all ROR positive & finite",
      np.isfinite(sig["ror"]).all() and (sig["ror"] > 0).all())
check("IC975 >= IC025 (valid credible interval) for all rows",
      (sig["IC975"] >= sig["IC025"]).all())

# 4. Corrected SMILES actually fixed ----------------------------------------
ds = pd.read_csv("antidepressants_dataset.csv")
cit = ds[ds.drug_name == "Citalopram"]["smiles"].iloc[0]
m = Chem.MolFromSmiles(cit)
check("citalopram SMILES parses", m is not None)
check("citalopram skeleton matches authoritative (WSEQXVZVJXJVFP)",
      m is not None and inchi.MolToInchiKey(m)[:14] == "WSEQXVZVJXJVFP",
      inchi.MolToInchiKey(m)[:14] if m else "n/a")
# para-halogen-phenyl SMARTS should now match citalopram (4-fluorophenyl)
phf = Chem.MolFromSmarts("c1ccc(cc1)[F,Cl,Br]")
check("citalopram now matches para_halogen_phenyl (fixed structure)",
      m is not None and m.HasSubstructMatch(phf))

# 5. TRS validation reproduces -----------------------------------------------
from scipy.stats import spearmanr
trs = pd.read_csv("drug_trs.csv")
rho, p = spearmanr(trs["TRS_fitted"], trs["observed_tox_index"])
check("TRS in-sample (circular) correlation reproduces (rho>0, p<0.05)",
      rho > 0 and p < 0.05, f"rho={rho:+.3f}, p={p:.4f}")
check("drug_trs.csv carries the leave-one-out column",
      "TRS_loocv" in trs.columns)
loo = pd.read_csv("trs_loocv.csv")
rho_loo, p_loo = spearmanr(loo["TRS_loocv"], loo["observed_tox_index"])
print(f"      [INFO] LOOCV (honest, out-of-sample): rho={rho_loo:+.3f}, "
      f"p={p_loo:.4f}  <- NOT expected to be significant; this is what makes "
      f"the in-sample number circular")

# 5b. Multicollinearity + FDR outputs exist ----------------------------------
import os
check("structural_feature_multicollinearity.csv written",
      os.path.exists("structural_feature_multicollinearity.csv"))
corr_cont = pd.read_csv("structure_toxicity_correlation_continuous.csv")
check("continuous structure-toxicity results carry BH-FDR columns",
      {"p_adj_bh", "significant_fdr05"}.issubset(corr_cont.columns))

# 5c. PubChem structure verification covers all 28 drugs ---------------------
if os.path.exists("pubchem_verification.csv"):
    pcv = pd.read_csv("pubchem_verification.csv")
    check("PubChem verification covers all 28 drugs",
          len(pcv) == len(ds))
    check("all 28 drug structures match PubChem skeleton",
          pcv["skeleton_match"].all(),
          f"{int(pcv['skeleton_match'].sum())}/{len(pcv)} match")

# 6. ChEMBL screen excludes our own drugs -----------------------------------
cand = pd.read_csv("chembl_candidates_screened.csv")
our_keys = set()
for s in ds["smiles"]:
    mm = Chem.MolFromSmiles(s)
    if mm:
        our_keys.add(inchi.MolToInchiKey(mm)[:14])
overlap = []
for s in cand["smiles"]:
    mm = Chem.MolFromSmiles(str(s))
    if mm and inchi.MolToInchiKey(mm)[:14] in our_keys:
        overlap.append(s)
check("no existing drug leaked into candidate list (incl. fixed citalopram)",
      len(overlap) == 0, f"{len(overlap)} overlaps")
check("passing candidates carry >=1 safeophore and 0 toxicophores",
      ((cand[cand["pass"]]["n_safe"] >= 1) &
       (cand[cand["pass"]]["n_tox"] == 0)).all())
check("candidate screen carries novelty (Tanimoto) + Veber/Ghose + approval columns",
      {"max_tanimoto_to_known", "novel_scaffold", "veber_ok", "ghose_ok",
       "already_approved_drug"}.issubset(cand.columns))

# 7. Within-class ranking covers all drugs with data ------------------------
rank = pd.read_csv("class_safety_ranking.csv")
check("ranking has every drug that has FAERS signals",
      set(rank["drug_name"]) == set(sig["drug_name"]),
      f"ranked {rank['drug_name'].nunique()} vs signalled {sig['drug_name'].nunique()}")

print("\n" + ("ALL CHECKS PASSED" if not issues
             else f"ISSUES: {issues}"))
