"""
Cheminformatics analysis for the antidepressant toxicity-structure study.
Covers Objectives 3-5: structural feature mapping, toxicophore/safeophore
identification, and the Toxicity Risk Score scaffold.

Run AFTER you have FAERS disproportionality signals (faers_signals.csv) so the
structure-toxicity correlation can be computed. The descriptor/substructure
parts run standalone with no FAERS data.

Usage:
    python cheminformatics_analysis.py
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Draw
from rdkit.Chem import Crippen, Lipinski


def benjamini_hochberg(pvals):
    """Benjamini-Hochberg FDR correction. Returns adjusted p-values (q-values),
    same order/length as input, with NaNs passed through untouched."""
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    n = len(pv)
    if n == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / (np.arange(1, n + 1))
    # enforce monotonicity (running minimum from the largest p-value down)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    adj = np.empty(n)
    adj[order] = q
    out[valid] = adj
    return out

# ----------------------------------------------------------------------
# 1. Load curated dataset
# ----------------------------------------------------------------------
df = pd.read_csv("antidepressants_dataset.csv")
df["mol"] = df["smiles"].apply(Chem.MolFromSmiles)
assert df["mol"].notna().all(), "Some SMILES failed to parse"


# ----------------------------------------------------------------------
# 2. Molecular descriptors (physicochemical context for tox)
# ----------------------------------------------------------------------
def compute_descriptors(mol):
    return pd.Series({
        "MW": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "RotBonds": Descriptors.NumRotatableBonds(mol),
        "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "FractionCSP3": rdMolDescriptors.CalcFractionCSP3(mol),
        "FormalCharge": Chem.GetFormalCharge(mol),
        "HeavyAtoms": mol.GetNumHeavyAtoms(),
    })


desc = df["mol"].apply(compute_descriptors)
df = pd.concat([df, desc], axis=1)


# ----------------------------------------------------------------------
# 3. Substructure / functional-group flags (candidate toxicophores)
#    SMARTS patterns for features mechanistically linked to organ toxicity.
#    These are HYPOTHESES to be tested against FAERS signals, not conclusions.
# ----------------------------------------------------------------------
SUBSTRUCTURE_SMARTS = {
    # --- cardiac-relevant (hERG / QT literature) ---
    "tertiary_amine":        "[NX3]([#6])([#6])[#6]",
    "basic_secondary_amine": "[NX3;H1]([#6])[#6]",
    # --- hepatic-relevant (reactive metabolite / structural alert) ---
    "aniline":               "[NX3;H2,H1]c1ccccc1",
    "para_halogen_phenyl":   "c1ccc(cc1)[F,Cl,Br]",
    "methylenedioxy":        "OCOc1ccccc1",               # MDP - bioactivation alert
    "thiophene":             "c1ccsc1",
    "furan":                 "c1ccoc1",
    # --- metabolic / mitochondrial-relevant ---
    "carboxylic_acid":       "C(=O)[OH]",
    "long_alkyl_chain":      "CCCCCC",                    # beta-oxidation substrate
    # --- generic alerts ---
    "trifluoromethyl":       "C(F)(F)F",
    "nitrile":               "C#N",
    "morpholine":            "N1CCOCC1",
    "piperazine":            "N1CCNCC1",
    "alpha_chloro_ketone":   "C(=O)C[Cl]",
    "propargylamine":        "C#CC[NX3]",                 # irreversible MAO alert
}

patterns = {k: Chem.MolFromSmarts(v) for k, v in SUBSTRUCTURE_SMARTS.items()}
for name, patt in patterns.items():
    if patt is None:
        raise ValueError(f"Bad SMARTS for {name}")
    df[f"has_{name}"] = df["mol"].apply(lambda m: int(m.HasSubstructMatch(patt)))


def is_tricyclic_dibenzo(m):
    """Fused tricyclic scaffold with >=2 aromatic 6-rings (TCA-like core).
    Ring-info based because fused-polycyclic SMARTS is brittle."""
    ri = m.GetRingInfo()
    arom6 = sum(1 for ring in ri.AtomRings()
                if len(ring) == 6
                and all(m.GetAtomWithIdx(a).GetIsAromatic() for a in ring))
    return int(arom6 >= 2 and rdMolDescriptors.CalcNumRings(m) >= 3)


df["has_tricyclic_core"] = df["mol"].apply(is_tricyclic_dibenzo)


# ----------------------------------------------------------------------
# 4. Pairwise structural similarity (Tanimoto on Morgan fingerprints)
#    Used to cluster drugs and sanity-check class groupings.
# ----------------------------------------------------------------------
from rdkit.Chem import AllChem
from rdkit import DataStructs

fps = [AllChem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=2048)
       for m in df["mol"]]
n = len(fps)
sim = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        sim[i, j] = DataStructs.TanimotoSimilarity(fps[i], fps[j])
sim_df = pd.DataFrame(sim, index=df["drug_name"], columns=df["drug_name"])


# ----------------------------------------------------------------------
# 5. Toxicophore / safeophore correlation (REQUIRES FAERS signals)
#    Merge structural flags with per-drug toxicity signal, then test which
#    substructures associate with high vs low signal per organ system.
# ----------------------------------------------------------------------
def correlate_with_signals(signal_csv="faers_signals.csv"):
    """
    signal_csv must have columns:
        drug_name, organ_system, ror, ror_lower, ror_upper, n_cases, signal(0/1)
    organ_system in {cardiac, hepatic, neurologic, metabolic}
    """
    sig = pd.read_csv(signal_csv)
    flag_cols = [c for c in df.columns if c.startswith("has_")]
    merged = sig.merge(df[["drug_name"] + flag_cols], on="drug_name")

    results = []
    for organ in merged["organ_system"].unique():
        sub = merged[merged["organ_system"] == organ]
        for fc in flag_cols:
            present = sub[sub[fc] == 1]["signal"]
            absent  = sub[sub[fc] == 0]["signal"]
            if len(present) == 0 or len(absent) == 0:
                continue
            # 2x2: feature present/absent vs signal yes/no -> Fisher exact
            from scipy.stats import fisher_exact
            a = int(((sub[fc] == 1) & (sub["signal"] == 1)).sum())
            b = int(((sub[fc] == 1) & (sub["signal"] == 0)).sum())
            c = int(((sub[fc] == 0) & (sub["signal"] == 1)).sum())
            d = int(((sub[fc] == 0) & (sub["signal"] == 0)).sum())
            try:
                odds, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
            except Exception:
                odds, p = np.nan, np.nan
            results.append({
                "organ_system": organ, "feature": fc.replace("has_", ""),
                "n_with_feature": a + b, "signals_with_feature": a,
                "odds_ratio": odds, "p_value": p,
                "role": ("toxicophore" if odds > 1 else "safeophore")
                         if not np.isnan(odds) else "indeterminate",
            })
    res = pd.DataFrame(results).sort_values(["organ_system", "p_value"])
    # Benjamini-Hochberg FDR across ALL organ x feature tests (one family).
    res["p_adj_bh"] = benjamini_hochberg(res["p_value"].values)
    res["significant_fdr05"] = res["p_adj_bh"] < 0.05
    res.to_csv("structure_toxicity_correlation.csv", index=False)
    return res


def correlate_continuous(signal_csv="faers_signals.csv", sig_df=None,
                         exclude_drug=None, save_csv=True):
    """
    Continuous-ROR version of the structure-toxicity test. Instead of collapsing
    each drug to a 0/1 signal (which discards information and produces degenerate
    2x2 tables on only ~27 drugs), compare the full distribution of log(ROR)
    between drugs that CARRY a structural feature and those that do not, per
    organ system, with a Mann-Whitney U test. Uses every drug, not just signals.

    Reports:
        median_ror_with / median_ror_without : effect direction/size on ROR scale
        rank_biserial : non-parametric effect size in [-1, 1]
        p_value       : Mann-Whitney U (two-sided)
    A feature whose carriers have systematically higher ROR (rank_biserial > 0)
    is a candidate toxicophore for that organ; lower => candidate safeophore.

    sig_df / exclude_drug let callers (LOOCV, sensitivity analyses) re-run this
    test on an in-memory subset of drugs without touching faers_signals.csv.
    """
    from scipy.stats import mannwhitneyu

    sig = sig_df.copy() if sig_df is not None else pd.read_csv(signal_csv)
    sig = sig[sig["ror"].notna() & np.isfinite(sig["ror"])]
    if exclude_drug is not None:
        excl = {exclude_drug} if isinstance(exclude_drug, str) else set(exclude_drug)
        sig = sig[~sig["drug_name"].isin(excl)]
    flag_cols = [c for c in df.columns if c.startswith("has_")]
    merged = sig.merge(df[["drug_name"] + flag_cols], on="drug_name")

    results = []
    for organ in merged["organ_system"].unique():
        sub = merged[merged["organ_system"] == organ]
        for fc in flag_cols:
            with_f = sub[sub[fc] == 1]["ror"].astype(float)
            without_f = sub[sub[fc] == 0]["ror"].astype(float)
            # need carriers on both sides to compare
            if len(with_f) < 2 or len(without_f) < 2:
                continue
            try:
                u, p = mannwhitneyu(with_f, without_f, alternative="two-sided")
                # rank-biserial effect size from U
                rb = 1 - (2 * u) / (len(with_f) * len(without_f))
                rb = -rb  # orient so >0 means carriers have HIGHER ror
            except Exception:
                u, p, rb = np.nan, np.nan, np.nan
            results.append({
                "organ_system": organ, "feature": fc.replace("has_", ""),
                "n_with_feature": int(len(with_f)),
                "median_ror_with": round(float(with_f.median()), 3),
                "median_ror_without": round(float(without_f.median()), 3),
                "rank_biserial": round(float(rb), 3),
                "p_value": round(float(p), 4),
                "role": ("toxicophore" if rb > 0 else "safeophore")
                        if not np.isnan(rb) else "indeterminate",
            })
    res = pd.DataFrame(results).sort_values(["organ_system", "p_value"])
    # Benjamini-Hochberg FDR across ALL organ x feature tests (one family of
    # 4 organs x 15 features = up to 60 comparisons), not per-organ.
    res["p_adj_bh"] = benjamini_hochberg(res["p_value"].values)
    res["significant_fdr05"] = res["p_adj_bh"] < 0.05
    if save_csv:
        res.to_csv("structure_toxicity_correlation_continuous.csv", index=False)
    return res


# ----------------------------------------------------------------------
# 5b. Multicollinearity among structural features (e.g. tricyclic core vs
#     tertiary amine both mark the TCA scaffold). Reported as phi coefficient
#     (Pearson correlation on 0/1 flags) across all 28 drugs, flagging any
#     |phi| >= 0.5 pair as materially redundant for TRS-weight interpretation.
# ----------------------------------------------------------------------
def multicollinearity_check(flag_threshold=0.5):
    flag_cols = [c for c in df.columns if c.startswith("has_")]
    mat = df[flag_cols].corr()  # Pearson on 0/1 columns == phi coefficient
    mat.to_csv("structural_feature_multicollinearity.csv")

    pairs = []
    for i, f1 in enumerate(flag_cols):
        for f2 in flag_cols[i + 1:]:
            phi = mat.loc[f1, f2]
            if pd.notna(phi) and abs(phi) >= flag_threshold:
                pairs.append({
                    "feature_1": f1.replace("has_", ""),
                    "feature_2": f2.replace("has_", ""),
                    "phi": round(float(phi), 3),
                })
    flagged = pd.DataFrame(pairs).sort_values(
        "phi", key=lambda s: s.abs(), ascending=False) if pairs else pd.DataFrame(
        columns=["feature_1", "feature_2", "phi"])
    flagged.to_csv("structural_feature_multicollinearity_flagged.csv", index=False)
    return mat, flagged


# ----------------------------------------------------------------------
# 6. Toxicity Risk Score scaffold (Objective 5)
#    Simple, transparent, additive score. Weights are PLACEHOLDERS until
#    fitted/validated against the FAERS signal correlation above.
# ----------------------------------------------------------------------
def toxicity_risk_score(mol, weights=None):
    """
    Transparent additive structure-based score. Replace placeholder weights
    with values derived from the FAERS correlation (e.g. log odds ratios).
    """
    if weights is None:
        # placeholder weights: +ve = toxicophore, -ve = safeophore
        weights = {
            "tricyclic_core": 1.5,
            "tertiary_amine": 0.8, "para_halogen_phenyl": 0.5,
            "methylenedioxy": 1.0, "alpha_chloro_ketone": 1.2,
            "long_alkyl_chain": 0.7, "propargylamine": 1.3,
            "morpholine": -0.8, "carboxylic_acid": -0.5,
            "nitrile": -0.3, "FractionCSP3": -1.0,  # continuous, scaled below
        }
    score = 0.0
    for feat, w in weights.items():
        if feat == "tricyclic_core":
            score += w * is_tricyclic_dibenzo(mol)
        elif feat in SUBSTRUCTURE_SMARTS:
            patt = patterns[feat]
            score += w * int(mol.HasSubstructMatch(patt))
        elif feat == "FractionCSP3":
            score += w * rdMolDescriptors.CalcFractionCSP3(mol)
    return score


df["TRS_placeholder"] = df["mol"].apply(toxicity_risk_score)


# ----------------------------------------------------------------------
# 7. Outputs
# ----------------------------------------------------------------------
if __name__ == "__main__":
    out_cols = (["drug_name", "class", "MW", "LogP", "TPSA", "FractionCSP3",
                 "TRS_placeholder"]
                + [c for c in df.columns if c.startswith("has_")])
    df[out_cols].to_csv("drug_features.csv", index=False)
    sim_df.to_csv("structural_similarity_matrix.csv")
    print("Wrote drug_features.csv and structural_similarity_matrix.csv")
    print(f"\n{len(df)} molecules, {len(SUBSTRUCTURE_SMARTS)} substructure flags")
    print("\nFeature prevalence across the set:")
    print(df[[c for c in df.columns if c.startswith("has_")]].sum()
          .sort_values(ascending=False).to_string())

    mat, flagged = multicollinearity_check()
    print(f"\nStructural-feature multicollinearity: {len(flagged)} pair(s) with "
          f"|phi| >= 0.5 (wrote structural_feature_multicollinearity.csv, "
          f"..._flagged.csv)")
    if len(flagged):
        print(flagged.to_string(index=False))

    import os
    if os.path.exists("faers_signals.csv"):
        print("\nFAERS signals found - structure-toxicity correlation (binary Fisher)")
        print(correlate_with_signals().to_string())
        print("\nStructure-toxicity correlation (continuous ROR, Mann-Whitney)")
        cc = correlate_continuous()
        print(cc.to_string())
        print(f"\nBH-FDR: {int(cc['significant_fdr05'].sum())} of {len(cc)} tests "
              f"survive q<0.05 (out of 4 organs x {len(SUBSTRUCTURE_SMARTS)+1} features)")
    else:
        print("\n[FAERS signals not yet present - run faers_disproportionality.py "
              "first, then re-run for Objectives 3-5 correlation.]")
