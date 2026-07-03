"""
Objective 5: fit a transparent, structure-based Toxicity Risk Score (TRS) from
the FAERS-derived structure-toxicity correlation, replacing the placeholder
weights in cheminformatics_analysis.py.

Method (deliberately simple and interpretable):
  - For each structural feature, its weight = MEAN rank-biserial effect size
    across the four organ systems (from structure_toxicity_correlation_continuous
    .csv). rank-biserial > 0 means carriers tend to have HIGHER ROR (toxicophore,
    positive weight); < 0 means lower ROR (safeophore, negative weight).
  - TRS(drug) = sum of weights of the structural features it carries.
    Higher TRS => predicted higher overall toxicity burden.

Validation:
  - Spearman correlation between fitted TRS and each drug's OBSERVED toxicity
    index (mean ln(ROR) over the 4 organ systems). A positive, significant
    correlation means the structure-only score tracks real-world signals.

Writes trs_weights.csv (the proposed scoring scheme) and drug_trs.csv
(per-drug fitted score vs observed toxicity).

Usage:
    python fit_trs.py
"""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# df already carries all has_<feature> flags + descriptors.
from cheminformatics_analysis import df, correlate_continuous

ORGANS = ["cardiac", "hepatic", "neurologic", "metabolic"]


def fit_weights(corr_csv="structure_toxicity_correlation_continuous.csv"):
    corr = pd.read_csv(corr_csv)
    # weight per feature = mean signed effect size across organ systems
    w = (corr.groupby("feature")["rank_biserial"].mean()
             .rename("weight").reset_index())
    w = w.sort_values("weight", ascending=False)
    w.to_csv("trs_weights.csv", index=False)
    return dict(zip(w["feature"], w["weight"]))


def trs_from_weights(weights):
    """Fitted TRS per drug (Series indexed like df) for an arbitrary weight dict."""
    flag_cols = [c for c in df.columns if c.startswith("has_")]
    trs = pd.Series(0.0, index=df.index)
    for c in flag_cols:
        feat = c.replace("has_", "")
        if feat in weights:
            trs += df[c] * weights[feat]
    return trs


def loocv(signal_csv="faers_signals.csv"):
    """
    Leave-one-drug-out cross-validation for the TRS.

    Fitting and validating the TRS on the same 27 drugs (as fit_weights() +
    the Spearman check in main() do) is circular: the feature weights are
    themselves derived from those drugs' ROR values, so the "validation"
    partly re-tests the data used to build the score.

    For each drug d in turn:
      1. Recompute the organ x feature rank-biserial effect sizes (and hence
         the feature weights) using every OTHER drug's ROR values only.
      2. Score d with those held-out weights (d's own ROR values never touch
         the weights that score it).
      3. Record (predicted TRS, observed toxicity index) for d.

    The out-of-sample Spearman correlation across all held-out predictions is
    the honest estimate of how well structure predicts toxicity; compare it
    to the in-sample (circular) rho printed by main().
    """
    sig_all = pd.read_csv(signal_csv)
    obs = observed_toxicity_index(signal_csv).reset_index()
    drugs = obs["drug_name"].tolist()

    rows = []
    for held_out in drugs:
        corr_loo = correlate_continuous(sig_df=sig_all, exclude_drug=held_out,
                                        save_csv=False)
        w_loo = (corr_loo.groupby("feature")["rank_biserial"].mean()).to_dict()
        trs_loo = trs_from_weights(w_loo)
        drug_mask = df["drug_name"] == held_out
        if not drug_mask.any():
            continue
        pred = float(trs_loo[drug_mask].iloc[0])
        actual = float(obs.loc[obs["drug_name"] == held_out, "observed_tox_index"].iloc[0])
        rows.append({"drug_name": held_out, "TRS_loocv": round(pred, 3),
                     "observed_tox_index": round(actual, 3)})

    out = pd.DataFrame(rows)
    rho, p = spearmanr(out["TRS_loocv"], out["observed_tox_index"])
    out.to_csv("trs_loocv.csv", index=False)
    return out, rho, p


def observed_toxicity_index(signal_csv="faers_signals.csv"):
    sig = pd.read_csv(signal_csv)
    wide = sig.pivot_table(index="drug_name", columns="organ_system",
                           values="ror", aggfunc="first").reindex(columns=ORGANS)
    return np.log(wide.clip(lower=1e-6)).mean(axis=1).rename("observed_tox_index")


def main():
    weights = fit_weights()
    print("Fitted TRS weights (per structural feature, + = toxicophore):")
    for f, wt in sorted(weights.items(), key=lambda kv: -kv[1]):
        print(f"  {f:<22} {wt:+.3f}")

    # Compute fitted TRS for each drug from its structural flags.
    flag_cols = [c for c in df.columns if c.startswith("has_")]
    trs = pd.Series(0.0, index=df.index)
    for c in flag_cols:
        feat = c.replace("has_", "")
        if feat in weights:
            trs += df[c] * weights[feat]
    scored = df[["drug_name"]].copy()
    scored["TRS_fitted"] = trs.round(3)

    # Join observed toxicity and validate.
    obs = observed_toxicity_index().reset_index()
    merged = scored.merge(obs, on="drug_name", how="inner")
    rho, p = spearmanr(merged["TRS_fitted"], merged["observed_tox_index"])

    # Leave-one-drug-out cross-validation (honest, non-circular estimate).
    loo, rho_loo, p_loo = loocv()
    merged = merged.merge(loo[["drug_name", "TRS_loocv"]], on="drug_name", how="left")

    merged = merged.sort_values("TRS_fitted", ascending=False)
    merged["TRS_fitted"] = merged["TRS_fitted"].round(3)
    merged["observed_tox_index"] = merged["observed_tox_index"].round(3)
    merged.to_csv("drug_trs.csv", index=False)

    print(f"\nIn-sample validation (CIRCULAR - weights were fit on these same "
          f"drugs): Spearman(TRS_fitted, observed_tox_index) "
          f"rho = {rho:+.3f}, p = {p:.4f}  (n = {len(merged)})")
    print(f"Leave-one-drug-out validation (honest, out-of-sample): "
          f"Spearman(TRS_loocv, observed_tox_index) "
          f"rho = {rho_loo:+.3f}, p = {p_loo:.4f}  (n = {len(loo)})")
    print("  (positive rho => structure-only score tracks real FAERS toxicity)\n")
    print("Per-drug fitted TRS vs observed toxicity (top/bottom):")
    print(merged.head(8).to_string(index=False))
    print("  ...")
    print(merged.tail(5).to_string(index=False))
    print("\nWrote trs_weights.csv, drug_trs.csv (incl. TRS_loocv column) and trs_loocv.csv")


if __name__ == "__main__":
    main()
