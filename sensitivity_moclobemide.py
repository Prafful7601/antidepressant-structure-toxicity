"""
Sensitivity analysis: does moclobemide's very high neurologic ROR (3.93) - a
clear outlier vs the rest of the class - drive the headline results (the
structure-toxicity associations and the TRS validation), or are they robust
to its removal?

Re-runs the structure-toxicity correlation and TRS fit/validation (both
in-sample and leave-one-out) with moclobemide excluded from the drug set,
and reports the fitted/LOOCV rho and the top structure-toxicity associations
side-by-side with the full-set numbers.

Usage:
    python sensitivity_moclobemide.py
"""

import pandas as pd
from scipy.stats import spearmanr

from cheminformatics_analysis import df, correlate_continuous
from fit_trs import trs_from_weights, observed_toxicity_index, loocv

EXCLUDE = "Moclobemide"


def main():
    sig_all = pd.read_csv("faers_signals.csv")
    sig_excl = sig_all[sig_all["drug_name"] != EXCLUDE]

    # ---- Full-set (already computed elsewhere, reload for the comparison) ----
    corr_full = pd.read_csv("structure_toxicity_correlation_continuous.csv")
    trs_full = pd.read_csv("drug_trs.csv")
    rho_full_insample, p_full_insample = spearmanr(
        trs_full["TRS_fitted"], trs_full["observed_tox_index"])

    # ---- Excluding moclobemide ----
    corr_excl = correlate_continuous(sig_df=sig_excl, save_csv=False)
    w_excl = corr_excl.groupby("feature")["rank_biserial"].mean().to_dict()
    trs_excl = trs_from_weights(w_excl)
    obs = pd.read_csv("faers_signals.csv")
    obs_excl = obs[obs["drug_name"] != EXCLUDE]
    wide = obs_excl.pivot_table(index="drug_name", columns="organ_system",
                                values="ror", aggfunc="first")
    import numpy as np
    obs_idx = np.log(wide.clip(lower=1e-6)).mean(axis=1).rename("observed_tox_index").reset_index()

    scored = df[["drug_name"]].copy()
    scored["TRS_fitted_excl_moclobemide"] = trs_excl.round(3)
    merged = scored.merge(obs_idx, on="drug_name", how="inner")
    merged = merged[merged["drug_name"] != EXCLUDE]
    rho_excl, p_excl = spearmanr(merged["TRS_fitted_excl_moclobemide"],
                                 merged["observed_tox_index"])

    # LOOCV with moclobemide excluded from the drug universe entirely.
    loo_excl_rows = []
    drugs = merged["drug_name"].tolist()
    for held_out in drugs:
        sig_loo = sig_excl[sig_excl["drug_name"] != held_out]
        corr_loo = correlate_continuous(sig_df=sig_loo, save_csv=False)
        w_loo = corr_loo.groupby("feature")["rank_biserial"].mean().to_dict()
        pred = float(trs_from_weights(w_loo)[df["drug_name"] == held_out].iloc[0])
        actual = float(merged.loc[merged["drug_name"] == held_out, "observed_tox_index"].iloc[0])
        loo_excl_rows.append({"drug_name": held_out, "TRS_loocv": pred,
                              "observed_tox_index": actual})
    loo_excl = pd.DataFrame(loo_excl_rows)
    rho_loo_excl, p_loo_excl = spearmanr(loo_excl["TRS_loocv"], loo_excl["observed_tox_index"])

    print("=" * 78)
    print(f"SENSITIVITY ANALYSIS: excluding {EXCLUDE} (neurologic ROR 3.93, the "
          f"strongest single signal in the dataset)")
    print("=" * 78)
    print(f"\nTRS validation, in-sample (circular) Spearman rho:")
    print(f"  full set (n={len(trs_full)}):            rho = {rho_full_insample:+.3f}, p = {p_full_insample:.4f}")
    print(f"  excl. {EXCLUDE} (n={len(merged)}):  rho = {rho_excl:+.3f}, p = {p_excl:.4f}")
    print(f"\nTRS validation, leave-one-out (honest) Spearman rho:")
    _, rho_loo_full, p_loo_full = loocv()
    print(f"  full set (n={len(trs_full)}):            rho = {rho_loo_full:+.3f}, p = {p_loo_full:.4f}")
    print(f"  excl. {EXCLUDE} (n={len(loo_excl)}):  rho = {rho_loo_excl:+.3f}, p = {p_loo_excl:.4f}")

    print(f"\nTop structure-toxicity associations, full set vs excl. {EXCLUDE} "
          f"(by |rank_biserial|, neurologic organ only - the organ moclobemide "
          f"dominates):")
    for label, c in [("full set", corr_full), (f"excl. {EXCLUDE}", corr_excl)]:
        top = (c[c["organ_system"] == "neurologic"]
               .reindex(c.columns, axis=1)
               .sort_values("p_value").head(5))
        print(f"\n  -- {label} --")
        print(top[["feature", "rank_biserial", "p_value", "role"]].to_string(index=False))

    out = merged.merge(loo_excl[["drug_name", "TRS_loocv"]], on="drug_name", how="left")
    out.to_csv("sensitivity_excl_moclobemide.csv", index=False)
    corr_excl.to_csv("structure_toxicity_correlation_continuous_excl_moclobemide.csv",
                     index=False)
    print(f"\nWrote sensitivity_excl_moclobemide.csv and "
          f"structure_toxicity_correlation_continuous_excl_moclobemide.csv")


if __name__ == "__main__":
    main()
