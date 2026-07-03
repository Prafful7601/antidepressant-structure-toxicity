"""
Objective 2: rank antidepressants by safety profile WITHIN each pharmacological
class (SSRI, TCA, Atypical, RIMA, SNRI).

Reads:
    faers_signals.csv          (per drug x organ: ror, signal, ...)
    antidepressants_dataset.csv (drug_name -> class, mechanism)

For each drug we summarise its four organ-system Reporting Odds Ratios into a
single, transparent toxicity index:

    toxicity_index = mean over {cardiac, hepatic, neurologic, metabolic} of ln(ROR)

    ln(ROR) is used so the scale is symmetric around 0:
        > 0  => drug is OVER-reported for toxicity vs other antidepressants (less safe)
        = 0  => reported at the class-average rate
        < 0  => UNDER-reported (relatively safer)

We also report n_signals = how many of the 4 organ systems were flagged as a
formal signal (ROR lower bound > 1 AND >= 3 cases).

Within each class, rank 1 = SAFEST (lowest toxicity index).

Writes class_safety_ranking.csv and prints a per-class table.

Usage:
    python rank_within_class.py
"""

import numpy as np
import pandas as pd

ORGANS = ["cardiac", "hepatic", "neurologic", "metabolic"]


def main():
    sig = pd.read_csv("faers_signals.csv")
    meta = pd.read_csv("antidepressants_dataset.csv")[
        ["drug_name", "class", "mechanism"]]

    # Wide table: one row per drug, one ROR column per organ system.
    wide = sig.pivot_table(index="drug_name", columns="organ_system",
                           values="ror", aggfunc="first")
    wide = wide.reindex(columns=ORGANS)

    # Signal count per drug (how many organ systems flagged).
    n_sig = (sig.assign(signal=pd.to_numeric(sig["signal"], errors="coerce"))
                .groupby("drug_name")["signal"].sum().rename("n_signals"))

    # Composite toxicity index = mean ln(ROR) across organ systems.
    ln = np.log(wide.clip(lower=1e-6))
    tox = ln.mean(axis=1).rename("toxicity_index")

    out = (wide.join(tox).join(n_sig)
              .reset_index()
              .merge(meta, on="drug_name", how="left"))
    out["n_signals"] = out["n_signals"].fillna(0).astype(int)

    # Round ROR columns + index for readability.
    for c in ORGANS:
        out[c] = out[c].round(2)
    out["toxicity_index"] = out["toxicity_index"].round(3)

    # Rank within class: 1 = safest (lowest toxicity index).
    out["safety_rank_in_class"] = (out.groupby("class")["toxicity_index"]
                                      .rank(method="min").astype(int))
    out = out.sort_values(["class", "toxicity_index"])

    cols = (["class", "safety_rank_in_class", "drug_name", "mechanism"]
            + ORGANS + ["n_signals", "toxicity_index"])
    out = out[cols]
    out.to_csv("class_safety_ranking.csv", index=False)

    # Pretty per-class print.
    print("Within-class safety ranking  (rank 1 = safest; toxicity_index = "
          "mean ln(ROR) over 4 organ systems)\n")
    for cls, grp in out.groupby("class"):
        print(f"=== {cls} ===")
        print(grp.drop(columns="class").to_string(index=False))
        safest = grp.iloc[0]
        riskiest = grp.iloc[-1]
        print(f"  -> safest: {safest['drug_name']} "
              f"(index {safest['toxicity_index']});  "
              f"highest-risk: {riskiest['drug_name']} "
              f"(index {riskiest['toxicity_index']})\n")

    missing = sorted(set(meta["drug_name"]) - set(out["drug_name"]))
    if missing:
        print(f"Not ranked (no FAERS signal data): {', '.join(missing)}")
    print("\nWrote class_safety_ranking.csv")


if __name__ == "__main__":
    main()
