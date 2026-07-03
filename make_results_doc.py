"""Generate Results_Tables_and_Charts.docx: key tables + charts from outputs."""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from docx import Document
from docx.shared import Inches, Pt

from make_charts import generate_all, CH

ORGANS = ["cardiac", "hepatic", "neurologic", "metabolic"]

sig = pd.read_csv("faers_signals.csv")
rank = pd.read_csv("class_safety_ranking.csv")
trs = pd.read_csv("drug_trs.csv")
loo = pd.read_csv("trs_loocv.csv")
wts = pd.read_csv("trs_weights.csv")
cand = pd.read_csv("chembl_candidates_screened.csv")
corr = pd.read_csv("structure_toxicity_correlation_continuous.csv")

# Charts (single source of truth: make_charts.py) - also returns rho/p so the
# tables/captions below stay in sync with what the figures actually show.
_stats = generate_all()
RHO_IN, PVAL_IN = _stats["rho_in"], _stats["p_in"]
RHO_LOO, PVAL_LOO = _stats["rho_loo"], _stats["p_loo"]

# ================================================================ DOC
doc = Document()
doc.add_heading("Results: Tables and Charts", level=0)
doc.add_paragraph("Companion data appendix to the project summary. Generated "
                  "from the final pipeline outputs.")

def add_table(df, title, note=""):
    doc.add_heading(title, level=1)
    if note:
        doc.add_paragraph(note)
    t = doc.add_table(rows=1, cols=len(df.columns))
    t.style = "Light Grid Accent 1"
    for i, c in enumerate(df.columns):
        run = t.rows[0].cells[i].paragraphs[0].add_run(str(c))
        run.bold = True; run.font.size = Pt(8)
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for i, v in enumerate(row):
            p = cells[i].paragraphs[0]; r = p.add_run(str(v)); r.font.size = Pt(8)

# Table 1: flagged signals (PRIMARY threshold: ROR lower 95% bound > 2)
flagged = (sig[sig["signal"] == 1]
           [["drug_name", "organ_system", "n_cases", "ror", "ror_lower",
             "ror_upper", "IC", "IC025", "IC975"]]
           .round(2).sort_values("ror", ascending=False))
add_table(flagged, "Table 1. Primary disproportionality signals (Objective 1)",
          f"{len(flagged)} signals (ROR lower 95% bound > 2 and >= 3 cases; "
          f"legacy ROR>1 criterion flags {int(sig['signal_ror1'].sum())}, "
          "reported for comparison only, not used for conclusions), ranked "
          "by ROR, with 95% CI/credible intervals.")

# Table 2: within-class ranking
rk2 = rank[["class", "safety_rank_in_class", "drug_name"] + ORGANS +
           ["n_signals", "toxicity_index"]].sort_values(
    ["class", "safety_rank_in_class"])
add_table(rk2, "Table 2. Within-class safety ranking (Objective 2)",
          "Rank 1 = safest in class (lowest toxicity index). n_signals uses "
          "the primary (ROR>2) threshold. Index is RELATIVE to other "
          "antidepressants, not an absolute risk.")

# Table 3: structure-tox (significant), with BH-FDR q-values
cs = (corr[corr["p_value"] < 0.10]
      [["organ_system", "feature", "rank_biserial", "p_value", "p_adj_bh", "role"]]
      .sort_values("p_value"))
add_table(cs, "Table 3. Structure-toxicity associations, p < 0.10 (Objective 3)",
          f"Continuous-ROR Mann-Whitney test; q = Benjamini-Hochberg FDR "
          f"across all {len(corr)} organ x feature tests. "
          f"{int(corr['significant_fdr05'].sum())}/{len(corr)} survive q<0.05. "
          "Exploratory leads only.")

# Table 4: TRS weights
add_table(wts.round(3), "Table 4. Toxicity Risk Score weights (Objective 5)",
          "Per-feature weight = mean effect size on ROR across organ systems. "
          "Underlying tests do not survive FDR correction (Table 3).")

# Table 4b: TRS validation, in-sample vs LOOCV
trs_merged = trs.merge(loo[["drug_name", "TRS_loocv"]], on="drug_name",
                       how="left", suffixes=("", "_dup"))
add_table(trs_merged[["drug_name", "TRS_fitted", "TRS_loocv", "observed_tox_index"]],
          "Table 4b. TRS validation: in-sample fit vs. leave-one-drug-out prediction",
          f"In-sample (circular): Spearman rho={RHO_IN:+.3f}, p={PVAL_IN:.4f}. "
          f"Leave-one-drug-out (honest): rho={RHO_LOO:+.3f}, p={PVAL_LOO:.4f} - "
          "the score does not generalise out-of-sample.")

# Table 5: top candidates
top = (cand[cand["pass"]].head(15)
       [["molecule_chembl_id", "via_target", "MW", "LogP",
         "max_tanimoto_to_known", "already_approved_drug", "safeophores", "TRS"]])
add_table(top, "Table 5. Top prioritized candidates (Objective 4)",
          f"{int(cand['pass'].sum())} of {len(cand)} screened molecules pass; "
          f"{int(cand[cand['pass']]['novel_scaffold'].sum())} are genuinely "
          f"novel (Tanimoto<0.85 to any reference drug); "
          f"{int(cand[cand['pass']]['already_approved_drug'].sum())} already "
          "an approved drug elsewhere. Top 15 by lowest predicted TRS shown.")

# Charts
doc.add_heading("Charts", level=1)
for img, cap in [
    ("trs_validation.png",
     f"Figure 1. TRS in-sample fit (rho={RHO_IN:+.2f}, circular) vs. "
     f"leave-one-drug-out honest prediction (rho={RHO_LOO:+.2f}) - the score "
     "does not generalise out-of-sample."),
    ("ror_heatmap.png",
     "Figure 2. Reporting Odds Ratio for every drug and organ system. Solid "
     "outline = primary signal (ROR lower CI>2); dashed = legacy ROR>1 only."),
    ("trs_weights.png",
     "Figure 3. Fitted structural-feature weights. Red = toxicophore, blue = "
     "safeophore (colorblind-safe palette). None survive FDR correction."),
    ("class_toxicity.png",
     "Figure 4. Per-drug toxicity index coloured by pharmacological class."),
    ("feature_stability.png",
     "Figure 5. Structural-feature effect direction by organ system - sign "
     "flips show organ-specific, not global, toxicophore/safeophore roles."),
]:
    doc.add_picture(f"{CH}/{img}", width=Inches(5.6))
    cp = doc.add_paragraph(cap); cp.runs[0].font.size = Pt(9); cp.runs[0].italic = True

doc.save("Results_Tables_and_Charts.docx")
print("Wrote Results_Tables_and_Charts.docx with 6 tables + 5 charts")
print("Charts saved in:", os.path.abspath(CH))
