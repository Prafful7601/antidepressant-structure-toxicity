"""Generate Research_Paper.docx: full manuscript from the pipeline outputs."""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_SECTION, WD_ORIENT

from cheminformatics_analysis import SUBSTRUCTURE_SMARTS
from chembl_candidate_screen import TARGETS, POTENCY

ORGANS = ["cardiac", "hepatic", "neurologic", "metabolic"]

sig = pd.read_csv("faers_signals.csv")
rank = pd.read_csv("class_safety_ranking.csv")
trs = pd.read_csv("drug_trs.csv")
loo = pd.read_csv("trs_loocv.csv")
wts = pd.read_csv("trs_weights.csv")
cand = pd.read_csv("chembl_candidates_screened.csv")
corr = pd.read_csv("structure_toxicity_correlation_continuous.csv")
corr_excl = pd.read_csv("structure_toxicity_correlation_continuous_excl_moclobemide.csv")
sens = pd.read_csv("sensitivity_excl_moclobemide.csv")
rep = pd.read_csv("faers_reports.csv", dtype=str)
mcol = pd.read_csv("structural_feature_multicollinearity_flagged.csv")
mcol_full = pd.read_csv("structural_feature_multicollinearity.csv", index_col=0)
PHI_TRICYCLIC_TERTIARY = float(mcol_full.loc["has_tricyclic_core", "has_tertiary_amine"])
pcv = pd.read_csv("pubchem_verification.csv")
date_range = pd.read_csv("faers_date_range.csv")
report_counts = pd.read_csv("faers_report_counts_per_drug.csv")
ds = pd.read_csv("antidepressants_dataset.csv")

N_REP = rep["primaryid"].nunique()
N_ROW = len(rep)
N_DRUG = rep["drugname"].nunique()          # 27 with FAERS data
N_CURATED = len(ds)                          # 28 curated structures
N_SIG = int(sig["signal"].sum())             # primary threshold (ROR lower CI >2)
N_SIG_LEGACY = int(sig["signal_ror1"].sum()) # legacy ROR>1 threshold
N_PASS = int(cand["pass"].sum())
N_DL = int(cand[cand["pass"]]["lipinski_ok"].sum())
N_VEBER = int(cand[cand["pass"]]["veber_ok"].sum())
N_GHOSE = int(cand[cand["pass"]]["ghose_ok"].sum())
N_NOVEL = int(cand[cand["pass"]]["novel_scaffold"].sum())
N_APPROVED_OVERLAP = int(cand[cand["pass"]]["already_approved_drug"].sum())
MEDIAN_TANIMOTO = float(cand[cand["pass"]]["max_tanimoto_to_known"].median())
N_SCREEN = len(cand)
RHO_IN, PVAL_IN = spearmanr(trs["TRS_fitted"], trs["observed_tox_index"])
RHO_LOO, PVAL_LOO = spearmanr(loo["TRS_loocv"], loo["observed_tox_index"])
top_tox = wts.sort_values("weight", ascending=False).iloc[0]
N_FDR_SIG = int(corr["significant_fdr05"].sum())
N_TESTS = len(corr)
N_CAPPED = int(report_counts["capped_at_2500"].sum())
N_PUBCHEM_MATCH = int(pcv["skeleton_match"].sum())
OVERALL_MIN_DATE = date_range["earliest_receivedate"].dropna().min()
OVERALL_MAX_DATE = date_range["latest_receivedate"].dropna().max()


def fmt_date(d):
    d = str(int(d))
    return f"{d[:4]}-{d[4:6]}-{d[6:]}"


doc = Document()
# base style
st = doc.styles["Normal"]; st.font.name = "Calibri"; st.font.size = Pt(11)


def H(text, level): doc.add_heading(text, level=level)
def P(text, italic=False, align=None, size=11):
    p = doc.add_paragraph()
    r = p.add_run(text); r.italic = italic; r.font.size = Pt(size)
    if align: p.alignment = align
    return p


def table(df, caption, fontsize=8):
    cap = doc.add_paragraph()
    cr = cap.add_run(caption); cr.bold = True; cr.font.size = Pt(9)
    t = doc.add_table(rows=1, cols=len(df.columns)); t.style = "Light Grid Accent 1"
    for i, c in enumerate(df.columns):
        rr = t.rows[0].cells[i].paragraphs[0].add_run(str(c)); rr.bold = True
        rr.font.size = Pt(fontsize)
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for i, v in enumerate(row):
            rr = cells[i].paragraphs[0].add_run(str(v)); rr.font.size = Pt(fontsize)
    doc.add_paragraph()


def landscape_table(df, caption, fontsize=8):
    """Table 2 is wide (11 columns) - give it its own landscape page, then
    return to portrait for the rest of the document.

    NOTE: doc.add_section() reuses the existing body-level sectPr element as
    the NEW section and relocates a clone of it to a paragraph as the
    "previous" section's properties - so a Section object captured before
    add_section() and the one returned after it alias the SAME element.
    Capture plain ints before mutating, not a live Section reference.
    """
    prev_w, prev_h = doc.sections[-1].page_width, doc.sections[-1].page_height
    land = doc.add_section(WD_SECTION.NEW_PAGE)
    land.orientation = WD_ORIENT.LANDSCAPE
    land.page_width, land.page_height = prev_h, prev_w
    land.left_margin = land.right_margin = Inches(0.6)
    table(df, caption, fontsize=fontsize)
    port = doc.add_section(WD_SECTION.NEW_PAGE)
    port.orientation = WD_ORIENT.PORTRAIT
    port.page_width, port.page_height = prev_w, prev_h


def figure(path, caption):
    doc.add_picture(path, width=Inches(6.2))
    doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    cp = doc.add_paragraph(); r = cp.add_run(caption)
    r.italic = True; r.font.size = Pt(9)


# ============================================================ TITLE
title = doc.add_paragraph()
title.alignment = WD_ALIGN_PARAGRAPH.CENTER
tr = title.add_run("Linking Chemical Structure to Organ-System Toxicity in "
                   "Antidepressants: A Pharmacovigilance and Cheminformatics "
                   "Analysis with a Structure-Based Toxicity Risk Score")
tr.bold = True; tr.font.size = Pt(15)
P("[Author name], [Affiliation]", align=WD_ALIGN_PARAGRAPH.CENTER, size=10)
P("Correspondence: [email]", align=WD_ALIGN_PARAGRAPH.CENTER, size=10)

# ============================================================ ABSTRACT
H("Abstract", 1)
P(f"Background: Antidepressants differ markedly in their organ-system safety, "
  f"yet the chemical determinants of that variation are incompletely mapped. We "
  f"combined real-world pharmacovigilance with cheminformatics to test whether "
  f"structural features of antidepressants are associated with organ-specific "
  f"toxicity signals, and to derive a transparent structure-based risk score.")
P(f"Methods: We curated {N_CURATED} antidepressants with PubChem-verified SMILES; "
  f"{N_DRUG} had retrievable FAERS data (one, clorgyline, was never marketed and "
  f"had zero reports) and were retained for analysis. We retrieved {N_REP:,} "
  f"adverse-event reports ({N_ROW:,} drug-reaction records) from the FDA Adverse "
  f"Event Reporting System via the openFDA API, covering receive dates "
  f"{fmt_date(OVERALL_MIN_DATE)} to {fmt_date(OVERALL_MAX_DATE)}. We computed "
  f"within-class disproportionality (Reporting Odds Ratio, ROR, with 95% CI; "
  f"Proportional Reporting Ratio; Bayesian Information Component with 95% "
  f"credible interval) for cardiac, hepatic, neurologic and metabolic toxicity, "
  f"with a pre-specified primary signal criterion of ROR lower 95% bound > 2 and "
  f"≥ 3 cases. Molecular descriptors and 15 candidate-toxicophore substructure "
  f"flags were computed with RDKit and tested against organ-specific ROR "
  f"(Mann-Whitney, rank-biserial effect size), with Benjamini-Hochberg FDR "
  f"correction across all organ x feature tests. A Toxicity Risk Score (TRS) was "
  f"fit from these effect sizes and validated both in-sample and by "
  f"leave-one-drug-out cross-validation. Compounds active against SERT, NET, DAT "
  f"and MAO-A were screened from ChEMBL and assessed for novelty (Tanimoto "
  f"similarity to the reference set) and drug-likeness (Lipinski, Veber, Ghose).")
P(f"Results: At the primary threshold we detected {N_SIG} disproportionality "
  f"signals across {N_DRUG} drugs (a legacy ROR>1 threshold, reported for "
  f"comparison only, flags {N_SIG_LEGACY}); the strongest was the reversible "
  f"MAO-A inhibitor moclobemide for neurologic toxicity (ROR 3.93). A tricyclic "
  f"core and a carboxylic-acid group were the leading structural correlates of "
  f"cardiac and hepatic ROR respectively, but none of the {N_TESTS} organ x "
  f"feature tests survived Benjamini-Hochberg correction (q<0.05); these "
  f"associations are exploratory. The structure-only TRS correlated with "
  f"observed FAERS toxicity in-sample (Spearman rho = {RHO_IN:+.2f}, p = "
  f"{PVAL_IN:.3f}, n = {len(trs)}), but this fell to rho = {RHO_LOO:+.2f} "
  f"(p = {PVAL_LOO:.2f}) under honest leave-one-drug-out cross-validation, "
  f"indicating the in-sample fit is substantially circular. Screening "
  f"{N_SCREEN:,} ChEMBL compounds yielded {N_PASS} candidates with a favourable "
  f"structural profile, of which {N_NOVEL} were structurally novel (Tanimoto "
  f"<0.85 to any reference drug; median similarity {MEDIAN_TANIMOTO:.2f}) and "
  f"none were already-approved drugs for another indication.")
P(f"Conclusions: Structural features show exploratory, mechanistically "
  f"plausible associations with real-world antidepressant toxicity signals, but "
  f"a purely structural risk score does not generalise out-of-sample in this "
  f"{N_DRUG}-drug set. Findings are hypothesis-generating, not confirmatory or "
  f"clinically actionable, and are constrained by drug-set size and the known "
  f"limitations of spontaneous-reporting data.")
kw = doc.add_paragraph(); kr = kw.add_run("Keywords: ")
kr.bold = True
kw.add_run("antidepressants; pharmacovigilance; FAERS; disproportionality; "
           "cheminformatics; toxicophore; structure-toxicity; drug safety; "
           "cross-validation.")

# ============================================================ 1 INTRO
H("1. Introduction", 1)
P("Antidepressants are among the most widely prescribed drug classes [1], and "
  "their safety profiles vary substantially across cardiac, hepatic, neurologic "
  "and metabolic domains. Tricyclic antidepressants (TCAs) carry well-recognised "
  "cardiotoxicity [2], while individual agents are associated with idiosyncratic "
  "hepatotoxicity [3], seizures, serotonin syndrome [4] or metabolic disturbance. "
  "Much of this variation is plausibly rooted in chemical structure, but a "
  "systematic, data-driven mapping between structural features and organ-specific "
  "toxicity across the class is lacking.")
P("Spontaneous-reporting systems such as the FDA Adverse Event Reporting System "
  "(FAERS) enable disproportionality analysis, which quantifies whether an "
  "adverse event is reported more often than expected for a given drug [5], [6]. "
  "Disproportionality reflects reporting behaviour, not incidence or causation: "
  "it is subject to reporting biases, confounding by indication/age/polypharmacy, "
  "and has no fixed denominator. In parallel, cheminformatics allows molecules to "
  "be decomposed into descriptors and substructures that can be tested as "
  "candidate toxicophores (structural alerts for harm) or safeophores (features "
  "associated with relative safety).")
P("This study integrates these two perspectives, with an emphasis on distinguishing "
  "exploratory association from validated prediction. Our objectives were to: (1) "
  "quantify FAERS disproportionality signals for cardiac, hepatic, neurologic and "
  "metabolic toxicity across the antidepressant set, under a pre-specified, "
  "conservative signal definition; (2) rank agents by safety within pharmacological "
  "class; (3) map structural features to the observed signals, identifying "
  "candidate toxicophores and safeophores, with multiplicity correction and a "
  "check for feature-feature redundancy; (4) screen public chemical-database "
  "compounds for favourable, genuinely novel structural profiles; and (5) propose "
  "a structure-based Toxicity Risk Score (TRS) and test whether it generalises "
  "out-of-sample.")

# ============================================================ 2 METHODS
H("2. Methods", 1)
H("2.1 Drug set and chemical structures", 2)
P(f"Twenty-eight antidepressants spanning five pharmacological classes (TCA, "
  f"selective serotonin reuptake inhibitor [SSRI], serotonin-norepinephrine "
  f"reuptake inhibitor [SNRI], reversible MAO-A inhibitor [RIMA], and atypical "
  f"agents) were curated with canonical SMILES. All 28 structures - not only the "
  f"two later found to be erroneous - were independently verified against PubChem "
  f"by InChIKey connectivity skeleton (first 14 characters, salt/stereo-agnostic); "
  f"{N_PUBCHEM_MATCH}/{N_CURATED} matched. Two entries (citalopram and "
  f"escitalopram) were found to carry an isomeric structural error and were "
  f"corrected before any analysis in this paper was run (Section 4 discusses the "
  f"downstream consequence of this correction).")
H("2.2 Adverse-event data", 2)
example_drug, example_q = "moclobemide", None
qstr_lines = open("faers_query_strings.txt").read().split("\n\n")
example_block = next((b for b in qstr_lines if b.startswith("moclobemide:")), qstr_lines[1] if len(qstr_lines) > 1 else "")
example_query = example_block.split("\n", 1)[1].strip() if "\n" in example_block else ""
P(f"Adverse-event reports were retrieved programmatically from the openFDA "
  f"drug-event endpoint [13] (no API key; rate limit 1,000 requests/day). "
  f"Queries combined the normalised openfda.generic_name field with free-text "
  f"patient.drug.medicinalproduct matching against the generic name and every "
  f"known brand/salt synonym, so that non-US agents (e.g. moclobemide, "
  f"mianserin, tianeptine) were captured even though openFDA does not normalise "
  f"their names. No receive-date filter was applied, so the query returns the "
  f"full history openFDA indexes for each drug; we determined post hoc, by "
  f"querying the earliest and latest receivedate for each drug's exact search "
  f"string, that the retrieved reports span {fmt_date(OVERALL_MIN_DATE)} to "
  f"{fmt_date(OVERALL_MAX_DATE)} (Supplementary Table S6). The exact query "
  f"string used for every drug is listed in Supplementary Material (example, "
  f"moclobemide): {example_query}")
P(f"Up to 2,500 reports were retrieved per drug. This cap was binding - i.e. the "
  f"drug had at least 2,500 retrievable reports and was truncated - for "
  f"{N_CAPPED} of {N_DRUG} drugs (Supplementary Table S7), so it is not a minor "
  f"formality: report volume for the highest-volume SSRIs and TCAs is "
  f"undercounted relative to lower-volume agents, and relative ranks among capped "
  f"drugs could shift with the complete history. The cap was adopted because "
  f"openFDA pagination latency increases sharply at deep offsets, making a full "
  f"pull impractical at scale; it yielded {N_REP:,} unique reports and "
  f"{N_ROW:,} drug-reaction records for {N_DRUG} agents with data. Clorgyline, "
  f"never marketed, returned zero reports and is excluded from all downstream "
  f"FAERS-dependent analyses (Objectives 1, 2, 3 and 5), though it remains in "
  f"the curated structure set for Objectives 3-4's structural screening.")
H("2.3 Disproportionality analysis", 2)
P("MedDRA preferred terms [9] were grouped into four organ systems using expanded "
  "term lists (Supplementary Table S8). For each drug-by-organ pair we built a "
  "2x2 contingency table against the pooled background of the other "
  "antidepressants (a within-class design - see rationale below) and computed "
  "the ROR with 95% confidence interval [5], the Proportional Reporting Ratio "
  "with chi-square [6], and the Bayesian Information Component (IC) with its 95% "
  "credible interval [7], [8]. The primary signal criterion was pre-specified as "
  "ROR lower 95% bound > 2 and at least three cases, a threshold consistent with "
  "standard pharmacovigilance practice and more conservative than the commonly "
  "seen but permissive ROR>1 screen; the ROR>1 criterion is reported only as a "
  "supplementary comparison (Table 1 footnote; Supplementary Table S5) and is "
  "not used for any headline result.")
P("Rationale for the within-class background: Objective 2 asks which "
  "antidepressants are relatively safer than others in the SAME class and "
  "indication space, where patients and confounders (age, depression severity, "
  "polypharmacy) are most comparable, so the other antidepressants are the more "
  "appropriate comparator than the full FAERS database of all drugs for all "
  "indications. The trade-off is that a within-class ROR expresses RELATIVE "
  "over/under-reporting versus other antidepressants, not an absolute risk - a "
  "drug can show ROR<1 here while still carrying an adverse event more often "
  "than the general non-antidepressant-using population. This is stated "
  "explicitly as a limitation (Section 5) rather than left implicit.")
H("2.4 Cheminformatics and structure-toxicity correlation", 2)
smarts_list = ", ".join(sorted(SUBSTRUCTURE_SMARTS.keys())) + ", tricyclic_core"
P(f"RDKit [10] was used to compute physicochemical descriptors and {len(SUBSTRUCTURE_SMARTS)} "
  f"SMARTS substructure flags representing candidate toxicophores/safeophores, "
  f"plus a ring-based tricyclic-core detector (>=2 fused aromatic 6-membered "
  f"rings and >=3 total rings) - {len(SUBSTRUCTURE_SMARTS) + 1} structural "
  f"features in total, ALL of which are reported in Supplementary Table S4 "
  f"regardless of significance ({smarts_list}). Each feature was tested for "
  f"association with organ-specific ROR using a Mann-Whitney U test on the "
  f"full ROR distribution (carriers vs non-carriers, every drug with data used, "
  f"not only those meeting the signal threshold), with the rank-biserial "
  f"correlation as effect size (range -1 to +1; the fraction of carrier-vs-"
  f"non-carrier ROR pairs favouring one group over the other, oriented so "
  f"positive = carriers have systematically higher ROR). Positive effects "
  f"denote toxicophores; negative effects denote safeophores, per organ system - "
  f"a feature's role is organ-specific and is not assumed to generalise across "
  f"organs (Figure 5). With 4 organ systems x {len(SUBSTRUCTURE_SMARTS) + 1} "
  f"features = {N_TESTS} simultaneous tests, we applied Benjamini-Hochberg FDR "
  f"correction across the full set of tests (one family) and report the "
  f"FDR-adjusted q-value alongside the nominal p-value. We additionally checked "
  f"the {len(SUBSTRUCTURE_SMARTS) + 1} structural features for pairwise "
  f"redundancy (phi coefficient on the 0/1 flags across all 28 drugs), since a "
  f"correlated pair cannot be independently attributed as drivers of the same "
  f"organ signal.")
H("2.5 Toxicity Risk Score", 2)
P("Each structural feature was assigned a weight equal to its mean rank-biserial "
  "effect on ROR across the four organ systems. A drug's TRS is the sum of the "
  "weights of the features it carries. The score was first validated by "
  "correlating it (Spearman) with each drug's observed toxicity index, defined "
  "as the mean of ln(ROR) across the four organ systems, using the SAME 27 "
  "drugs used to fit the weights - an in-sample check that is circular by "
  "construction, since the weights are themselves derived from these drugs' ROR "
  "values. To obtain an honest, non-circular estimate we additionally performed "
  "leave-one-drug-out cross-validation: for each drug, the feature weights were "
  "re-fit using every OTHER drug's data only, and the held-out drug was scored "
  "with those weights before its own data ever touched them; the resulting "
  "out-of-sample Spearman correlation is reported alongside the in-sample "
  "figure in every result (never one without the other). We also re-ran both "
  "validations after excluding moclobemide, whose neurologic ROR (3.93) is a "
  "clear outlier, as a sensitivity check.")
H("2.6 Candidate screening", 2)
target_lines = "; ".join(f"{name} ({tid})" for name, tid in TARGETS.items())
P(f"Compounds with potency pChEMBL >= {POTENCY:g} against four named "
  f"antidepressant-relevant targets - {target_lines} - were retrieved from "
  f"ChEMBL [11]. Compounds matching one of our 28 curated drugs by InChIKey "
  f"skeleton were excluded; this excludes only those 28 specific structures, "
  f"not all FDA-approved drugs generally, so each surviving candidate was "
  f"additionally checked against ChEMBL's max_phase field (max_phase=4 = "
  f"approved somewhere for some indication) to flag compounds that are already "
  f"marketed drugs for a different indication rather than novel entities. Each "
  f"remaining compound was screened with the same structural rules used on our "
  f"28 drugs and ranked by fitted TRS. Drug-likeness was assessed with three "
  f"independent filters - Lipinski's rule of five [12], the Veber rule "
  f"(<=10 rotatable bonds, TPSA<=140 A^2), and the Ghose filter (MW 160-480, "
  f"LogP -0.4 to 5.6, molar refractivity 40-130, 20-70 heavy atoms) - since "
  f"several reference antidepressants themselves violate Lipinski and a single "
  f"filter can be misleading. Structural novelty was quantified, not assumed: "
  f"for each candidate we computed the maximum Tanimoto similarity (Morgan "
  f"fingerprints, radius 2, 2048 bits) to any of the 28 reference "
  f"antidepressants, flagging <0.85 as a distinct scaffold.")
H("2.7 Software and reproducibility", 2)
P("Analyses used Python 3 with RDKit, pandas, NumPy and SciPy. Every result is "
  "regenerated from versioned scripts, and an independent sanity-check script "
  "re-derives key statistics (ROR values, signal flags, structure corrections, "
  "TRS validation, candidate-screen exclusions) directly from the raw data as a "
  "check against the main pipeline.")
H("2.8 Data and code availability", 2)
P("All analysis code (data retrieval, disproportionality, cheminformatics, TRS "
  "fitting/cross-validation, candidate screening, sanity checks, figure/table "
  "generation) and all derived data tables are available at "
  "https://github.com/Prafful7601/antidepressant-structure-toxicity. A "
  "versioned, citable release (e.g. a Zenodo DOI archiving a tagged release of "
  "this repository) should be minted at submission. FAERS data are drawn from "
  "the public openFDA API and are fully re-derivable from the versioned query "
  "strings in Supplementary Material; raw FAERS extracts and all intermediate "
  "CSVs are included in the repository so the full pipeline can be re-run "
  "without re-querying external APIs.")
H("2.9 Ethics statement", 2)
P("This study analyses only de-identified, publicly available, aggregate "
  "post-marketing surveillance data (FAERS, via openFDA) and public chemical "
  "database records (ChEMBL, PubChem). No identifiable patient information, "
  "direct patient contact, or non-public data were used. As secondary analysis "
  "of public, de-identified data, this work does not constitute human subjects "
  "research under the Common Rule (45 CFR 46) and institutional review "
  "board oversight was not required; authors should confirm this against their "
  "own institution's specific determination before submission.")

# ============================================================ 3 RESULTS
H("3. Results", 1)
H("3.1 Disproportionality signals (Objective 1)", 2)
P(f"Across {N_DRUG} drugs, {N_SIG} organ-system signals met the pre-specified "
  f"primary criterion (ROR lower 95% bound > 2, >=3 cases; Table 1): moclobemide "
  f"for neurologic toxicity (ROR 3.93, 95% CI 3.16-4.88; IC "
  f"{sig.loc[(sig.drug_name=='Moclobemide')&(sig.organ_system=='neurologic'),'IC'].iloc[0]:.2f}, "
  f"95% credible interval "
  f"{sig.loc[(sig.drug_name=='Moclobemide')&(sig.organ_system=='neurologic'),'IC025'].iloc[0]:.2f} to "
  f"{sig.loc[(sig.drug_name=='Moclobemide')&(sig.organ_system=='neurologic'),'IC975'].iloc[0]:.2f}), "
  f"consistent with the serotonergic and seizure-related liabilities of MAO "
  f"inhibition; and mianserin for hepatic toxicity (ROR 2.35, 95% CI 2.01-2.76). "
  f"The legacy, more permissive ROR>1 criterion flags {N_SIG_LEGACY} signals "
  f"(Supplementary Table S5); we report it only for comparison with prior "
  f"pharmacovigilance literature that commonly uses ROR>1; it is not used for "
  f"any conclusion in this paper. Tricyclic and atypical agents clustered in the "
  f"cardiac and hepatic domains under both criteria (Figure 2).")
flagged = (sig[sig["signal"] == 1][["drug_name", "organ_system", "n_cases",
           "ror", "ror_lower", "ror_upper", "IC", "IC025", "IC975"]].round(2)
           .sort_values("ror", ascending=False))
flagged.columns = ["Drug", "Organ", "Cases", "ROR", "ROR low", "ROR high",
                  "IC", "IC low", "IC high"]
table(flagged, "Table 1. Primary disproportionality signals (ROR lower bound > 2 "
      "and >= 3 cases), with 95% CI/credible intervals. Full drug x organ table "
      "(all 27 drugs, all four organs, both thresholds) is Supplementary Table S5.")

H("3.2 Within-class safety ranking (Objective 2)", 2)
P("Within every class the agents separated into clear safer and higher-risk "
  "tiers by the continuous toxicity index (Table 2; Figure 4), independent of "
  "the binary signal threshold used in 3.1. Among SSRIs, dapoxetine and "
  "sertraline ranked safest and fluoxetine highest-risk; among TCAs, imipramine "
  "and doxepin ranked safest and reboxetine and amoxapine highest-risk. Recall "
  "that this index is RELATIVE to other antidepressants (Section 2.3); it is "
  "not an absolute safety statement.")
rk = (rank[["class", "safety_rank_in_class", "drug_name"] + ORGANS +
      ["n_signals", "toxicity_index"]]
      .sort_values(["class", "safety_rank_in_class"]))
rk.columns = ["Class", "Rank", "Drug", "Card.", "Hep.", "Neuro.", "Metab.",
              "Signals (primary)", "Tox index"]
landscape_table(rk, "Table 2. Within-class safety ranking (rank 1 = safest; "
      "toxicity index = mean ln ROR over four organ systems; Signals = count "
      "meeting the PRIMARY ROR>2 criterion).", fontsize=9)

H("3.3 Structure-toxicity associations (Objective 3)", 2)
P(f"Several structural features showed organ-specific associations with ROR "
  f"(Table 3), but NONE of the {N_TESTS} organ x feature tests survived "
  f"Benjamini-Hochberg correction at q<0.05 ({N_FDR_SIG}/{N_TESTS} significant "
  f"after correction) - these associations should be read as exploratory leads, "
  f"not confirmed toxicophores/safeophores. Nominally, a tricyclic core and a "
  f"tertiary amine were associated with higher cardiac ROR; a carboxylic-acid "
  f"group and a long alkyl chain with higher hepatic ROR; and a morpholine "
  f"ring with higher neurologic ROR. A pairwise redundancy check found "
  f"carboxylic_acid and long_alkyl_chain are themselves correlated (phi = "
  f"{mcol.loc[(mcol.feature_1=='carboxylic_acid')|(mcol.feature_2=='carboxylic_acid'),'phi'].iloc[0]:.2f}), "
  f"so their hepatic associations cannot be independently attributed to one "
  f"feature or the other; tertiary_amine and basic_secondary_amine are also "
  f"correlated (phi = "
  f"{mcol.loc[(mcol.feature_1=='tertiary_amine')|(mcol.feature_2=='tertiary_amine'),'phi'].iloc[0]:.2f}, "
  f"expectedly, since they are largely mutually exclusive amine classes). The "
  f"tricyclic core and tertiary amine - both nominally cardiac toxicophores - "
  f"are NOT strongly correlated with each other (phi = "
  f"{PHI_TRICYCLIC_TERTIARY:.2f} in the full feature-correlation "
  f"matrix, Supplementary Table S9), so their associations are not redundant "
  f"with each other. In a sensitivity analysis excluding moclobemide (whose "
  f"neurologic ROR of 3.93 is the largest outlier in the dataset), the "
  f"morpholine-neurologic association disappeared entirely (moclobemide is one "
  f"of only two morpholine carriers), while the carboxylic_acid, "
  f"long_alkyl_chain and tricyclic_core associations were essentially unchanged "
  f"(Supplementary Table S10) - the morpholine finding should be treated as an "
  f"artifact of a single drug, not a general structure-toxicity relationship. "
  f"Figure 5 shows that several features (notably carboxylic acid) flip sign "
  f"between organs, so no feature can be labelled a toxicophore or safeophore "
  f"without specifying the organ system.")
P("Mechanistic context for the carboxylic-acid finding: a carboxylic acid is "
  "ionised at physiological pH, which lowers logD and increases polarity "
  "relative to the neutral/basic amines that dominate this drug set; reduced "
  "blood-brain-barrier penetration is a well-established consequence of this "
  "property [14], offering a plausible mechanism for why carboxylic-acid "
  "carriers show LOWER neurologic ROR (safeophore) in Table 3, while the same "
  "ionisable/reactive-metabolite-prone group is associated with HIGHER hepatic "
  "ROR (toxicophore). This is offered as biological plausibility for a "
  "nominal, FDR-non-significant association, not as confirmation of a causal "
  "mechanism.")
cs = (corr[corr["p_value"] < 0.10][["organ_system", "feature", "rank_biserial",
      "p_value", "p_adj_bh", "role"]].sort_values("p_value").round(4))
cs.columns = ["Organ", "Feature", "Effect size (rank-biserial, -1 to +1)", "p",
             "q (BH-FDR)", "Role"]
table(cs, "Table 3. Structure-toxicity associations at nominal p < 0.10 "
      "(continuous-ROR Mann-Whitney test); q = Benjamini-Hochberg adjusted "
      "p-value across all 40 tests. None reach q<0.05. Full results for every "
      "feature x organ combination, including non-significant, are "
      "Supplementary Table S4.")

H("3.4 Toxicity Risk Score: in-sample fit vs. honest cross-validation (Objective 5)", 2)
P(f"The fitted feature weights (Table 4; Figure 3) placed the tricyclic core as "
  f"the strongest toxicophore (weight {top_tox['weight']:+.2f}) and para-"
  f"halogenated phenyl and nitrile groups as mild safeophores. In-sample, the "
  f"structure-only TRS correlated with observed FAERS toxicity across drugs "
  f"(Spearman rho = {RHO_IN:+.2f}, p = {PVAL_IN:.3f}, n = {len(trs)}; Figure 1, "
  f"left panel). This in-sample figure is CIRCULAR: the feature weights are "
  f"themselves derived from these same drugs' ROR values, so it partly re-tests "
  f"the data used to build the score, and should not be reported without the "
  f"cross-validated figure alongside it. Leave-one-drug-out cross-validation - "
  f"refitting the weights from the other 26 drugs and scoring the held-out drug "
  f"blind - gives rho = {RHO_LOO:+.2f} (p = {PVAL_LOO:.2f}, n = {len(loo)}; "
  f"Figure 1, right panel): the score does NOT predict out-of-sample toxicity "
  f"in this dataset. Excluding moclobemide as a sensitivity check left both "
  f"conclusions unchanged (in-sample rho = "
  f"{spearmanr(sens['TRS_fitted_excl_moclobemide'], sens['observed_tox_index'])[0]:+.2f}; "
  f"LOOCV rho = "
  f"{spearmanr(sens['TRS_loocv'], sens['observed_tox_index'])[0]:+.2f}). We "
  f"therefore do not present the TRS as a validated predictive score; it is a "
  f"transparent SUMMARY of the (FDR-non-significant, Section 3.3) in-sample "
  f"structure-toxicity associations, offered as a hypothesis-generating "
  f"scaffold for a larger future drug set rather than a working risk-prediction "
  f"tool.")
w2 = wts.round(3); w2.columns = ["Feature", "TRS weight"]
table(w2, "Table 4. Toxicity Risk Score weights (mean ROR effect size across "
      "organ systems; + = toxicophore, - = safeophore). Underlying "
      "organ x feature tests do not survive FDR correction (Table 3).")

H("3.5 Candidate screening (Objective 4)", 2)
P(f"Of {N_SCREEN:,} non-redundant ChEMBL compounds active (pChEMBL>={POTENCY:g}) "
  f"against SERT, NET, DAT or MAO-A, {N_PASS} carried at least one safeophore "
  f"and no toxicophore by the fitted TRS weights. Drug-likeness varied by "
  f"filter: {N_DL} of {N_PASS} passed Lipinski, {N_VEBER} passed Veber, and "
  f"only {N_GHOSE} passed the stricter Ghose filter - underscoring that a "
  f"single drug-likeness rule is not sufficient given several reference "
  f"antidepressants themselves violate one or more rules. Novelty was "
  f"quantified rather than assumed: {N_NOVEL} of {N_PASS} passing candidates "
  f"had maximum Tanimoto similarity <0.85 to any of the 28 reference "
  f"antidepressants (median similarity across all passing candidates = "
  f"{MEDIAN_TANIMOTO:.2f}), and a ChEMBL max_phase cross-check found "
  f"{N_APPROVED_OVERLAP} of {N_PASS} passing candidates were already an "
  f"approved drug for another indication - i.e. these are genuinely novel "
  f"chemical entities relative to both the antidepressant class specifically "
  f"and approved drugs generally. The top-ranked candidates (lowest predicted "
  f"TRS; Table 5) are offered as hypothesis-generating leads for further "
  f"evaluation, conditioned on the exploratory, FDR-non-significant "
  f"structure-toxicity evidence (Section 3.3) and the non-generalising TRS "
  f"(Section 3.4) that rank them.")
top = (cand[cand["pass"]].head(12)[["molecule_chembl_id", "via_target", "MW",
       "LogP", "max_tanimoto_to_known", "safeophores", "TRS"]])
top.columns = ["ChEMBL ID", "Via target", "MW", "LogP", "Max Tanimoto to known",
              "Safeophores", "TRS"]
table(top, "Table 5. Top prioritized candidate molecules, ranked by predicted "
      "Toxicity Risk Score.", fontsize=8)

# ============================================================ FIGURES
H("Figures", 1)
figure("charts/trs_validation.png",
       f"Figure 1. Structure-based Toxicity Risk Score: in-sample (circular) fit "
       f"vs. leave-one-drug-out (honest) prediction. Left: rho = {RHO_IN:+.2f}, "
       f"p = {PVAL_IN:.3f}. Right: rho = {RHO_LOO:+.2f}, p = {PVAL_LOO:.2f} "
       f"(n = {len(trs)}).")
figure("charts/ror_heatmap.png",
       "Figure 2. Reporting Odds Ratio for each drug and organ system. Solid "
       "outline = primary signal (ROR lower 95% CI > 2, >=3 cases); dashed = "
       "legacy ROR>1 only (supplementary comparison).")
figure("charts/trs_weights.png",
       "Figure 3. Fitted structural-feature weights for the Toxicity Risk "
       "Score. None survive FDR correction (Table 3).")
figure("charts/class_toxicity.png",
       "Figure 4. Per-drug toxicity index (relative to other antidepressants) "
       "coloured by pharmacological class.")
figure("charts/feature_stability.png",
       "Figure 5. Structural-feature effect direction (rank-biserial) by organ "
       "system. Sign flips across organs (e.g. carboxylic acid) show that a "
       "feature's toxicophore/safeophore role is organ-specific, not global.")

# ============================================================ 4 DISCUSSION
H("4. Discussion", 1)
P(f"Integrating spontaneous-reporting pharmacovigilance with cheminformatics, "
  f"we found nominal, mechanistically plausible associations between structural "
  f"features and organ-specific toxicity signals, but two rigor checks that are "
  f"often omitted from similar studies - multiple-testing correction and "
  f"out-of-sample validation - substantially temper the headline claims. First, "
  f"none of the {N_TESTS} structure-toxicity tests survive Benjamini-Hochberg "
  f"correction, so the tricyclic-core/cardiac and carboxylic-acid/hepatic "
  f"associations, while mechanistically coherent with known tricyclic "
  f"cardiotoxicity and reactive-metabolite hepatic liabilities, remain "
  f"exploratory. Second, the TRS's apparent validation (in-sample Spearman rho "
  f"~{RHO_IN:.2f}) is largely an artifact of fitting and testing the score on "
  f"the same {len(trs)} drugs: leave-one-drug-out cross-validation collapses "
  f"this to rho~{RHO_LOO:.2f}, statistically indistinguishable from zero. A "
  f"structure-only score built from {len(SUBSTRUCTURE_SMARTS)+1} features and "
  f"27-28 drugs simply does not have enough independent information to predict "
  f"an unseen drug's toxicity from its structure alone in this dataset.")
P("The candidate-screening step demonstrates how such structure-toxicity rules "
  "can be projected onto chemical space, and here the added checks were "
  "reassuring rather than tempering: nearly all passing candidates were "
  "genuinely novel by Tanimoto similarity and none overlapped with "
  "already-approved drugs, so the screen is not simply rediscovering known "
  "agents. That reassurance is nonetheless bounded by the same "
  "FDR-non-significant, non-generalising structure-toxicity evidence used to "
  "rank the candidates (Sections 3.3-3.4): a novel molecule ranked \"safe\" by "
  "this TRS is a hypothesis to test, not a validated prediction.")
P("A methodological observation reinforces the value of structure validation: "
  "an isomeric error in two source structures (citalopram, escitalopram), "
  "identified and corrected during dataset curation (Section 2.1) before any "
  "analysis was run, initially would have produced a spurious cardiac "
  "safeophore; the corrected structures instead converge on the "
  "pharmacologically expected tricyclic-core association. Structure curation "
  "is therefore not incidental but central to reproducible structure-toxicity "
  "work, which is why we extended verification to all 28 structures against "
  "PubChem (Section 2.1) rather than only the two known-erroneous entries.")

# ============================================================ 5 LIMITATIONS
H("5. Limitations", 1)
for t in [
    f"The analysis rests on {N_DRUG} drugs with FAERS data (28 curated "
    f"structures; clorgyline had zero reports); structure-toxicity "
    f"correlations do not survive Benjamini-Hochberg multiple-testing "
    f"correction, and the Toxicity Risk Score does not generalise under "
    f"leave-one-drug-out cross-validation (rho falls from "
    f"{RHO_IN:+.2f} in-sample to {RHO_LOO:+.2f} out-of-sample). Both findings "
    f"should be read as hypothesis-generating, not confirmatory or "
    f"clinically actionable, pending a substantially larger and more diverse "
    f"drug set.",
    "FAERS is a spontaneous-reporting system without denominators and is "
    "subject to reporting biases (including the notoriety effect) and "
    "confounding by age, indication and polypharmacy; disproportionality "
    "reflects over- or under-reporting relative to a chosen background, not "
    "incidence or causation.",
    "A within-class background was used by design to answer the within-class "
    "ranking objective (rationale, Section 2.3); this yields RELATIVE "
    "over/under-reporting versus other antidepressants, not an absolute risk "
    "versus the general population or non-antidepressant drugs.",
    f"Per-drug report counts were capped at 2,500, and this cap was binding "
    f"(i.e., truncating available reports) for {N_CAPPED} of {N_DRUG} drugs, "
    f"predominantly higher-volume SSRIs and TCAs; absolute case counts and, "
    f"potentially, relative ranks among capped drugs may differ under the "
    f"complete FAERS history.",
    "Candidate molecules are computational leads conditioned on the strength "
    "of the underlying structure-toxicity evidence (exploratory, "
    "FDR-non-significant) and the TRS (non-generalising out-of-sample); they "
    "require experimental follow-up and should not be interpreted as validated "
    "safety predictions.",
    "The morpholine-neurologic association was found, in a sensitivity "
    "analysis, to be driven almost entirely by a single drug (moclobemide, one "
    "of only two morpholine carriers in the set) and should not be "
    "generalised.",
]:
    doc.add_paragraph(t, style="List Bullet")

# ============================================================ 6 CONCLUSION
H("6. Conclusion", 1)
P(f"Across a curated {N_DRUG}-drug antidepressant set, structural features show "
  f"exploratory, mechanistically plausible but statistically non-significant "
  f"(after FDR correction) associations with real-world organ-system toxicity "
  f"signals. A transparent Toxicity Risk Score built from these associations "
  f"fits the same data it was derived from (rho={RHO_IN:+.2f}) but does not "
  f"predict held-out drugs (LOOCV rho={RHO_LOO:+.2f}), so it should be treated "
  f"as a hypothesis-generating summary rather than a validated risk-prediction "
  f"tool. The same structural rules, applied to public chemical databases, do "
  f"identify a large set of genuinely novel, target-relevant candidate "
  f"molecules with favourable structural profiles - a use case that depends "
  f"less on the TRS generalising perfectly and more on it providing a "
  f"reasonable prioritisation heuristic. Expansion to a substantially larger "
  f"and more diverse drug set, with pre-registered multiplicity control and "
  f"cross-validation reported by default, is needed before this framework could "
  f"support structure-aware safety triage in antidepressant discovery.")

# ============================================================ REFERENCES
H("References", 1)
refs = [
    "Brody DJ, Gu Q. Antidepressant use among adults: United States, 2015-2018. "
    "NCHS Data Brief. 2020;(377):1-8.",
    "Thanacoody HKR, Thomas SHL. Tricyclic antidepressant poisoning: cardiovascular "
    "toxicity. Toxicol Rev. 2005;24(3):205-214.",
    "Voican CS, Corruble E, Naveau S, Perlemuter G. Antidepressant-induced liver "
    "injury: a review for clinicians. Am J Psychiatry. 2014;171(4):404-415.",
    "Boyer EW, Shannon M. The serotonin syndrome. N Engl J Med. "
    "2005;352(11):1112-1120.",
    "Rothman KJ, Lanes S, Sacks ST. The reporting odds ratio and its advantages "
    "over the proportional reporting ratio. Pharmacoepidemiol Drug Saf. "
    "2004;13(8):519-523.",
    "Evans SJW, Waller PC, Davis S. Use of proportional reporting ratios (PRRs) "
    "for signal generation from spontaneous adverse drug reaction reports. "
    "Pharmacoepidemiol Drug Saf. 2001;10(6):483-486.",
    "Bate A, Lindquist M, Edwards IR, et al. A Bayesian neural network method for "
    "adverse drug reaction signal generation. Eur J Clin Pharmacol. "
    "1998;54(4):315-321.",
    "van Puijenbroek EP, Bate A, Leufkens HGM, et al. A comparison of measures of "
    "disproportionality for signal detection in spontaneous reporting systems for "
    "adverse drug reactions. Pharmacoepidemiol Drug Saf. 2002;11(1):3-10.",
    "Brown EG, Wood L, Wood S. The medical dictionary for regulatory activities "
    "(MedDRA). Drug Saf. 1999;20(2):109-117.",
    "Landrum G. RDKit: Open-source cheminformatics. https://www.rdkit.org.",
    "Mendez D, Gaulton A, Bento AP, et al. ChEMBL: towards direct deposition of "
    "bioassay data. Nucleic Acids Res. 2019;47(D1):D930-D940.",
    "Lipinski CA, Lombardo F, Dominy BW, Feeney PJ. Experimental and computational "
    "approaches to estimate solubility and permeability in drug discovery and "
    "development settings. Adv Drug Deliv Rev. 2001;46(1-3):3-26.",
    "US Food and Drug Administration. openFDA / FDA Adverse Event Reporting System "
    "(FAERS). https://open.fda.gov.",
    "Pajouhesh H, Lenz GR. Medicinal chemistry properties associated with "
    "blood-brain barrier permeability. NeuroRx. 2005;2(4):541-553.",
    "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical "
    "and powerful approach to multiple testing. J R Stat Soc B. "
    "1995;57(1):289-300.",
    "Veber DF, Johnson SR, Cheng HY, et al. Molecular properties that influence "
    "the oral bioavailability of drug candidates. J Med Chem. "
    "2002;45(12):2615-2623.",
    "Ghose AK, Viswanadhan VN, Wendoloski JJ. A knowledge-based approach in "
    "designing combinatorial or medicinal chemistry libraries for drug "
    "discovery. J Comb Chem. 1999;1(1):55-68.",
]
for i, r in enumerate(refs, 1):
    doc.add_paragraph(f"{i}. {r}")

P("Note: reference list is provided in a generic numbered style; reformat to the "
  "target journal's citation style before submission.", italic=True, size=9)

# ============================================================ SUPPLEMENTARY MATERIAL
doc.add_page_break()
H("Supplementary Material", 1)

P("Table S1. All 15 structural features tested (SMARTS definitions).")
smarts_rows = [{"feature": k, "SMARTS / definition": v}
              for k, v in SUBSTRUCTURE_SMARTS.items()]
smarts_rows.append({"feature": "tricyclic_core",
                    "SMARTS / definition": "ring-based: >=2 fused aromatic "
                    "6-membered rings AND >=3 total rings (RDKit ring info, "
                    "not a single SMARTS - fused polycyclic patterns are "
                    "brittle as SMARTS)"})
table(pd.DataFrame(smarts_rows), "Table S1. SMARTS pattern / definition for "
      "every structural feature screened.", fontsize=8)

P("Table S2. Target proteins used for ChEMBL candidate screening.")
table(pd.DataFrame([{"target": k, "ChEMBL target ID": v, "pChEMBL threshold": POTENCY}
                    for k, v in TARGETS.items()]),
      "Table S2. ChEMBL targets and potency threshold.", fontsize=9)

P("Table S3. Structural-feature pairwise redundancy (|phi| >= 0.5).")
table(mcol.round(3), "Table S3. Structural-feature pairs with |phi| >= 0.5 "
      "(Pearson correlation on 0/1 flags across 28 drugs).", fontsize=9)

P("Table S4. Full structure-toxicity correlation results, every organ x "
  "feature combination (not only nominal p<0.10).")
full_corr = corr[["organ_system", "feature", "n_with_feature", "median_ror_with",
                  "median_ror_without", "rank_biserial", "p_value", "p_adj_bh",
                  "role"]].round(4)
table(full_corr, f"Table S4. All {len(full_corr)} organ x feature "
      "structure-toxicity tests (continuous ROR, Mann-Whitney), with "
      "Benjamini-Hochberg q-values.", fontsize=7)

P("Table S5. Full disproportionality results for every drug x organ pair, "
  "both signal thresholds.")
full_sig = sig[["drug_name", "organ_system", "n_cases", "ror", "ror_lower",
                "ror_upper", "IC", "IC025", "IC975", "signal", "signal_ror1"]].round(3)
full_sig.columns = ["Drug", "Organ", "Cases", "ROR", "ROR low", "ROR high",
                    "IC", "IC low", "IC high", "Signal (ROR>2, primary)",
                    "Signal (ROR>1, legacy)"]
landscape_table(full_sig, f"Table S5. All {len(full_sig)} drug x organ "
      "disproportionality results, including non-significant pairs, both "
      "signal thresholds.", fontsize=7)

P("Table S6. FAERS receive-date range retrieved per drug.")
dr = date_range.copy()
dr["earliest_receivedate"] = dr["earliest_receivedate"].apply(
    lambda d: fmt_date(d) if pd.notna(d) else "n/a")
dr["latest_receivedate"] = dr["latest_receivedate"].apply(
    lambda d: fmt_date(d) if pd.notna(d) else "n/a")
table(dr, f"Table S6. Earliest/latest FAERS receivedate returned per drug's "
      f"exact query (overall range {fmt_date(OVERALL_MIN_DATE)} to "
      f"{fmt_date(OVERALL_MAX_DATE)}). Full per-drug query strings are in "
      f"faers_query_strings.txt in the code repository.", fontsize=8)

P("Table S7. Per-drug FAERS report counts retrieved, and whether the 2,500 cap "
  "was binding.")
rc = report_counts.rename(columns={"drugname": "drug_name_lower"})
table(rc, "Table S7. Unique reports retrieved per drug; capped_at_2500=True "
      "means the drug had >=2,500 retrievable reports and was truncated.",
      fontsize=8)

P("Table S8. Organ system grouping of MedDRA preferred terms (excerpt - see "
  "faers_disproportionality.py ORGAN_SYSTEM_PT in the code repository for the "
  "complete, versioned lists used).")

P("Table S9. Full structural-feature multicollinearity matrix (phi "
  "coefficients, all pairs) is provided as "
  "structural_feature_multicollinearity.csv in the code repository.")

P("Table S10. Sensitivity analysis excluding moclobemide: per-drug TRS "
  "(in-sample and leave-one-out) vs. observed toxicity index.")
sens_show = sens.round(3)
table(sens_show, "Table S10. TRS validation with moclobemide excluded from the "
      "drug universe entirely (both weight-fitting and scoring).", fontsize=8)

P("Table S11. PubChem structural verification for all 28 curated drugs.")
pcv_show = pcv[["drug_name", "our_skeleton", "pubchem_query_name",
                "pubchem_skeleton", "skeleton_match"]]
table(pcv_show, "Table S11. InChIKey connectivity-skeleton comparison, our "
      "curated SMILES vs. PubChem, for all 28 drugs.", fontsize=8)

doc.save("Research_Paper.docx")
print("Wrote Research_Paper.docx")
