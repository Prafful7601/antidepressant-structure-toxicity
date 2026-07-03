"""Generate Research_Paper_IEEE.docx: IEEE two-column formatted manuscript."""
import pandas as pd
from scipy.stats import spearmanr
from docx import Document
from docx.shared import Inches, Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_SECTION
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

from chembl_candidate_screen import TARGETS, POTENCY
from cheminformatics_analysis import SUBSTRUCTURE_SMARTS
N_FEATURES = len(SUBSTRUCTURE_SMARTS) + 1

ORGANS = ["cardiac", "hepatic", "neurologic", "metabolic"]
FONT = "Times New Roman"

sig = pd.read_csv("faers_signals.csv")
rank = pd.read_csv("class_safety_ranking.csv")
trs = pd.read_csv("drug_trs.csv")
loo = pd.read_csv("trs_loocv.csv")
wts = pd.read_csv("trs_weights.csv")
cand = pd.read_csv("chembl_candidates_screened.csv")
corr = pd.read_csv("structure_toxicity_correlation_continuous.csv")
rep = pd.read_csv("faers_reports.csv", dtype=str)
date_range = pd.read_csv("faers_date_range.csv")
ds = pd.read_csv("antidepressants_dataset.csv")
pcv = pd.read_csv("pubchem_verification.csv")
report_counts = pd.read_csv("faers_report_counts_per_drug.csv")

N_REP = rep["primaryid"].nunique(); N_ROW = len(rep)
N_DRUG = rep["drugname"].nunique()
N_CURATED = len(ds)
N_SIG = int(sig["signal"].sum())               # primary (ROR lower CI > 2)
N_SIG_LEGACY = int(sig["signal_ror1"].sum())   # legacy ROR>1
N_PASS = int(cand["pass"].sum()); N_DL = int(cand[cand["pass"]]["lipinski_ok"].sum())
N_NOVEL = int(cand[cand["pass"]]["novel_scaffold"].sum())
N_APPROVED_OVERLAP = int(cand[cand["pass"]]["already_approved_drug"].sum())
N_SCREEN = len(cand)
N_TESTS = len(corr)
N_FDR_SIG = int(corr["significant_fdr05"].sum())
N_CAPPED = int(report_counts["capped_at_2500"].sum())
N_PUBCHEM_MATCH = int(pcv["skeleton_match"].sum())
RHO_IN, PVAL_IN = spearmanr(trs["TRS_fitted"], trs["observed_tox_index"])
RHO_LOO, PVAL_LOO = spearmanr(loo["TRS_loocv"], loo["observed_tox_index"])
top_tox = wts.sort_values("weight", ascending=False).iloc[0]
OVERALL_MIN_DATE = date_range["earliest_receivedate"].dropna().min()
OVERALL_MAX_DATE = date_range["latest_receivedate"].dropna().max()


def fmt_date(d):
    d = str(int(d))
    return f"{d[:4]}-{d[4:6]}-{d[6:]}"

doc = Document()
normal = doc.styles["Normal"]
normal.font.name = FONT; normal.font.size = Pt(10)


def set_cols(section, n, space_twips=360):
    sectPr = section._sectPr
    cols = sectPr.find(qn("w:cols"))
    if cols is None:
        cols = OxmlElement("w:cols"); sectPr.append(cols)
    cols.set(qn("w:num"), str(n)); cols.set(qn("w:space"), str(space_twips))


def set_margins(section):
    section.top_margin = Inches(0.75); section.bottom_margin = Inches(1.0)
    section.left_margin = Inches(0.625); section.right_margin = Inches(0.625)


def newsec(ncols):
    s = doc.add_section(WD_SECTION.CONTINUOUS)
    set_margins(s); set_cols(s, ncols)
    return s


def body(text, justify=True, size=10, italic=False, bold=False, space_after=6):
    p = doc.add_paragraph()
    r = p.add_run(text); r.font.size = Pt(size); r.font.name = FONT
    r.italic = italic; r.bold = bold
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY if justify else WD_ALIGN_PARAGRAPH.LEFT
    p.paragraph_format.space_after = Pt(space_after)
    return p


def section_head(text):  # IEEE level-1: centered, small caps look (uppercase)
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(text); r.font.size = Pt(10); r.font.name = FONT; r.bold = False
    p.paragraph_format.space_before = Pt(10); p.paragraph_format.space_after = Pt(4)


def sub_head(text):  # IEEE level-2: italic, left
    p = doc.add_paragraph()
    r = p.add_run(text); r.italic = True; r.font.size = Pt(10); r.font.name = FONT
    p.paragraph_format.space_before = Pt(4); p.paragraph_format.space_after = Pt(2)


def caption(text, size=8, center=True):
    p = doc.add_paragraph()
    if center: p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(text); r.font.size = Pt(size); r.font.name = FONT
    return p


def ieee_table(df, label, title, fontsize=7):
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(label); r.bold = False; r.font.size = Pt(8); r.font.name = FONT
    p2 = doc.add_paragraph(); p2.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r2 = p2.add_run(title); r2.font.size = Pt(8); r2.font.name = FONT
    p2.paragraph_format.space_after = Pt(2)
    t = doc.add_table(rows=1, cols=len(df.columns)); t.style = "Table Grid"
    t.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for i, c in enumerate(df.columns):
        rr = t.rows[0].cells[i].paragraphs[0].add_run(str(c))
        rr.bold = True; rr.font.size = Pt(fontsize); rr.font.name = FONT
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for i, v in enumerate(row):
            rr = cells[i].paragraphs[0].add_run(str(v))
            rr.font.size = Pt(fontsize); rr.font.name = FONT
    doc.add_paragraph().paragraph_format.space_after = Pt(4)


def ieee_figure(path, cap, width=3.3):
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run().add_picture(path, width=Inches(width))
    caption(cap, size=8)


# ===================== SECTION 1: full-width title block =====================
set_margins(doc.sections[0]); set_cols(doc.sections[0], 1)
t = doc.add_paragraph(); t.alignment = WD_ALIGN_PARAGRAPH.CENTER
tr = t.add_run("Linking Chemical Structure to Organ-System Toxicity in "
               "Antidepressants: A Pharmacovigilance and Cheminformatics Analysis "
               "with a Structure-Based Toxicity Risk Score")
tr.bold = True; tr.font.size = Pt(20); tr.font.name = FONT
for line, sz in [("[Author Name]", 11), ("[Department, Institution]", 10),
                 ("[email address]", 10)]:
    pp = doc.add_paragraph(); pp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    rr = pp.add_run(line); rr.font.size = Pt(sz); rr.font.name = FONT

# ===================== SECTION 2: two-column body ============================
newsec(2)

# Abstract (IEEE: bold "Abstract—", italic body, 9pt)
ap = doc.add_paragraph()
ar = ap.add_run("Abstract—"); ar.bold = True; ar.italic = True
ar.font.size = Pt(9); ar.font.name = FONT
abs = (f"Antidepressants differ markedly in organ-system safety, yet the chemical "
       f"determinants are incompletely mapped. We combined real-world "
       f"pharmacovigilance with cheminformatics to test whether structural features "
       f"are associated with organ-specific toxicity, and to derive a transparent "
       f"structure-based risk score. We curated {N_CURATED} antidepressants "
       f"(PubChem-verified, {N_PUBCHEM_MATCH}/{N_CURATED} matched); {N_DRUG} had "
       f"retrievable FAERS data (clorgyline, never marketed, had none). We "
       f"retrieved {N_REP:,} adverse-event reports ({N_ROW:,} drug-reaction "
       f"records) from FAERS via openFDA (receive dates "
       f"{fmt_date(OVERALL_MIN_DATE)} to {fmt_date(OVERALL_MAX_DATE)}), and "
       f"computed within-class disproportionality (ROR with 95% CI, PRR, "
       f"Information Component with 95% credible interval) with a pre-specified "
       f"primary signal criterion of ROR lower 95% bound > 2 and >=3 cases. RDKit "
       f"substructure flags were tested against ROR (Mann-Whitney, rank-biserial "
       f"effect size) with Benjamini-Hochberg FDR correction across all "
       f"organ x feature tests; a Toxicity Risk Score (TRS) was fit and validated "
       f"both in-sample and by leave-one-drug-out cross-validation; and ChEMBL "
       f"compounds were screened and checked for novelty (Tanimoto) and "
       f"drug-likeness (Lipinski, Veber, Ghose). We detected {N_SIG} primary "
       f"signals across {N_DRUG} drugs (legacy ROR>1 threshold: {N_SIG_LEGACY}, "
       f"reported for comparison only). None of {N_TESTS} organ x feature tests "
       f"survived FDR correction (q<0.05); nominally, a tricyclic core and a "
       f"carboxylic-acid group were the leading cardiac/hepatic correlates. The "
       f"structure-only TRS correlated with observed toxicity in-sample "
       f"(Spearman ρ = {RHO_IN:+.2f}, p = {PVAL_IN:.3f}, n = {len(trs)}) but fell "
       f"to ρ = {RHO_LOO:+.2f} (p = {PVAL_LOO:.2f}) under honest leave-one-out "
       f"cross-validation. Screening {N_SCREEN:,} compounds yielded {N_PASS} "
       f"candidates with a favourable profile ({N_DL} Lipinski-ok, {N_NOVEL} "
       f"genuinely novel by Tanimoto similarity, {N_APPROVED_OVERLAP} already an "
       f"approved drug elsewhere). Structural features show exploratory, "
       f"mechanistically plausible but statistically non-significant "
       f"associations with real-world toxicity; the TRS does not generalise "
       f"out-of-sample; findings are hypothesis-generating, not confirmatory.")
ar2 = ap.add_run(abs); ar2.italic = True; ar2.font.size = Pt(9); ar2.font.name = FONT
ap.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

ip = doc.add_paragraph()
ir = ip.add_run("Index Terms—"); ir.bold = True; ir.italic = True
ir.font.size = Pt(9); ir.font.name = FONT
ir2 = ip.add_run("antidepressants, pharmacovigilance, FAERS, disproportionality, "
                 "cheminformatics, toxicophore, drug safety.")
ir2.italic = True; ir2.font.size = Pt(9); ir2.font.name = FONT

# ---- I. INTRODUCTION ----
section_head("I.  Introduction")
body("Antidepressants are among the most widely prescribed drug classes [1], and "
     "their safety profiles vary substantially across cardiac, hepatic, neurologic "
     "and metabolic domains. Tricyclic antidepressants (TCAs) carry recognised "
     "cardiotoxicity [2], while individual agents are linked to idiosyncratic "
     "hepatotoxicity [3], seizures, serotonin syndrome [4] or metabolic "
     "disturbance. Much of this variation is plausibly rooted in chemical "
     "structure, yet a systematic, data-driven mapping between structural features "
     "and organ-specific toxicity across the class is lacking.")
body("Spontaneous-reporting systems such as the FDA Adverse Event Reporting System "
     "(FAERS) enable disproportionality analysis, which quantifies whether an event "
     "is reported more often than expected for a drug [5], [6]. In parallel, "
     "cheminformatics decomposes molecules into descriptors and substructures that "
     "can be tested as candidate toxicophores (alerts for harm) or safeophores "
     "(features associated with relative safety).")
body("This study integrates both perspectives. Our objectives were to: (1) "
     "quantify FAERS disproportionality signals for cardiac, hepatic, neurologic "
     "and metabolic toxicity; (2) rank agents by safety within pharmacological "
     "class; (3) map structural features to the signals, identifying toxicophores "
     "and safeophores; (4) screen public chemical-database compounds for favourable "
     "structural profiles; and (5) propose and validate a structure-based Toxicity "
     "Risk Score (TRS).")

# ---- II. METHODS ----
section_head("II.  Methods")
sub_head("A.  Drug Set and Chemical Structures")
body(f"Twenty-eight antidepressants spanning five classes (TCA, SSRI, SNRI, RIMA, "
     f"atypical) were curated with canonical SMILES. All 28 structures (not only "
     f"the two later found erroneous) were independently verified against "
     f"PubChem by InChIKey connectivity skeleton; {N_PUBCHEM_MATCH}/{N_CURATED} "
     f"matched. Two entries (citalopram, escitalopram) carried an isomeric error "
     f"and were corrected before analysis.")
sub_head("B.  Adverse-Event Data")
body(f"Reports were retrieved from the openFDA drug-event endpoint [13]. To capture "
     f"non-US agents, queries combined the normalised generic-name field with "
     f"free-text product- and brand-name matching; no date filter was applied, "
     f"and we determined post hoc that retrieved reports span "
     f"{fmt_date(OVERALL_MIN_DATE)} to {fmt_date(OVERALL_MAX_DATE)} (exact "
     f"per-drug query strings in Supplementary Material). Up to 2,500 reports "
     f"were retrieved per drug; this cap was binding for {N_CAPPED} of {N_DRUG} "
     f"drugs (mostly high-volume SSRIs/TCAs), yielding {N_REP:,} unique reports "
     f"and {N_ROW:,} drug-reaction records for {N_DRUG} agents with data "
     f"(clorgyline, never marketed, had none).")
sub_head("C.  Disproportionality Analysis")
body("MedDRA preferred terms [9] were grouped into four organ systems. For each "
     "drug-by-organ pair we built a 2x2 table against the pooled background of the "
     "other antidepressants - a within-class design, chosen because Objective 2 "
     "compares agents within the same indication/confounder space, at the cost of "
     "expressing RELATIVE, not absolute, over/under-reporting - and computed the "
     "ROR with 95% CI [5], PRR [6], and the Bayesian Information Component with "
     "95% credible interval [7], [8]. The pre-specified PRIMARY signal criterion "
     "was ROR lower 95% bound > 2 and >=3 cases (standard practice); the more "
     "permissive legacy ROR>1 criterion is reported only for comparison.")
sub_head("D.  Cheminformatics and Structure-Toxicity Correlation")
body(f"RDKit [10] computed physicochemical descriptors and 15 SMARTS substructure "
     f"flags (all reported regardless of significance, Supplementary Table S1) "
     f"plus a ring-based tricyclic-core detector. Each of 16 features was tested "
     f"against organ-specific ROR (every drug with data, not only signals) using "
     f"a Mann-Whitney test with rank-biserial effect size, per organ (a feature's "
     f"role is organ-specific, not assumed global - Fig. 5). With 4 organs x 16 "
     f"features = {N_TESTS} tests, Benjamini-Hochberg FDR correction was applied "
     f"across the full family. Features were additionally checked pairwise for "
     f"redundancy (phi coefficient on 0/1 flags).")
sub_head("E.  Toxicity Risk Score")
body("Each feature was weighted by its mean rank-biserial effect on ROR across the "
     "four organ systems; a drug's TRS is the sum of weights of features it "
     "carries. Weights are fit on the same drugs used for in-sample validation, "
     "which is circular by construction; we additionally performed leave-one-"
     "drug-out cross-validation (weights refit on the other 26 drugs, held-out "
     "drug scored blind) for an honest, non-circular estimate, and repeated both "
     "after excluding moclobemide (neurologic ROR 3.93, a clear outlier) as a "
     "sensitivity check.")
sub_head("F.  Candidate Screening")
target_list = ", ".join(f"{k} ({v})" for k, v in TARGETS.items())
body(f"Compounds with potency pChEMBL >= {POTENCY:g} against four named targets "
     f"({target_list}) were retrieved from ChEMBL [11]. Compounds matching our 28 "
     f"drugs by InChIKey skeleton were excluded (only those 28 structures, not "
     f"all approved drugs generally - remaining candidates were additionally "
     f"checked against ChEMBL max_phase to flag drugs already approved for other "
     f"indications). Each candidate was screened with the same rules, ranked by "
     f"TRS, assessed for drug-likeness by Lipinski [12], Veber, and Ghose, and "
     f"scored for structural novelty (max Tanimoto similarity, Morgan "
     f"fingerprints, to the 28 reference drugs).")
sub_head("G.  Data, Code and Ethics")
body("All code and derived data are intended for a public, versioned repository "
     "(e.g. GitHub/Zenodo) at publication [Repository URL / DOI to be inserted]. "
     "This study uses only de-identified, public aggregate data (FAERS, ChEMBL, "
     "PubChem) and does not constitute human subjects research; authors should "
     "confirm this against their institution's own determination.")

# ---- III. RESULTS ----
section_head("III.  Results")
sub_head("A.  Disproportionality Signals")
body(f"Across {N_DRUG} drugs, {N_SIG} signals met the primary criterion (Table I): "
     f"moclobemide, neurologic (ROR 3.93; 691 cases), consistent with the "
     f"liabilities of MAO inhibition; and mianserin, hepatic (ROR 2.35). The "
     f"legacy ROR>1 criterion flags {N_SIG_LEGACY} (Supplementary Table S5, "
     f"comparison only). Tricyclics clustered in cardiac/hepatic domains under "
     f"both criteria (Fig. 2).")
sub_head("B.  Within-Class Safety Ranking")
body("Within every class agents separated into safer and higher-risk tiers by the "
     "continuous toxicity index (Table II; Fig. 4), independent of the binary "
     "threshold above. Among SSRIs, dapoxetine and sertraline ranked safest and "
     "fluoxetine highest-risk; among TCAs, imipramine and doxepin ranked safest, "
     "reboxetine and amoxapine highest-risk. This index is relative to other "
     "antidepressants, not an absolute risk.")
sub_head("C.  Structure-Toxicity Associations")
body(f"A tricyclic core and a tertiary amine nominally associated with cardiac "
     f"ROR; a carboxylic-acid group and a long alkyl chain with hepatic ROR; a "
     f"morpholine ring with neurologic ROR (Table III) - but NONE of {N_TESTS} "
     f"organ x feature tests survived FDR correction ({N_FDR_SIG}/{N_TESTS} at "
     f"q<0.05): exploratory only. carboxylic_acid/long_alkyl_chain are "
     f"themselves correlated (phi≈0.80) and cannot be independently attributed; "
     f"tricyclic_core/tertiary_amine are not (phi≈0.21). Excluding moclobemide, "
     f"the morpholine association disappears entirely (moclobemide is 1 of only "
     f"2 morpholine carriers) while the others are essentially unchanged "
     f"(Supplementary Table S10) - morpholine's association is a single-drug "
     f"artifact. Feature roles are organ-specific, not global (Fig. 5).")
sub_head("D.  Toxicity Risk Score: Circular Fit vs. Honest Validation")
body(f"The fitted weights (Table IV; Fig. 3) placed the tricyclic core as the "
     f"strongest toxicophore (weight {top_tox['weight']:+.2f}). In-sample, the "
     f"structure-only TRS correlated with observed FAERS toxicity (Spearman ρ = "
     f"{RHO_IN:+.2f}, p = {PVAL_IN:.3f}, n = {len(trs)}; Fig. 1 left) - but this "
     f"is CIRCULAR, since the weights are derived from these same drugs. "
     f"Leave-one-drug-out cross-validation gives ρ = {RHO_LOO:+.2f} "
     f"(p = {PVAL_LOO:.2f}; Fig. 1 right): the score does NOT predict held-out "
     f"drugs. Excluding moclobemide left both conclusions unchanged. We "
     f"therefore present the TRS as a transparent summary of exploratory "
     f"associations, not a validated predictive score.")
sub_head("E.  Candidate Screening")
body(f"Of {N_SCREEN:,} ChEMBL compounds, {N_PASS} carried a safeophore and no "
     f"toxicophore ({N_DL} Lipinski-ok; Table V). Novelty was quantified: "
     f"{N_NOVEL}/{N_PASS} were structurally novel (Tanimoto <0.85 to any "
     f"reference drug) and {N_APPROVED_OVERLAP}/{N_PASS} were already an "
     f"approved drug elsewhere - genuinely novel candidates, offered as "
     f"hypothesis-generating leads conditioned on the exploratory, "
     f"non-generalising evidence above.")

# ---- IV. DISCUSSION ----
section_head("IV.  Discussion")
body(f"Integrating spontaneous-reporting pharmacovigilance with cheminformatics, "
     f"we found nominal, mechanistically plausible structure-toxicity "
     f"associations, but two checks often omitted elsewhere - FDR correction and "
     f"out-of-sample validation - substantially temper the claims. None of "
     f"{N_TESTS} tests survive FDR correction, and the TRS's apparent validation "
     f"(ρ~{RHO_IN:.2f} in-sample) collapses under leave-one-out cross-validation "
     f"(ρ~{RHO_LOO:.2f}), indicating {N_FEATURES} features and {N_DRUG}-"
     f"{N_CURATED} drugs do not carry enough independent information to "
     f"predict an unseen drug's toxicity from structure alone here.")
body("The candidate screen fared better under scrutiny: nearly all passing "
     "candidates were genuinely novel by Tanimoto similarity and none were "
     "already-approved drugs, so the screen is not simply rediscovering known "
     "agents - though that reassurance is still bounded by the exploratory "
     "evidence used to rank them. An isomeric error in two source structures "
     "(citalopram, escitalopram), found and corrected during curation before any "
     "analysis, initially would have produced a spurious cardiac safeophore; "
     "correcting it, the analysis instead converges on the expected "
     "tricyclic-core association - underscoring that structure curation is "
     "central to reproducible work, which is why we verified all 28 structures "
     "against PubChem, not only the two known-erroneous ones.")

# ---- V. LIMITATIONS ----
section_head("V.  Limitations")
body(f"The analysis rests on {N_DRUG} drugs with FAERS data; correlations do not "
     f"survive FDR correction, and the TRS does not generalise under "
     f"leave-one-drug-out cross-validation (ρ {RHO_IN:+.2f} in-sample vs. "
     f"{RHO_LOO:+.2f} out-of-sample) - both hypothesis-generating, not "
     f"confirmatory or clinically actionable. FAERS lacks denominators and is "
     f"subject to reporting biases and confounding by age, indication and "
     f"polypharmacy; disproportionality reflects over/under-reporting relative "
     f"to a chosen background, not incidence or causation. A within-class "
     f"background was used by design (Section II-C), yielding relative, not "
     f"absolute, risk. Per-drug counts were capped at 2,500, binding for "
     f"{N_CAPPED}/{N_DRUG} drugs (mostly high-volume agents). Candidate "
     f"molecules are computational leads requiring experimental follow-up. The "
     f"morpholine-neurologic association is a single-drug (moclobemide) "
     f"artifact and should not be generalised.")

# ---- VI. CONCLUSION ----
section_head("VI.  Conclusion")
body(f"Across a curated {N_DRUG}-drug antidepressant set, structural features "
     f"show exploratory, mechanistically plausible but FDR-non-significant "
     f"associations with real-world organ-system toxicity. A transparent "
     f"Toxicity Risk Score fits the data it was derived from "
     f"(ρ={RHO_IN:+.2f}) but does not predict held-out drugs "
     f"(LOOCV ρ={RHO_LOO:+.2f}), and should be treated as a hypothesis-"
     f"generating summary, not a validated prediction tool. Applied to public "
     f"chemical databases, the same rules do identify a large set of "
     f"genuinely novel, target-relevant candidates. A substantially larger, "
     f"more diverse drug set, with pre-registered multiplicity control and "
     f"cross-validation reported by default, is needed before this framework "
     f"could support structure-aware safety triage in antidepressant discovery.")

# ===================== full-width tables and figures =========================
newsec(1)
flagged = (sig[sig["signal"] == 1][["drug_name", "organ_system", "n_cases",
           "ror", "ror_lower", "ror_upper"]].round(2)
           .sort_values("ror", ascending=False))
flagged.columns = ["Drug", "Organ", "Cases", "ROR", "ROR low", "ROR high"]
ieee_table(flagged, "TABLE I",
           "PRIMARY DISPROPORTIONALITY SIGNALS, ROR LOWER CI > 2 (RANKED BY ROR)")

rk = (rank[["class", "safety_rank_in_class", "drug_name"] + ORGANS +
      ["n_signals", "toxicity_index"]].sort_values(["class", "safety_rank_in_class"]))
rk.columns = ["Class", "Rank", "Drug", "Card", "Hep", "Neuro", "Metab",
              "Sig", "Tox idx"]
ieee_table(rk, "TABLE II", "WITHIN-CLASS SAFETY RANKING (RANK 1 = SAFEST)")

cs = (corr[corr["p_value"] < 0.10][["organ_system", "feature", "rank_biserial",
      "p_value", "p_adj_bh", "role"]].sort_values("p_value").round(4))
cs.columns = ["Organ", "Feature", "Effect", "p", "q (FDR)", "Role"]
ieee_table(cs, "TABLE III",
           "STRUCTURE-TOXICITY ASSOCIATIONS, p<0.10 (NONE SURVIVE q<0.05)")

w2 = wts.round(3); w2.columns = ["Feature", "TRS weight"]
ieee_table(w2, "TABLE IV", "TOXICITY RISK SCORE WEIGHTS")

top = (cand[cand["pass"]].head(12)[["molecule_chembl_id", "via_target", "MW",
       "LogP", "max_tanimoto_to_known", "safeophores", "TRS"]])
top.columns = ["ChEMBL ID", "Target", "MW", "LogP", "Max Tanimoto", "Safeophores", "TRS"]
ieee_table(top, "TABLE V", "TOP PRIORITIZED CANDIDATE MOLECULES")

ieee_figure("charts/trs_validation.png",
            f"Fig. 1.  TRS in-sample fit (ρ={RHO_IN:+.2f}) vs. leave-one-drug-out "
            f"honest prediction (ρ={RHO_LOO:+.2f}), n={len(trs)}.", 4.6)
ieee_figure("charts/ror_heatmap.png",
            "Fig. 2.  Reporting Odds Ratio by drug and organ. Solid outline = "
            "primary signal (ROR lower CI>2); dashed = legacy ROR>1.", 3.6)
ieee_figure("charts/trs_weights.png",
            "Fig. 3.  Fitted structural-feature weights for the TRS. None "
            "survive FDR correction.", 4.4)
ieee_figure("charts/class_toxicity.png",
            "Fig. 4.  Per-drug toxicity index coloured by class.", 4.6)
ieee_figure("charts/feature_stability.png",
            "Fig. 5.  Structural-feature effect direction by organ system - "
            "sign flips show organ-specific, not global, roles.", 4.2)

# ===================== references (two-column) ==============================
newsec(2)
section_head("References")
refs = [
    "D. J. Brody and Q. Gu, “Antidepressant use among adults: United States, "
    "2015–2018,” NCHS Data Brief, no. 377, pp. 1–8, Sep. 2020.",
    "H. K. R. Thanacoody and S. H. L. Thomas, “Tricyclic antidepressant "
    "poisoning: cardiovascular toxicity,” Toxicol. Rev., vol. 24, no. 3, "
    "pp. 205–214, 2005.",
    "C. S. Voican, E. Corruble, S. Naveau, and G. Perlemuter, "
    "“Antidepressant-induced liver injury: a review for clinicians,” Am. J. "
    "Psychiatry, vol. 171, no. 4, pp. 404–415, 2014.",
    "E. W. Boyer and M. Shannon, “The serotonin syndrome,” N. Engl. J. Med., "
    "vol. 352, no. 11, pp. 1112–1120, 2005.",
    "K. J. Rothman, S. Lanes, and S. T. Sacks, “The reporting odds ratio and "
    "its advantages over the proportional reporting ratio,” Pharmacoepidemiol. "
    "Drug Saf., vol. 13, no. 8, pp. 519–523, 2004.",
    "S. J. W. Evans, P. C. Waller, and S. Davis, “Use of proportional reporting "
    "ratios (PRRs) for signal generation from spontaneous adverse drug reaction "
    "reports,” Pharmacoepidemiol. Drug Saf., vol. 10, no. 6, pp. 483–486, "
    "2001.",
    "A. Bate et al., “A Bayesian neural network method for adverse drug "
    "reaction signal generation,” Eur. J. Clin. Pharmacol., vol. 54, no. 4, "
    "pp. 315–321, 1998.",
    "E. P. van Puijenbroek et al., “A comparison of measures of "
    "disproportionality for signal detection in spontaneous reporting systems for "
    "adverse drug reactions,” Pharmacoepidemiol. Drug Saf., vol. 11, no. 1, "
    "pp. 3–10, 2002.",
    "E. G. Brown, L. Wood, and S. Wood, “The medical dictionary for regulatory "
    "activities (MedDRA),” Drug Saf., vol. 20, no. 2, pp. 109–117, 1999.",
    "G. Landrum, “RDKit: Open-source cheminformatics.” [Online]. "
    "Available: https://www.rdkit.org",
    "D. Mendez et al., “ChEMBL: towards direct deposition of bioassay "
    "data,” Nucleic Acids Res., vol. 47, no. D1, pp. D930–D940, 2019.",
    "C. A. Lipinski, F. Lombardo, B. W. Dominy, and P. J. Feeney, “Experimental "
    "and computational approaches to estimate solubility and permeability in drug "
    "discovery and development settings,” Adv. Drug Deliv. Rev., vol. 46, "
    "no. 1–3, pp. 3–26, 2001.",
    "U.S. Food and Drug Administration, “openFDA: FDA Adverse Event Reporting "
    "System (FAERS).” [Online]. Available: https://open.fda.gov",
    "H. Pajouhesh and G. R. Lenz, “Medicinal chemistry properties associated "
    "with blood-brain barrier permeability,” NeuroRx, vol. 2, no. 4, "
    "pp. 541–553, 2005.",
    "Y. Benjamini and Y. Hochberg, “Controlling the false discovery rate: a "
    "practical and powerful approach to multiple testing,” J. R. Stat. Soc. B, "
    "vol. 57, no. 1, pp. 289–300, 1995.",
    "D. F. Veber et al., “Molecular properties that influence the oral "
    "bioavailability of drug candidates,” J. Med. Chem., vol. 45, no. 12, "
    "pp. 2615–2623, 2002.",
    "A. K. Ghose, V. N. Viswanadhan, and J. J. Wendoloski, “A knowledge-based "
    "approach in designing combinatorial or medicinal chemistry libraries for "
    "drug discovery,” J. Comb. Chem., vol. 1, no. 1, pp. 55–68, 1999.",
]
for i, r in enumerate(refs, 1):
    p = doc.add_paragraph()
    p.paragraph_format.left_indent = Inches(0.17)
    p.paragraph_format.first_line_indent = Inches(-0.17)
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    run = p.add_run(f"[{i}]\t{r}"); run.font.size = Pt(8); run.font.name = FONT

doc.save("Research_Paper_IEEE.docx")
print("Wrote Research_Paper_IEEE.docx")
