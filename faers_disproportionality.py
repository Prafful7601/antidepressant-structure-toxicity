"""
FAERS disproportionality analysis for the antidepressant toxicity study.
Covers Objectives 1-2: quantify cardiac/hepatic/neurologic/metabolic signals
across the 26 antidepressants and rank drugs by safety within class.

This script assumes you have already:
  1. Downloaded FAERS quarterly ASCII files from
     https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html
     (DEMO, DRUG, REAC tables for each quarter)
  2. Loaded and concatenated them into a single tidy table `faers_reports.csv`
     with at minimum these columns:
        primaryid, drugname, pt   (pt = MedDRA Preferred Term of the reaction)
     One row per (report x drug x reaction).

If you do NOT yet have FAERS loaded, see load_faers_ascii() below for a loader,
and the MedDRA grouping in ORGAN_SYSTEM_PT for how reactions map to the four
organ systems. The PT lists here are STARTERS — expand them with the full
MedDRA SMQ / HLGT terms for a publishable analysis.

Disproportionality measures computed per (drug, organ_system):
  - ROR (Reporting Odds Ratio) + 95% CI
  - PRR (Proportional Reporting Ratio) + chi-square
  - IC  (Information Component, Bayesian) + 95% credible interval (IC025, IC975)
Primary signal criterion: ROR lower 95% bound > 2 AND n_cases >= 3 (standard
pharmacovigilance practice; stricter than the permissive ROR>1 screen). The
legacy ROR>1 flag is retained as `signal_ror1` for supplementary comparison
only - it is NOT used for any headline result.

Usage:
    python faers_disproportionality.py
"""

import pandas as pd
import numpy as np

# ----------------------------------------------------------------------
# Drug name normalization: map FAERS free-text drugnames to our 26 drugs.
# FAERS has brand names, salts, misspellings - extend these lists.
# ----------------------------------------------------------------------
DRUG_SYNONYMS = {
    "Moclobemide": ["moclobemide", "manerix", "aurorix"],
    "Clorgyline": ["clorgyline", "clorgiline"],
    "Fluoxetine": ["fluoxetine", "prozac", "sarafem", "fluoxetine hcl"],
    "Fluvoxamine": ["fluvoxamine", "luvox", "faverin"],
    "Paroxetine": ["paroxetine", "paxil", "seroxat", "paroxetine hcl"],
    "Sertraline": ["sertraline", "zoloft", "lustral", "sertraline hcl"],
    "Citalopram": ["citalopram", "celexa", "cipramil"],
    "Escitalopram": ["escitalopram", "lexapro", "cipralex"],
    "Dapoxetine": ["dapoxetine", "priligy"],
    "Trazodone": ["trazodone", "desyrel", "oleptro"],
    "Mianserin": ["mianserin", "tolvon"],
    "Mirtazapine": ["mirtazapine", "remeron", "zispin"],
    "Bupropion": ["bupropion", "wellbutrin", "zyban", "bupropion hcl"],
    "Tianeptine": ["tianeptine", "stablon", "coaxil", "tatinol"],
    "Amineptine": ["amineptine", "survector"],
    "Atomoxetine": ["atomoxetine", "strattera"],
    "Imipramine": ["imipramine", "tofranil"],
    "Amitriptyline": ["amitriptyline", "elavil", "endep"],
    "Trimipramine": ["trimipramine", "surmontil"],
    "Doxepin": ["doxepin", "sinequan", "silenor"],
    "Dothiepin": ["dothiepin", "dosulepin", "prothiaden"],
    "Clomipramine": ["clomipramine", "anafranil"],
    "Desipramine": ["desipramine", "norpramin"],
    "Nortriptyline": ["nortriptyline", "pamelor", "aventyl"],
    "Amoxapine": ["amoxapine", "asendin"],
    "Reboxetine": ["reboxetine", "edronax"],
    "Venlafaxine": ["venlafaxine", "effexor"],
    "Duloxetine": ["duloxetine", "cymbalta"],
}
SYN_TO_DRUG = {syn: drug for drug, syns in DRUG_SYNONYMS.items() for syn in syns}


# ----------------------------------------------------------------------
# MedDRA PT -> organ system mapping (STARTER lists; expand for publication).
# For a real study, use the MedDRA SMQ (Standardised MedDRA Queries) or
# group by HLGT/SOC. These illustrative PTs let the pipeline run today.
# ----------------------------------------------------------------------
ORGAN_SYSTEM_PT = {
    "cardiac": [
        # --- repolarisation / QT (SMQ: Torsade de pointes/QT prolongation) ---
        "qt prolongation", "electrocardiogram qt prolonged", "long qt syndrome",
        "torsade de pointes", "ventricular tachycardia", "ventricular arrhythmia",
        "ventricular fibrillation", "ventricular extrasystoles",
        "ventricular flutter", "cardiac fibrillation",
        # --- rate / rhythm ---
        "tachycardia", "sinus tachycardia", "supraventricular tachycardia",
        "bradycardia", "sinus bradycardia", "arrhythmia", "palpitations",
        "atrial fibrillation", "atrial flutter", "extrasystoles",
        "heart rate increased", "heart rate irregular",
        "atrioventricular block", "bundle branch block",
        # --- ischaemia / pump failure / arrest ---
        "myocardial infarction", "acute myocardial infarction",
        "myocardial ischaemia", "angina pectoris", "cardiac arrest",
        "cardio-respiratory arrest", "cardiac failure",
        "cardiac failure congestive", "cardiac failure acute",
        "cardiomyopathy", "myocarditis", "pericarditis",
        "ejection fraction decreased", "electrocardiogram st segment depression",
        "electrocardiogram abnormal", "sudden cardiac death", "sudden death",
    ],
    "hepatic": [
        # --- injury / failure (SMQ: Drug-related hepatic disorders) ---
        "hepatotoxicity", "drug-induced liver injury", "hepatitis",
        "hepatitis acute", "hepatitis fulminant", "hepatitis cholestatic",
        "hepatitis toxic", "autoimmune hepatitis", "hepatocellular injury",
        "liver injury", "hepatic failure", "acute hepatic failure",
        "hepatic necrosis", "hepatic cirrhosis", "hepatic steatosis",
        "hepatic function abnormal", "hepatic disorder", "hepatic enzyme increased",
        # --- cholestasis / jaundice ---
        "cholestasis", "jaundice", "cholestatic jaundice", "jaundice cholestatic",
        "hyperbilirubinaemia", "blood bilirubin increased",
        # --- lab abnormalities ---
        "alanine aminotransferase increased", "aspartate aminotransferase increased",
        "gamma-glutamyltransferase increased", "blood alkaline phosphatase increased",
        "liver function test abnormal", "liver function test increased",
        "transaminases increased", "hepatic enzyme abnormal", "ascites",
    ],
    "neurologic": [
        # --- serotonergic / NMS ---
        "serotonin syndrome", "neuroleptic malignant syndrome",
        # --- seizures ---
        "seizure", "convulsion", "generalised tonic-clonic seizure",
        "status epilepticus", "partial seizures", "grand mal convulsion",
        # --- movement / EPS ---
        "tremor", "extrapyramidal disorder", "dystonia", "akathisia",
        "dyskinesia", "tardive dyskinesia", "myoclonus", "bradykinesia",
        "parkinsonism", "restless legs syndrome", "muscle rigidity",
        "movement disorder",
        # --- consciousness / cognition / other CNS ---
        "encephalopathy", "coma", "loss of consciousness", "syncope",
        "dizziness", "headache", "paraesthesia", "ataxia", "confusional state",
        "disturbance in attention", "memory impairment", "cognitive disorder",
        "somnolence", "sedation", "hypoaesthesia", "speech disorder",
        "cerebrovascular accident", "balance disorder",
    ],
    "metabolic": [
        # --- weight / glycaemic ---
        "weight increased", "weight decreased", "obesity", "increased appetite",
        "decreased appetite", "hyperglycaemia", "blood glucose increased",
        "diabetes mellitus", "type 2 diabetes mellitus", "hypoglycaemia",
        "blood glucose decreased", "glucose tolerance impaired",
        "diabetic ketoacidosis", "glycosylated haemoglobin increased",
        # --- electrolytes / SIADH ---
        "hyponatraemia", "syndrome of inappropriate antidiuretic hormone secretion",
        "blood sodium decreased", "hypokalaemia", "hyperkalaemia",
        # --- lipids / metabolic syndrome ---
        "hypercholesterolaemia", "hypertriglyceridaemia", "dyslipidaemia",
        "blood cholesterol increased", "blood triglycerides increased",
        "metabolic syndrome", "hyperlipidaemia",
        # --- prolactin / endocrine commonly linked ---
        "hyperprolactinaemia", "blood prolactin increased",
    ],
}
PT_TO_ORGAN = {pt: org for org, pts in ORGAN_SYSTEM_PT.items() for pt in pts}


# ----------------------------------------------------------------------
# Optional loader for raw FAERS ASCII quarterly files.
# ----------------------------------------------------------------------
def load_faers_ascii(quarter_dirs):
    """
    quarter_dirs: list of paths, each containing DEMOyyQq.txt, DRUGyyQq.txt,
    REACyyQq.txt ('$'-delimited). Returns tidy DataFrame: primaryid, drugname, pt.
    """
    import glob, os
    frames = []
    for d in quarter_dirs:
        drug_f = glob.glob(os.path.join(d, "DRUG*.txt"))[0]
        reac_f = glob.glob(os.path.join(d, "REAC*.txt"))[0]
        drug = pd.read_csv(drug_f, sep="$", dtype=str, on_bad_lines="skip")
        reac = pd.read_csv(reac_f, sep="$", dtype=str, on_bad_lines="skip")
        drug.columns = [c.lower() for c in drug.columns]
        reac.columns = [c.lower() for c in reac.columns]
        m = drug[["primaryid", "drugname"]].merge(
            reac[["primaryid", "pt"]], on="primaryid")
        frames.append(m)
    out = pd.concat(frames, ignore_index=True)
    out["drugname"] = out["drugname"].str.lower().str.strip()
    out["pt"] = out["pt"].str.lower().str.strip()
    return out


# ----------------------------------------------------------------------
# Disproportionality statistics
# ----------------------------------------------------------------------
def disproportionality(contingency):
    """
    contingency: DataFrame with columns a,b,c,d per (drug, organ_system) where
       a = reports with drug AND event
       b = reports with drug, NOT event
       c = reports with event, NOT drug
       d = reports with neither
    Returns the same frame with ROR, PRR, IC and bounds appended.
    """
    df = contingency.copy()
    a, b, c, d = df["a"], df["b"], df["c"], df["d"]

    # ROR
    df["ROR"] = (a / b) / (c / d)
    se_lnror = np.sqrt(1/a + 1/b + 1/c + 1/d)
    df["ROR025"] = np.exp(np.log(df["ROR"]) - 1.96 * se_lnror)
    df["ROR975"] = np.exp(np.log(df["ROR"]) + 1.96 * se_lnror)

    # PRR + chi-square (Yates)
    df["PRR"] = (a / (a + b)) / (c / (c + d))
    n = a + b + c + d
    expected_a = (a + b) * (a + c) / n
    df["chi2"] = (np.abs(a - expected_a) - 0.5) ** 2 / expected_a  # simplified

    # IC (Bayesian Confidence Propagation Neural Network style) + 95% credible interval
    e_a = (a + b) * (a + c) / n
    df["IC"] = np.log2((a + 0.5) / (e_a + 0.5))
    ic_halfwidth = 3.3 * (a + 0.5) ** (-0.5)  # approx normal-on-log2 95% half-width
    df["IC025"] = df["IC"] - ic_halfwidth   # lower bound
    df["IC975"] = df["IC"] + ic_halfwidth   # upper bound

    return df


def build_contingency(reports):
    """
    reports: tidy DataFrame primaryid, drug (mapped), organ_system.
    De-duplicate to one (report, drug, organ) then build 2x2 per pair.
    """
    total_reports = reports["primaryid"].nunique()
    rows = []
    drugs = reports["drug"].dropna().unique()
    organs = list(ORGAN_SYSTEM_PT.keys())
    for drug in drugs:
        for organ in organs:
            a = reports[(reports["drug"] == drug) &
                        (reports["organ_system"] == organ)]["primaryid"].nunique()
            drug_total = reports[reports["drug"] == drug]["primaryid"].nunique()
            organ_total = reports[reports["organ_system"] == organ]["primaryid"].nunique()
            b = drug_total - a
            c = organ_total - a
            d = total_reports - a - b - c
            rows.append({"drug_name": drug, "organ_system": organ,
                         "a": a, "b": b, "c": c, "d": d, "n_cases": a})
    return pd.DataFrame(rows)


ROR_THRESH = 2.0  # primary signal threshold on the ROR lower 95% bound


def run(reports_csv="faers_reports.csv", min_cases=3, ror_thresh=ROR_THRESH):
    reports = pd.read_csv(reports_csv, dtype=str)
    reports["drug"] = reports["drugname"].map(SYN_TO_DRUG)
    reports["organ_system"] = reports["pt"].map(PT_TO_ORGAN)
    reports = reports.dropna(subset=["drug", "organ_system"])

    cont = build_contingency(reports)
    # guard against zero cells (continuity correction)
    for col in ["a", "b", "c", "d"]:
        cont[col] = cont[col].replace(0, 0.5)
    stats = disproportionality(cont)

    # Primary signal: ROR lower 95% bound > 2 and >= 3 cases (standard practice).
    stats["signal"] = ((stats["ROR025"] > ror_thresh) &
                       (stats["n_cases"] >= min_cases)).astype(int)
    # Legacy ROR>1 flag, kept only for supplementary/sensitivity comparison.
    stats["signal_ror1"] = ((stats["ROR025"] > 1) &
                            (stats["n_cases"] >= min_cases)).astype(int)

    out = stats[["drug_name", "organ_system", "n_cases", "ROR", "ROR025",
                 "ROR975", "PRR", "chi2", "IC", "IC025", "IC975",
                 "signal", "signal_ror1"]]
    out = out.rename(columns={"ROR025": "ror_lower", "ROR975": "ror_upper",
                              "ROR": "ror"})
    out.to_csv("faers_signals.csv", index=False)
    print(f"Wrote faers_signals.csv: {out['signal'].sum()} signals at ROR>{ror_thresh:g} "
          f"(vs {out['signal_ror1'].sum()} at legacy ROR>1) "
          f"across {out['drug_name'].nunique()} drugs")
    return out


if __name__ == "__main__":
    import os
    if os.path.exists("faers_reports.csv"):
        print(run().to_string())
    else:
        print("No faers_reports.csv found.")
        print("Download FAERS quarterly ASCII files, then either:")
        print("  - load with load_faers_ascii([list_of_quarter_dirs]) and save "
              "as faers_reports.csv, or")
        print("  - build faers_reports.csv with columns: primaryid, drugname, pt")
        print("\nThen re-run this script.")
