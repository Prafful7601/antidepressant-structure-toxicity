# Structure-Toxicity Analysis of Antidepressants

Repository: https://github.com/Prafful7601/antidepressant-structure-toxicity

Author: Prafful Gupta (first author)

Pharmacovigilance (FAERS) + cheminformatics (RDKit/ChEMBL/PubChem) pipeline
linking antidepressant chemical structure to organ-system toxicity signals,
with a structure-based Toxicity Risk Score (TRS) validated by leave-one-drug-out
cross-validation.

This repository holds the reproducible analysis pipeline (code + data +
figures) behind the manuscript. The manuscript itself (drafts, results/summary
write-ups, and their `.docx` generator scripts) is kept out of this public
repo to avoid pre-disclosing unpublished results ahead of submission; it is
maintained locally alongside this repo and regenerated from the same CSVs.

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Pipeline (run in this order)

1. `fetch_faers_openfda.py` - pulls FAERS adverse-event reports from the
   openFDA API for all 28 drugs -> `faers_reports.csv`. (Optional:
   `topup_faers.py` to re-fetch any drug truncated by a transient network
   error.) **Slow / network-bound; `faers_reports.csv` is already included so
   you don't have to re-run this to reproduce downstream results.**
2. `faers_disproportionality.py` - ROR/PRR/IC disproportionality per
   drug x organ system -> `faers_signals.csv`.
3. `rank_within_class.py` - within-pharmacological-class safety ranking ->
   `class_safety_ranking.csv`.
4. `cheminformatics_analysis.py` - RDKit descriptors/substructure flags,
   structure-toxicity correlation (Fisher + Mann-Whitney, both with
   Benjamini-Hochberg FDR correction), and structural-feature
   multicollinearity check -> `drug_features.csv`,
   `structure_toxicity_correlation*.csv`,
   `structural_feature_multicollinearity*.csv`.
5. `fit_trs.py` - fits the Toxicity Risk Score and validates it both
   in-sample (circular) and by leave-one-drug-out cross-validation (honest)
   -> `trs_weights.csv`, `drug_trs.csv`, `trs_loocv.csv`.
6. `chembl_candidate_screen.py` - screens ChEMBL for novel candidate
   molecules (Tanimoto novelty, Lipinski/Veber/Ghose, ChEMBL approval-status
   overlap check) -> `chembl_candidates_screened.csv`.
7. `sensitivity_moclobemide.py` - sensitivity analysis re-running the
   structure-toxicity correlation and TRS fit/LOOCV with moclobemide
   excluded -> `sensitivity_excl_moclobemide.csv`.
8. `pubchem_verify.py` - independently verifies all 28 curated SMILES against
   PubChem by InChIKey skeleton -> `pubchem_verification.csv`.
9. `faers_date_range.py` - determines the actual FAERS receive-date range
   captured and dumps the exact openFDA query string per drug ->
   `faers_date_range.csv`, `faers_query_strings.txt`.
10. `sanity_check.py` - independently re-derives key statistics directly from
    raw data and flags any inconsistency with the pipeline's own outputs.
11. `make_charts.py` - generates all figures (`charts/*.png`) from the CSVs
    above; this is the single source of truth for every chart, including
    those used in the (locally-maintained, not published here) manuscript.

Every script reads its inputs from CSVs in this directory and is safe to
re-run independently once its upstream CSVs exist. The manuscript's own
`.docx` generator scripts are not part of this repo (see above) but consume
exactly the CSVs and charts produced by the steps above.

## Key methodological notes

- **Primary disproportionality signal**: ROR lower 95% bound > 2 and >= 3
  cases (standard pharmacovigilance practice). The more permissive legacy
  ROR>1 criterion is retained only as a supplementary comparison
  (`signal_ror1` column in `faers_signals.csv`) and is not used for any
  headline result.
- **Multiple testing**: all structure-toxicity tests (4 organ systems x 16
  structural features) are corrected with Benjamini-Hochberg FDR; see the
  `p_adj_bh` / `significant_fdr05` columns.
- **TRS validation**: `drug_trs.csv` and `trs_loocv.csv` report both the
  in-sample fit (circular - weights are derived from the same drugs) and the
  leave-one-drug-out cross-validated prediction (honest). Report both
  together; never the in-sample number alone.
- **Clorgyline** is in the curated 28-drug structure set but was never
  marketed and returned zero FAERS reports, so it is excluded from all
  FAERS-dependent analyses (27 drugs with data).

## Data availability

FAERS data are retrieved from the public [openFDA API](https://open.fda.gov)
with no API key required; the exact query string used for every drug is in
`faers_query_strings.txt`, and the effective date range retrieved is in
`faers_date_range.csv`. ChEMBL and PubChem data are retrieved from their
respective public REST APIs. All intermediate and final CSVs referenced above
are included in this repository so the full analysis can be reproduced
without re-querying any external API.

## License

Code is released under the MIT License (see `LICENSE`). Data derived from
FAERS/openFDA, ChEMBL, and PubChem remain subject to those sources' own
terms of use.
