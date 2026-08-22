-- All distinct free-text drug values in hosp.prescriptions, with their
-- associated formulary_drug_cd/gsn/ndc combinations and row/admission
-- counts. This is the raw material the drug matcher (see
-- src/mimic4_pipeline/drug_mapping.py) is run against -- inspect actual
-- data, never assume spellings.
select
    drug,
    formulary_drug_cd,
    gsn,
    ndc,
    route,
    form_rx,
    count(*)               as n_prescriptions,
    count(distinct hadm_id) as n_admissions
from mimiciv_hosp.prescriptions
group by all
order by drug, formulary_drug_cd;
