-- Full prescription-order rows for the drug-text values the matcher
-- accepted. Parameters (substituted by
-- src/mimic4_pipeline/build_drug_mapping.py):
--   {matched_drug_values} -- comma-separated, single-quoted drug values
select
    subject_id,
    hadm_id,
    pharmacy_id,
    drug,
    drug_type,
    formulary_drug_cd,
    gsn,
    ndc,
    route,
    form_rx,
    dose_val_rx,
    dose_unit_rx,
    doses_per_24_hrs,
    starttime,
    stoptime
from mimiciv_hosp.prescriptions
where drug in ({matched_drug_values})
order by subject_id, hadm_id, starttime;
