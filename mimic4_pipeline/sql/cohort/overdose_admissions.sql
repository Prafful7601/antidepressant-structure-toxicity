-- hadm_ids carrying an antidepressant-poisoning diagnosis, per the code set
-- resolved in poisoning_diagnosis_codes.sql (kept in sync with it -- run
-- that file separately to review the resolved code list; both are driven
-- by the same runtime lookup against d_icd_diagnoses, not a hardcoded list).
select distinct d.hadm_id
from mimiciv_hosp.diagnoses_icd d
where
    (
        d.icd_version = 10
        and (d.icd_code like 'T430%' or d.icd_code like 'T431%' or d.icd_code like 'T432%')
        and d.icd_code in (
            select icd_code from mimiciv_hosp.d_icd_diagnoses
            where icd_version = 10 and long_title like 'Poisoning by%'
        )
    )
    or (d.icd_version = 9 and d.icd_code like '9690%');
