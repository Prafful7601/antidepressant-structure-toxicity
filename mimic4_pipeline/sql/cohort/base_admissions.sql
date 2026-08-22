-- Starting population for the CONSORT flow: every hospital admission.
select
    subject_id,
    hadm_id,
    admittime,
    dischtime,
    deathtime,
    admission_type,
    hospital_expire_flag
from mimiciv_hosp.admissions;
