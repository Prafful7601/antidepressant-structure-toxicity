-- Propensity-model covariates for the Tricyclic vs SSRI+SNRI matched
-- comparison: admission-level age from mimic-code's own age concept (not
-- hand-computed from anchor_year -- avoids re-deriving something mimic-code
-- already gets right), sex, an emergency-vs-elective admission flag, and
-- every individual Charlson comorbidity component (not just the summary
-- index used elsewhere in this pipeline -- a richer covariate set is
-- standard practice for a propensity model, and mimic-code's own concept
-- already computes each component, so there's no reason to collapse to the
-- index here).
select
    a.hadm_id,
    ag.age,
    (p.gender = 'M')::int as is_male,
    (a.admission_type in ('EW EMER.', 'URGENT', 'DIRECT EMER.', 'AMBULATORY OBSERVATION'))::int as is_emergency_admission,
    c.myocardial_infarct,
    c.congestive_heart_failure,
    c.peripheral_vascular_disease,
    c.cerebrovascular_disease,
    c.dementia,
    c.chronic_pulmonary_disease,
    c.rheumatic_disease,
    c.peptic_ulcer_disease,
    c.mild_liver_disease,
    c.diabetes_without_cc,
    c.diabetes_with_cc,
    c.paraplegia,
    c.renal_disease,
    c.malignant_cancer,
    c.severe_liver_disease,
    c.metastatic_solid_tumor,
    c.aids
from mimiciv_hosp.admissions a
join mimiciv_hosp.patients p using (subject_id)
left join mimiciv_derived.age ag using (hadm_id)
left join mimiciv_derived.charlson c using (hadm_id);
