-- Pre-existing liver disease covariate, from mimic-code's own Charlson
-- comorbidity concept (mimic-iv/concepts_duckdb/comorbidity/charlson.sql),
-- already materialized into mimiciv_derived.charlson by build_mimic.sh.
-- Reused rather than re-derived from ICD codes by hand, per the brief's
-- instruction (task 6) to use mimic-code concept scripts where available --
-- pulled forward here since task 3 needs it as a covariate too.
select
    hadm_id,
    mild_liver_disease,
    severe_liver_disease,
    (mild_liver_disease = 1 or severe_liver_disease = 1) as has_liver_disease,
    charlson_comorbidity_index
from mimiciv_derived.charlson;
