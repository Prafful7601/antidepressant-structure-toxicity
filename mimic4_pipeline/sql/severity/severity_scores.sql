-- One row per hospital admission with every ICU-day severity score
-- mimic-code provides, joined via that admission's *first* ICU stay (9 of
-- the demo's 275 admissions have more than one ICU stay -- the standard
-- choice is the first, since it reflects acuity closest to admission).
-- Admissions with no ICU stay at all (a substantial fraction in a
-- hospital-wide cohort -- 147 of 275 in the demo) get NULLs here, not a
-- fabricated score; Charlson (admission-level, not ICU-dependent) is the
-- covariate that covers the whole cohort and is joined separately in
-- severity.py.
--
-- Concept tables used, all materialized by mimic-code's own
-- mimic-iv/concepts_duckdb/duckdb.sql during scripts/build_database.sh
-- (task 1) -- ported, not reimplemented, per the brief's task-6
-- instruction:
--   firstday/first_day_sofa.sql  -> mimiciv_derived.first_day_sofa
--   score/sapsii.sql             -> mimiciv_derived.sapsii
--   score/apsiii.sql             -> mimiciv_derived.apsiii
--   score/sirs.sql               -> mimiciv_derived.sirs
--   score/oasis.sql              -> mimiciv_derived.oasis
--   score/lods.sql                -> mimiciv_derived.lods
with first_icu_stay as (
    select hadm_id, stay_id, intime
    from (
        select
            hadm_id, stay_id, intime,
            row_number() over (partition by hadm_id order by intime asc) as rn
        from mimiciv_icu.icustays
    )
    where rn = 1
)
select
    f.hadm_id,
    f.stay_id,
    fs.sofa           as first_day_sofa,
    sp.sapsii,
    sp.sapsii_prob,
    ap.apsiii,
    ap.apsiii_prob,
    si.sirs,
    oa.oasis,
    oa.oasis_prob,
    lo.lods
from first_icu_stay f
left join mimiciv_derived.first_day_sofa fs using (stay_id)
left join mimiciv_derived.sapsii sp using (stay_id)
left join mimiciv_derived.apsiii ap using (stay_id)
left join mimiciv_derived.sirs si using (stay_id)
left join mimiciv_derived.oasis oa using (stay_id)
left join mimiciv_derived.lods lo using (stay_id);
