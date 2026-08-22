-- Task 5 investigation: how populated is the one QTc item found in
-- icu.chartevents (itemid 224359, "QTc", category "Routine Vital Signs")?
-- Demo result (2026-08-23): 14 rows total, 4 distinct patients, 4 distinct
-- admissions -- out of 668,862 chartevents rows and 140 icustays in the
-- whole demo. Too sparse to be a usable endpoint on its own; see the
-- README's task 5 section for the full writeup and recommendation.
select
    count(*)                          as n_rows,
    count(distinct subject_id)        as n_patients,
    count(distinct hadm_id)           as n_admissions,
    min(charttime)                    as min_charttime,
    max(charttime)                    as max_charttime
from mimiciv_icu.chartevents
where itemid = 224359;
