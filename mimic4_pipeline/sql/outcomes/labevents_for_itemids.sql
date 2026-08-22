-- Pulls labevents rows for a resolved set of itemids, joined to admittime
-- for baseline/peak timing. `valuenum` is MIMIC's own numeric-parsed
-- version of the raw (sometimes non-numeric, e.g. ">1000" or "ERROR")
-- `value` free-text field; filtering on valuenum is not null is exactly
-- the brief's "handle the known issue that value can contain non-numeric
-- strings" instruction -- no extra parsing needed on our end.
--
-- Parameters (substituted by src/mimic4_pipeline/outcomes.py):
--   {itemid_list} -- comma-separated itemid integers
-- ref_range_upper/lower are carried through so the trajectory shape can
-- express the extreme value as a multiple of its own upper limit of normal
-- (ULN) -- the standard Hy's-Law-style way of expressing a clinically
-- meaningful hepatic signal (>=3x ULN for ALT/AST, >=2x ULN for
-- bilirubin), rather than an arbitrary raw-unit cutoff. Confirmed against
-- the demo: 100% populated and consistent for all 4 resolved itemids
-- (sodium has two slightly different ranges, likely an assay/instrument
-- change over the demo's date range -- used per-draw, not a fixed
-- constant, for exactly that reason).
select
    le.subject_id,
    le.hadm_id,
    le.itemid,
    a.admittime,
    le.charttime,
    le.valuenum,
    le.ref_range_upper,
    le.ref_range_lower
from mimiciv_hosp.labevents le
join mimiciv_hosp.admissions a using (hadm_id)
where le.itemid in ({itemid_list})
  and le.valuenum is not null
  and le.hadm_id is not null;
