-- Given a set of formulary_drug_cd / gsn / ndc values already confirmed (by
-- drug-text matching) to belong to a curated antidepressant, find every
-- OTHER prescriptions row that shares one of those codes -- i.e. the same
-- product ordered under a drug-text spelling the text matcher didn't
-- recognize. formulary_drug_cd/gsn/ndc are more stable identifiers than
-- free text (brief, section 2), so this is a recovery/cross-check pass,
-- not the primary matching method.
--
-- Parameters (substituted by src/mimic4_pipeline/build_drug_mapping.py):
--   {formulary_prefixes} -- SQL boolean expression over formulary_drug_cd
--   {gsn_list}            -- comma-separated, single-quoted gsn values
--   {matched_drug_values} -- comma-separated, single-quoted drug values already matched by text
select distinct drug, formulary_drug_cd, gsn, ndc
from mimiciv_hosp.prescriptions
where ({formulary_prefixes} or gsn in ({gsn_list}))
  and drug not in ({matched_drug_values})
order by drug;
