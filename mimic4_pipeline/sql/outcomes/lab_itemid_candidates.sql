-- Candidate d_labitems rows for the 4 outcome analytes (ALT, AST, total
-- bilirubin, serum sodium), found by label pattern -- NOT a hardcoded
-- itemid list (brief, task 4: "do not hardcode itemid integers from
-- memory"). Deliberately over-inclusive across fluid/category so
-- specimen-type duplicates (e.g. blood chemistry vs whole-blood-gas
-- sodium) are visible for review rather than silently resolved; the
-- fluid/category restriction and the check for which candidates actually
-- have data happen in Python (resolve_lab_itemids() in
-- src/mimic4_pipeline/outcomes.py), which prints its reasoning.
--
-- Note the AST label is misspelled "Asparate Aminotransferase" in
-- MIMIC-IV's own d_labitems dictionary (in addition to the correctly
-- spelled "Aspartate..." variant, which does not actually occur) --
-- discovered by querying the real dictionary, not assumed.
select 'ALT' as analyte, itemid, label, fluid, category
from mimiciv_hosp.d_labitems
where lower(label) like '%alanine aminotransferase%'
union all
select 'AST' as analyte, itemid, label, fluid, category
from mimiciv_hosp.d_labitems
where lower(label) like '%aspartate aminotransferase%' or lower(label) like '%asparate aminotransferase%'
union all
select 'Total Bilirubin' as analyte, itemid, label, fluid, category
from mimiciv_hosp.d_labitems
where lower(label) like 'bilirubin, total%'
union all
select 'Sodium' as analyte, itemid, label, fluid, category
from mimiciv_hosp.d_labitems
where lower(label) like 'sodium%' and fluid = 'Blood'
order by analyte, fluid, category, itemid;
