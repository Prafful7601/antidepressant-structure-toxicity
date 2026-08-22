-- Task 5 investigation: search icu.d_items for anything QT/QTc/conduction-
-- interval related, rather than assuming where (or whether) it exists.
select itemid, label, abbreviation, category, unitname, param_type
from mimiciv_icu.d_items
where lower(label) like '%qt%' or lower(label) like '%qrs%' or lower(label) like '%pr interval%'
   or lower(label) like '%conduction%' or lower(abbreviation) like '%qt%'
order by category, label;
