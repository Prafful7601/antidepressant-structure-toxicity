-- Resolves the actual ICD codes for antidepressant poisoning/overdose from
-- the dictionary at runtime -- never hardcode these from memory (same
-- principle as lab itemid resolution in task 4).
--
-- ICD-10-CM: T43.0 (tricyclic/tetracyclic), T43.1 (MAOI), and T43.2
-- (other/unspecified antidepressants, which is where SSRI/SNRI live --
-- T43.21 = SNRI, T43.22 = SSRI) subcodes. Deliberately narrower than the
-- full "T43.x range" the brief mentions: T43.3-T43.9 cover
-- antipsychotic/psychostimulant/other-psychotropic poisoning, a different
-- confound than the antidepressant-overdose concern this flag targets.
-- Also excludes T43's "Adverse effect of..." and "Underdosing of..."
-- subcodes, which share the T43.0/.1/.2 numeric prefix but are not
-- overdoses at all -- distinguished only by long_title text, not the code
-- prefix, which is exactly why this is resolved from the dictionary
-- instead of a hardcoded regex on the code alone. (Verified against the
-- demo: one admission carries T43505S, "Adverse effect of unspecified
-- antipsychotics and neuroleptics, sequela" -- caught by a naive T43% LIKE
-- filter, correctly excluded by this one.)
--
-- ICD-9-CM: 969.0x ("Poisoning by antidepressant" and its per-class
-- subcodes 96900-96909) -- ICD-9 has no separate adverse-effect/underdosing
-- split in this range, so no long_title filter is needed there.
select icd_code, icd_version, long_title
from mimiciv_hosp.d_icd_diagnoses
where
    (
        icd_version = 10
        and (icd_code like 'T430%' or icd_code like 'T431%' or icd_code like 'T432%')
        and long_title like 'Poisoning by%'
    )
    or (icd_version = 9 and icd_code like '9690%')
order by icd_version, icd_code;
