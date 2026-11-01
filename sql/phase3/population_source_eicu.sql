-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
CREATE TEMP TABLE population_ids(sid bigint);
\copy population_ids FROM '@CUFF_WORKSPACE@/phase3/_restricted/population_ids_eicu.csv' CSV HEADER
BEGIN TRANSACTION READ ONLY;
\copy (SELECT p.patientunitstayid AS sid,p.unittype,p.apacheadmissiondx FROM population_ids k JOIN public.patient p ON p.patientunitstayid=k.sid ORDER BY p.patientunitstayid) TO '@CUFF_WORKSPACE@/phase3/_restricted/population_source_eicu.csv' CSV HEADER
\copy (SELECT unittype AS careunit,count(*) AS stays FROM public.patient GROUP BY unittype ORDER BY unittype) TO '@CUFF_WORKSPACE@/phase3/population_dictionary_eicu_units.csv' CSV HEADER
\copy (SELECT apacheadmissiondx AS diagnosis,count(*) AS stays FROM public.patient WHERE apacheadmissiondx ILIKE '%sepsis%' GROUP BY apacheadmissiondx ORDER BY apacheadmissiondx) TO '@CUFF_WORKSPACE@/phase3/population_dictionary_eicu_sepsis.csv' CSV HEADER
ROLLBACK;
