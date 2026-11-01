-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
CREATE TEMP TABLE audit_population_ids(sid bigint);
\copy audit_population_ids FROM '@CUFF_WORKSPACE@/phase3/_restricted/population_audit_ids_eicu.csv' CSV HEADER
BEGIN TRANSACTION READ ONLY;
\copy (SELECT k.sid,p.unittype,p.apacheadmissiondx FROM audit_population_ids k LEFT JOIN public.patient p ON p.patientunitstayid=k.sid ORDER BY k.sid) TO '@CUFF_WORKSPACE@/phase3/_restricted/population_audit_raw_eicu.csv' CSV HEADER
ROLLBACK;
