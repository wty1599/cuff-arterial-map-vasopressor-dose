-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
CREATE TEMP TABLE audit_population_ids(sid bigint);
\copy audit_population_ids FROM '@CUFF_WORKSPACE@/phase3/_restricted/population_audit_ids_mimic.csv' CSV HEADER
BEGIN TRANSACTION READ ONLY;
\copy (SELECT k.sid,i.first_careunit,s.stay_id AS sepsis_stay_id,s.sepsis3,s.suspected_infection_time,s.sofa_time,s.sofa_score FROM audit_population_ids k LEFT JOIN mimiciv_icu.icustays i ON i.stay_id=k.sid LEFT JOIN mimiciv_derived.sepsis3 s ON s.stay_id=k.sid ORDER BY k.sid,s.suspected_infection_time) TO '@CUFF_WORKSPACE@/phase3/_restricted/population_audit_raw_mimic.csv' CSV HEADER
ROLLBACK;
