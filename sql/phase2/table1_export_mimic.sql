-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
CREATE TEMP TABLE table1_ids (sid bigint);
\copy table1_ids FROM '@CUFF_WORKSPACE@/phase2/_restricted/table1_ids_mimic.csv' CSV HEADER
BEGIN TRANSACTION READ ONLY;
\copy (SELECT i.stay_id AS sid,i.subject_id,i.hadm_id,i.first_careunit,a.admission_type,a.admission_location,a.hospital_expire_flag,p.gender,s.sofa,s.respiration,s.coagulation,s.liver,s.cardiovascular,s.cns,s.renal FROM table1_ids k JOIN mimiciv_icu.icustays i ON i.stay_id=k.sid JOIN mimiciv_hosp.admissions a ON a.hadm_id=i.hadm_id JOIN mimiciv_hosp.patients p ON p.subject_id=i.subject_id LEFT JOIN mimiciv_derived.first_day_sofa s ON s.stay_id=i.stay_id) TO '@CUFF_WORKSPACE@/phase2/_restricted/table1_mimic_source.csv' CSV HEADER
ROLLBACK;
