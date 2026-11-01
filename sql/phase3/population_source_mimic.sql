-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
CREATE TEMP TABLE population_ids(sid bigint);
\copy population_ids FROM '@CUFF_WORKSPACE@/phase3/_restricted/population_ids_mimic.csv' CSV HEADER
BEGIN TRANSACTION READ ONLY;
\copy (SELECT i.stay_id AS sid,i.first_careunit,count(s.stay_id) sepsis_rows,count(*) FILTER(WHERE s.sepsis3 IS TRUE) sepsis_positive_rows,count(*) FILTER(WHERE s.sepsis3 IS FALSE) sepsis_negative_rows,count(*) FILTER(WHERE s.stay_id IS NOT NULL AND s.sepsis3 IS NULL) sepsis_null_rows FROM population_ids k JOIN mimiciv_icu.icustays i ON i.stay_id=k.sid LEFT JOIN mimiciv_derived.sepsis3 s ON s.stay_id=i.stay_id GROUP BY i.stay_id,i.first_careunit ORDER BY i.stay_id) TO '@CUFF_WORKSPACE@/phase3/_restricted/population_source_mimic.csv' CSV HEADER
\copy (SELECT first_careunit AS careunit,count(*) AS stays FROM mimiciv_icu.icustays GROUP BY first_careunit ORDER BY first_careunit) TO '@CUFF_WORKSPACE@/phase3/population_dictionary_mimic_units.csv' CSV HEADER
\copy (SELECT sepsis3,count(*) AS records,count(DISTINCT stay_id) AS stays FROM mimiciv_derived.sepsis3 GROUP BY sepsis3 ORDER BY sepsis3) TO '@CUFF_WORKSPACE@/phase3/population_dictionary_mimic_sepsis.csv' CSV HEADER
ROLLBACK;
