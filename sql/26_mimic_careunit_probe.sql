-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
SET default_transaction_read_only=on;
SET max_parallel_workers_per_gather=2;
\copy (SELECT careunit,count(*) n_source_transfers FROM mimiciv_hosp.transfers GROUP BY careunit ORDER BY careunit) TO '@CUFF_WORKSPACE@/phase0/mimic_transfer_careunit_dictionary.csv' CSV HEADER
\copy (SELECT t.transfer_id,t.hadm_id,t.careunit,t.intime,t.outtime FROM mimiciv_hosp.transfers t WHERE t.intime IS NOT NULL AND EXISTS(SELECT 1 FROM mimiciv_icu.icustays i WHERE i.hadm_id=t.hadm_id AND t.intime<i.outtime AND (t.outtime>i.intime OR t.outtime IS NULL))) TO '@CUFF_WORKSPACE@/phase0/_restricted/mimic_transfer_careunit_probe.csv' CSV HEADER
