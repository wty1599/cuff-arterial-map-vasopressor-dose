-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
SET default_transaction_read_only=on;
SET max_parallel_workers_per_gather=2;
\copy (SELECT itemid,count(*) AS n_rows FROM mimiciv_icu.chartevents WHERE itemid IN (220050,220051,220052,220056,220058,220179,220180,220181,220224,220227,220235,223751,223752,223830,224167,224289,224291,224404,224643,224828,225210,225309,225310,225312,225341,225556,225575,225698,225722,225737,225986,226008,226107,226861,227239,227240,227242,227243,227292,227537,227538,227539,227543,227546,227547,228022,228023,228026,229853) GROUP BY itemid) TO '@CUFF_WORKSPACE@/phase0/mimic_bp_dictionary_source_counts.csv' CSV HEADER
\copy (SELECT * FROM mimiciv_hosp.d_labitems WHERE label ~* 'creatin|kreatin') TO '@CUFF_WORKSPACE@/phase0/mimic_creatinine_dictionary.csv' CSV HEADER
\copy (SELECT EXISTS(SELECT 1 FROM mimiciv_hosp.labevents l JOIN mimiciv_hosp.d_labitems d USING(itemid) WHERE d.label ~* '^creatinine$' AND d.fluid='Blood' AND l.valuenum IS NOT NULL AND l.charttime IS NOT NULL) AS timed_serum_creatinine_rows_exist) TO '@CUFF_WORKSPACE@/phase0/mimic_creatinine_presence.csv' CSV HEADER
