-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
SET default_transaction_read_only = on;
SET max_parallel_workers_per_gather=2;
\copy (SELECT drugname, count(*) AS rows FROM public.infusiondrug GROUP BY drugname ORDER BY drugname) TO '@CUFF_WORKSPACE@/phase0/eicu_drug_dictionary_all.csv' CSV HEADER
\copy (SELECT nursingchartcelltypecat,nursingchartcelltypevallabel,nursingchartcelltypevalname, count(*) AS rows FROM public.nursecharting WHERE concat_ws(' ',nursingchartcelltypevallabel,nursingchartcelltypevalname) ~* '(rhythm|blood.?pressure|\mMAP\M|systol|diastol|non.?invasive|arterial)' GROUP BY 1,2,3 ORDER BY 1,2,3) TO '@CUFF_WORKSPACE@/phase0/eicu_nursing_dictionary_candidates.csv' CSV HEADER
\copy (SELECT 'icustays' AS stage,count(*) AS n FROM public.patient) TO '@CUFF_WORKSPACE@/phase0/inventory_counts_eicu.csv' CSV HEADER
