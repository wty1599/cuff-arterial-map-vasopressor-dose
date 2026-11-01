-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
SET default_transaction_read_only=on;
\copy (SELECT labname,labmeasurenameinterface,labmeasurenamesystem,count(*) AS n_rows,count(labresultoffset) AS timed_rows FROM public.lab WHERE labname ~* 'creatin|kreatin' GROUP BY 1,2,3) TO '@CUFF_WORKSPACE@/phase0/eicu_creatinine_dictionary_presence.csv' CSV HEADER
