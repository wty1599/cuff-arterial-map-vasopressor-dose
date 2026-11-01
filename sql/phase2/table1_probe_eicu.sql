-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
BEGIN READ ONLY;
\copy (SELECT table_schema,table_name,column_name,data_type FROM information_schema.columns WHERE table_name IN ('apachepatientresult','apacheapsvar','apachepredvar')) TO '@CUFF_WORKSPACE@/phase2/table1_schema_eicu.csv' CSV HEADER
ROLLBACK;
