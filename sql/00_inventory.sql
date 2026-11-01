-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
SET default_transaction_read_only = on;
SELECT current_database() AS database, version() AS version;
SELECT table_schema, table_name, column_name, data_type
FROM information_schema.columns
WHERE table_schema NOT IN ('pg_catalog','information_schema')
AND table_name IN ('patients','admissions','icustays','d_items','inputevents','chartevents','omr','labevents','diagnoses_icd','d_icd_diagnoses','patient','infusiondrug','vitalperiodic','vitalaperiodic','nursecharting','diagnosis','lab','treatment')
ORDER BY table_schema,table_name,ordinal_position;
SELECT schemaname,tablename,indexdef FROM pg_indexes
WHERE tablename IN ('chartevents','inputevents','vitalperiodic','vitalaperiodic','infusiondrug','nursecharting');
