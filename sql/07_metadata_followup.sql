-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
SET default_transaction_read_only=on;
\copy (SELECT table_name,column_name,data_type FROM information_schema.columns WHERE table_schema='public' AND table_name IN ('nurseassessment','nursecare','physicalexam','pastHistory','pasthistory')) TO '@CUFF_WORKSPACE@/phase0/eicu_assessment_schema.csv' CSV HEADER
\copy (SELECT nursingchartcelltypecat,nursingchartcelltypevallabel,nursingchartcelltypevalname,count(*) AS n_rows FROM public.nursecharting GROUP BY 1,2,3) TO '@CUFF_WORKSPACE@/phase0/eicu_nursing_dictionary_all.csv' CSV HEADER
\copy (SELECT treatmentstring,count(*) AS n_rows FROM public.treatment GROUP BY treatmentstring) TO '@CUFF_WORKSPACE@/phase0/eicu_treatment_dictionary_all.csv' CSV HEADER
