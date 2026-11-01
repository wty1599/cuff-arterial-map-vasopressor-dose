-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
SET default_transaction_read_only=on;
\copy (SELECT cellattributepath,celllabel,cellattribute,count(*) AS n_rows FROM public.nurseassessment GROUP BY 1,2,3) TO '@CUFF_WORKSPACE@/phase0/eicu_assessment_dictionary_all.csv' CSV HEADER
\copy (SELECT patientunitstayid,nurseassessoffset,celllabel,cellattribute,cellattributevalue FROM public.nurseassessment WHERE concat_ws(' ',cellattributepath,celllabel,cellattribute) ~* 'rhythm|arrhyth|cardiac') TO '@CUFF_WORKSPACE@/phase0/_restricted/eicu_rhythm_assessment.csv' CSV HEADER
