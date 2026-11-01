-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\copy (SELECT patientunitstayid::text sid, count(*)>0 diagnosis_covered FROM public.diagnosis GROUP BY patientunitstayid) TO '@CUFF_WORKSPACE@/phase1/_restricted/eicu_diagnosis_coverage.csv' CSV HEADER
