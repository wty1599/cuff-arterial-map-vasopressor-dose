-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\copy (SELECT i.stay_id::text sid, count(d.seq_num)>0 diagnosis_covered FROM mimiciv_icu.icustays i LEFT JOIN mimiciv_hosp.diagnoses_icd d USING(hadm_id) GROUP BY i.stay_id) TO '@CUFF_WORKSPACE@/phase1/_restricted/mimic_diagnosis_coverage.csv' CSV HEADER
