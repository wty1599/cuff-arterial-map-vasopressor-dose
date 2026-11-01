-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
SET default_transaction_read_only = on;
\copy (SELECT * FROM mimiciv_icu.d_items ORDER BY itemid) TO '@CUFF_WORKSPACE@/phase0/mimic_dictionary_all.csv' CSV HEADER
\copy (SELECT itemid,label,abbreviation,category,unitname,linksto FROM mimiciv_icu.d_items WHERE concat_ws(' ',label,abbreviation) ~* '(arterial|\mART\M|\mABP\M|\mNBP\M|\mNIBP\M|non.?invasive|blood.?pressure|systolic|diastolic|norepine|noradren|levophed|epinephr|adrenalin|phenylephr|neosyneph|dopamin|vasopress|pitressin|angiotensin|giapreza|weight|height|rhythm)' ORDER BY itemid) TO '@CUFF_WORKSPACE@/phase0/mimic_dictionary_candidates.csv' CSV HEADER
\copy (SELECT 'icustays' AS stage,count(*) AS n FROM mimiciv_icu.icustays UNION ALL SELECT 'patients',count(*) FROM mimiciv_hosp.patients) TO '@CUFF_WORKSPACE@/phase0/inventory_counts_mimic.csv' CSV HEADER
