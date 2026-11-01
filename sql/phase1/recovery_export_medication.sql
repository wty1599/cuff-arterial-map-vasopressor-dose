-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
CREATE TEMP TABLE recovery_stays (sid bigint);
\copy recovery_stays FROM '@CUFF_WORKSPACE@/phase1/_restricted/recovery_candidate_stays.csv' CSV HEADER
BEGIN TRANSACTION READ ONLY;
\copy (SELECT m.* FROM public.medication m JOIN recovery_stays s ON s.sid=m.patientunitstayid WHERE m.drugname ~* '(norepi|noradren|levophed|epineph|adrenalin|dopamin|phenyleph|neo.?syneph|neosyneph|vasopress|pitressin|angiotensin|giapreza)') TO '@CUFF_WORKSPACE@/phase1/_restricted/recovery_medication_source.csv' CSV HEADER
ROLLBACK;
