-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
CREATE TEMP TABLE recovery_audit_ids (infusiondrugid bigint);
\copy recovery_audit_ids FROM '@CUFF_WORKSPACE@/phase1/_restricted/recovery_audit_ids.csv' CSV HEADER
BEGIN TRANSACTION READ ONLY;
\copy (SELECT i.*,p.admissionweight FROM public.infusiondrug i JOIN recovery_audit_ids a USING(infusiondrugid) JOIN public.patient p USING(patientunitstayid)) TO '@CUFF_WORKSPACE@/phase1/_restricted/recovery_audit_infusion_source.csv' CSV HEADER
\copy (SELECT DISTINCT m.* FROM public.medication m JOIN public.infusiondrug i ON m.patientunitstayid=i.patientunitstayid JOIN recovery_audit_ids a ON a.infusiondrugid=i.infusiondrugid WHERE m.drugstartoffset<=i.infusionoffset AND m.drugstopoffset>i.infusionoffset AND m.drugname ~* '(norepi|noradren|levophed|epineph|adrenalin|dopamin|phenyleph|neo.?syneph|neosyneph|vasopress|pitressin|angiotensin|giapreza)') TO '@CUFF_WORKSPACE@/phase1/_restricted/recovery_audit_medication_source.csv' CSV HEADER
ROLLBACK;
