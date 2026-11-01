-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
CREATE TEMP TABLE a(sid text,patient text,sex text,hospital_death text,severity text,admission_type text,max_nee_observed text,max_nee_time text,max_nee_status text);
\copy a FROM '@CUFF_WORKSPACE@/phase2/_restricted/table1_audit_expected_eicu.csv' CSV HEADER
CREATE TEMP TABLE dm(itemid text,drug text);
\copy dm FROM '@CUFF_WORKSPACE@/phase2/_restricted/table1_audit_drugmap_eicu.csv' CSV HEADER
BEGIN TRANSACTION READ ONLY;
\copy (SELECT p.patientunitstayid AS sid,p.gender,p.hospitaldischargestatus,p.unitstaytype,p.admissionweight,r.apachescore,r.apacheversion FROM a JOIN public.patient p ON p.patientunitstayid=a.sid::bigint LEFT JOIN public.apachepatientresult r ON r.patientunitstayid=p.patientunitstayid AND r.apacheversion='IVa') TO '@CUFF_WORKSPACE@/phase2/_restricted/table1_audit_source_eicu_baseline.csv' CSV HEADER
\copy (SELECT x.*,dm.drug FROM a JOIN public.infusiondrug x ON x.patientunitstayid=a.sid::bigint JOIN dm ON x.drugname=dm.itemid) TO '@CUFF_WORKSPACE@/phase2/_restricted/table1_audit_source_eicu_drugs.csv' CSV HEADER
ROLLBACK;
