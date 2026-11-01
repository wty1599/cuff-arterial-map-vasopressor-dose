-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
CREATE TEMP TABLE table1_ids (sid bigint);
\copy table1_ids FROM '@CUFF_WORKSPACE@/phase2/_restricted/table1_ids_eicu.csv' CSV HEADER
BEGIN TRANSACTION READ ONLY;
\copy (SELECT p.patientunitstayid AS sid,p.gender,p.unittype,p.unitadmitsource,p.unitstaytype,p.hospitaldischargestatus,r.apachepatientresultsid,r.apacheversion,r.apachescore,r.acutephysiologyscore,r.actualhospitalmortality FROM table1_ids k JOIN public.patient p ON p.patientunitstayid=k.sid LEFT JOIN public.apachepatientresult r ON r.patientunitstayid=p.patientunitstayid) TO '@CUFF_WORKSPACE@/phase2/_restricted/table1_eicu_source.csv' CSV HEADER
\copy (SELECT apacheversion,count(*) AS records FROM public.apachepatientresult GROUP BY apacheversion ORDER BY apacheversion) TO '@CUFF_WORKSPACE@/phase2/table1_apache_version_dictionary.csv' CSV HEADER
ROLLBACK;
