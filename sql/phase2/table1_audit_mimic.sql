-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
CREATE TEMP TABLE a(sid text,patient text,sex text,hospital_death text,severity text,admission_type text,max_nee_observed text,max_nee_time text,max_nee_status text);
\copy a FROM '@CUFF_WORKSPACE@/phase2/_restricted/table1_audit_expected_mimic.csv' CSV HEADER
CREATE TEMP TABLE dm(itemid text,drug text);
\copy dm FROM '@CUFF_WORKSPACE@/phase2/_restricted/table1_audit_drugmap_mimic.csv' CSV HEADER
BEGIN TRANSACTION READ ONLY;
\copy (SELECT i.stay_id AS sid,p.gender,h.hospital_expire_flag,h.admission_type,s.sofa FROM a JOIN mimiciv_icu.icustays i ON i.stay_id=a.sid::bigint JOIN mimiciv_hosp.patients p USING(subject_id) JOIN mimiciv_hosp.admissions h ON h.hadm_id=i.hadm_id LEFT JOIN mimiciv_derived.first_day_sofa s ON s.stay_id=i.stay_id) TO '@CUFF_WORKSPACE@/phase2/_restricted/table1_audit_source_mimic_baseline.csv' CSV HEADER
\copy (SELECT x.stay_id AS sid,x.orderid,x.itemid,x.starttime,x.endtime,x.rate,x.rateuom,x.patientweight,x.ordercategorydescription,x.statusdescription,dm.drug FROM a JOIN mimiciv_icu.inputevents x ON x.stay_id=a.sid::bigint JOIN dm ON x.itemid::text=dm.itemid WHERE x.rate>0 AND x.ordercategorydescription='Continuous Med' AND x.statusdescription NOT IN ('Rewritten','Bolus')) TO '@CUFF_WORKSPACE@/phase2/_restricted/table1_audit_source_mimic_drugs.csv' CSV HEADER
ROLLBACK;
