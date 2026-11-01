-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
SET default_transaction_read_only=on;
SET max_parallel_workers_per_gather=2;
\copy (SELECT i.*,p.gender,p.anchor_age+(extract(year from i.intime)::int-p.anchor_year) AS age,p.dod,a.hospital_expire_flag,a.deathtime,a.dischtime,row_number() OVER(PARTITION BY i.subject_id ORDER BY i.intime,i.stay_id) AS patient_icu_rank FROM mimiciv_icu.icustays i JOIN mimiciv_hosp.patients p USING(subject_id) JOIN mimiciv_hosp.admissions a USING(hadm_id)) TO '@CUFF_WORKSPACE@/phase0/_restricted/mimic_stays.csv' CSV HEADER
\copy (SELECT stay_id,itemid,starttime,endtime,rate,rateuom,amount,amountuom,patientweight,orderid,linkorderid,ordercategoryname,ordercategorydescription,statusdescription FROM mimiciv_icu.inputevents WHERE itemid IN ('221289','221662','221749','221906','222315','229617','229630','229631','229632','229709','229764')) TO '@CUFF_WORKSPACE@/phase0/_restricted/mimic_drugs.csv' CSV HEADER
\copy (SELECT stay_id,charttime,itemid,value,valuenum,valueuom FROM mimiciv_icu.chartevents WHERE itemid IN ('220050','220051','220052','220179','220180','220181','225309','225310','225312','220048','224639','226512','226531','226707','226730')) TO '@CUFF_WORKSPACE@/phase0/_restricted/mimic_bp_cov.csv' CSV HEADER
\copy (SELECT o.* FROM mimiciv_hosp.omr o WHERE result_name ~* 'weight|height') TO '@CUFF_WORKSPACE@/phase0/_restricted/mimic_omr.csv' CSV HEADER
\copy (SELECT DISTINCT x.hadm_id FROM mimiciv_hosp.diagnoses_icd x JOIN mimiciv_hosp.d_icd_diagnoses d USING(icd_code,icd_version) WHERE d.long_title ~* 'atrial fibrillation') TO '@CUFF_WORKSPACE@/phase0/_restricted/mimic_af_dx.csv' CSV HEADER
\copy (SELECT p.stay_id,p.itemid,p.starttime,p.endtime FROM mimiciv_icu.procedureevents p JOIN mimiciv_icu.d_items d USING(itemid) WHERE d.label ~* 'arterial.*line|arterial.*catheter|art.*line') TO '@CUFF_WORKSPACE@/phase0/_restricted/mimic_arterial_line.csv' CSV HEADER
