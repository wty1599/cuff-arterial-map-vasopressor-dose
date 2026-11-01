-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
SET client_min_messages=warning;
SET TIME ZONE 'UTC';
CREATE TEMP TABLE audit_sample ("cid" text,"sid" text,"t" text,"cuff_map" text,"win" text,"art_map" text,"art_values" text,"nearest_diff" text,"nearest_art_t" text,"art_ends0or5_fraction" text,"art_min_t" text,"art_max_t" text,"pid" text,"is_primary" text,"adult" text,"hospital" text,"bmi" text,"af" text,"af_source" text,"nee" text,"nee_status" text,"age" text,"weight_raw" text,"height_raw" text,"database" text,"patient" text,"bmi_raw" text,"nee_raw" text,"weight_flag" text,"weight" text,"height_flag" text,"height" text,"bmi_flag" text,"nee_flag" text,"dose_weight_flag" text,"quality_any" text,"dose_stratum" text,"bmi40" text,"sex" text,"af_observed" text,"af_model" text,"af_subgroup" text,"audit_hash" text,"episode_start" text,"episode_stop" text,"audit_id" text,"expected_phase0_weight" text,"expected_phase0_height" text);
\copy audit_sample FROM '@CUFF_WORKSPACE@/phase3/icm_revision4/analysis/_restricted/audit_sample_mimic.csv' CSV HEADER
CREATE TEMP TABLE audit_bpmap ("itemid" text,"label" text,"abbreviation" text,"category" text,"unitname" text,"linksto" text,"classification" text,"component" text,"family" text,"selected" text,"rows" text,"human_review" text,"count_source" text,"count_status" text);
\copy audit_bpmap FROM '@CUFF_WORKSPACE@/phase0/map_bp_mimic.csv' CSV HEADER
CREATE TEMP TABLE audit_drugmap ("itemid" text,"label" text,"drug" text,"unitname" text,"linksto" text,"selected" text,"rows" text,"human_review" text,"count_source" text,"count_status" text);
\copy audit_drugmap FROM '@CUFF_WORKSPACE@/phase0/map_vaso_mimic.csv' CSV HEADER
CREATE TEMP TABLE audit_covmap ("itemid" text,"label" text,"unitname" text);
\copy audit_covmap FROM '@CUFF_WORKSPACE@/phase0/map_covariates_mimic.csv' CSV HEADER
BEGIN READ ONLY;
SET LOCAL statement_timeout=0;
\copy (SELECT i.stay_id,i.subject_id,i.hadm_id,i.intime,i.outtime             FROM mimiciv_icu.icustays i JOIN (SELECT DISTINCT sid::bigint AS stay_id FROM audit_sample) a USING(stay_id)) TO '@CUFF_WORKSPACE@/phase3/icm_revision4/analysis/_restricted/audit_source_mimic_stays.csv' CSV HEADER
\copy (SELECT a.audit_id,c.stay_id,extract(epoch FROM c.charttime) AS source_t,             c.itemid,m.classification AS modality,m.component,m.family,c.valuenum AS value,c.valueuom             FROM audit_sample a JOIN mimiciv_icu.chartevents c ON c.stay_id=a.sid::bigint             JOIN (SELECT DISTINCT itemid,classification,component,family FROM audit_bpmap WHERE selected='1') m               ON c.itemid=m.itemid::bigint             WHERE c.charttime >= timestamp '1970-01-01'+(a.t::double precision-300)*interval '1 second'               AND c.charttime <= timestamp '1970-01-01'+(a.t::double precision+300)*interval '1 second'               AND (m.classification='arterial' OR (m.classification='cuff' AND                  c.charttime=timestamp '1970-01-01'+a.t::double precision*interval '1 second'))) TO '@CUFF_WORKSPACE@/phase3/icm_revision4/analysis/_restricted/audit_source_mimic_bp.csv' CSV HEADER
\copy (SELECT x.stay_id,m.drug,x.itemid,x.starttime,x.endtime,x.rate,x.rateuom,             x.patientweight,x.orderid,x.linkorderid,x.ordercategoryname,x.ordercategorydescription,x.statusdescription             FROM mimiciv_icu.inputevents x JOIN (SELECT DISTINCT sid::bigint AS stay_id FROM audit_sample) a USING(stay_id)             JOIN (SELECT DISTINCT itemid,drug FROM audit_drugmap WHERE selected='1') m ON x.itemid=m.itemid::bigint) TO '@CUFF_WORKSPACE@/phase3/icm_revision4/analysis/_restricted/audit_source_mimic_drugs.csv' CSV HEADER
\copy (SELECT c.stay_id,extract(epoch FROM c.charttime) AS source_t,c.itemid,             c.valuenum AS value,c.valueuom,m.label             FROM mimiciv_icu.chartevents c JOIN (SELECT DISTINCT sid::bigint AS stay_id FROM audit_sample) a USING(stay_id)             JOIN audit_covmap m ON c.itemid=m.itemid::bigint             WHERE m.label IN ('Admission Weight (Kg)','Admission Weight (lbs.)','Daily Weight') OR lower(m.label) LIKE '%height%' ) TO '@CUFF_WORKSPACE@/phase3/icm_revision4/analysis/_restricted/audit_source_mimic_weights.csv' CSV HEADER
\copy (SELECT o.subject_id,o.chartdate,o.seq_num,o.result_name,o.result_value             FROM mimiciv_hosp.omr o JOIN (SELECT DISTINCT i.subject_id FROM mimiciv_icu.icustays i                 JOIN (SELECT DISTINCT sid::bigint AS stay_id FROM audit_sample) a USING(stay_id)) s USING(subject_id)             WHERE o.result_name ~* 'weight' ) TO '@CUFF_WORKSPACE@/phase3/icm_revision4/analysis/_restricted/audit_source_mimic_omr.csv' CSV HEADER
COMMIT;
