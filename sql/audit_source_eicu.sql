-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
SET client_min_messages=warning;
SET TIME ZONE 'UTC';
CREATE TEMP TABLE audit_sample ("audit_id" text,"sid" text,"t" text,"cuff_map" text,"art_map" text,"nee" text,"nee_status" text,"nearest_diff" text,"nearest_art_t" text,"art_values" text,"episode_start" text,"episode_stop" text);
\copy audit_sample FROM '@CUFF_WORKSPACE@/phase0/_restricted/audit_input_eicu.csv' CSV HEADER
CREATE TEMP TABLE audit_bpmap ("itemid" text,"label" text,"classification" text,"component" text,"unitname" text,"selected" text,"rows" text,"human_review" text);
\copy audit_bpmap FROM '@CUFF_WORKSPACE@/phase0/map_bp_eicu.csv' CSV HEADER
CREATE TEMP TABLE audit_drugmap ("itemid" text,"label" text,"drug" text,"unitname" text,"selected" text,"rows" text,"human_review" text);
\copy audit_drugmap FROM '@CUFF_WORKSPACE@/phase0/map_vaso_eicu.csv' CSV HEADER
BEGIN READ ONLY;
SET LOCAL statement_timeout=0;
\copy (SELECT a.audit_id,v.patientunitstayid AS stay_id,                 v.observationoffset::double precision*60 AS source_t,'arterial' AS modality,                 z.component,'art1' AS family,                 z.value,'mmHg' AS valueuom FROM audit_sample a                 JOIN public.vitalperiodic v ON v.patientunitstayid=a.sid::bigint                 CROSS JOIN LATERAL (VALUES ('map',v.systemicmean::double precision),('sbp',v.systemicsystolic::double precision),('dbp',v.systemicdiastolic::double precision)) z(component,value) WHERE v.observationoffset::double precision*60 BETWEEN a.t::double precision-300 AND a.t::double precision+300 UNION ALL SELECT a.audit_id,v.patientunitstayid AS stay_id,                 v.observationoffset::double precision*60 AS source_t,'cuff' AS modality,                 z.component,'cuff' AS family,                 z.value,'mmHg' AS valueuom FROM audit_sample a                 JOIN public.vitalaperiodic v ON v.patientunitstayid=a.sid::bigint                 CROSS JOIN LATERAL (VALUES ('map',v.noninvasivemean::double precision),('sbp',v.noninvasivesystolic::double precision),('dbp',v.noninvasivediastolic::double precision)) z(component,value) WHERE v.observationoffset::double precision*60=a.t::double precision) TO '@CUFF_WORKSPACE@/phase0/_restricted/audit_source_eicu_bp.csv' CSV HEADER
\copy (SELECT p.patientunitstayid AS stay_id,p.unitdischargeoffset,p.admissionweight             FROM public.patient p JOIN (SELECT DISTINCT sid::bigint AS patientunitstayid FROM audit_sample) a USING(patientunitstayid)) TO '@CUFF_WORKSPACE@/phase0/_restricted/audit_source_eicu_stays.csv' CSV HEADER
\copy (SELECT x.infusiondrugid,x.patientunitstayid AS stay_id,x.infusionoffset,             x.drugname,x.drugrate,x.infusionrate,x.drugamount,x.volumeoffluid,m.drug             FROM public.infusiondrug x JOIN (SELECT DISTINCT sid::bigint AS patientunitstayid FROM audit_sample) a USING(patientunitstayid)             JOIN (SELECT DISTINCT itemid,drug FROM audit_drugmap WHERE selected='1') m ON x.drugname=m.itemid) TO '@CUFF_WORKSPACE@/phase0/_restricted/audit_source_eicu_drugs.csv' CSV HEADER
COMMIT;
