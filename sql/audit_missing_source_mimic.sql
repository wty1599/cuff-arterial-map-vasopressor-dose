-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
SET client_min_messages=warning;
SET TIME ZONE 'UTC';
CREATE TEMP TABLE audit_missing_sample ("sid" text,"target" text,"episode_start" text,"episode_stop" text);
\copy audit_missing_sample FROM '@CUFF_WORKSPACE@/phase0/_restricted/missing_check_input_mimic.csv' CSV HEADER
ANALYZE audit_missing_sample;
BEGIN READ ONLY;
SET LOCAL statement_timeout=0;
\copy (SELECT * FROM ( SELECT itemid,label,abbreviation,category,unitname,  CASE  WHEN lower(label) ~ '(non[ -]?invasive|non invasive|\mnibp\M|\mnbp\M)'   AND lower(label) ~ '(mean|\mmap\M)' THEN 'cuff'  WHEN lower(label) ~ '(arterial.*(blood pressure|pressure)|\mart[ .-]*bp\M|\mabp\M)'   AND lower(label) ~ '(mean|\mmap\M)' THEN 'arterial'  ELSE NULL END AS target FROM mimiciv_icu.d_items WHERE linksto='chartevents'  AND lower(label) !~ '(alarm|target|iabp|pulmonary|\mpap\M|streef)' ) d WHERE target IS NOT NULL ORDER BY target,itemid) TO '@CUFF_WORKSPACE@/phase0/_restricted/audit_missing_source_mimic_dictionary.csv' CSV HEADER
\copy (WITH intervals AS (SELECT DISTINCT sid::bigint AS sid,target,episode_start::double precision AS a,                 episode_stop::double precision AS b FROM audit_missing_sample),             pt AS (SELECT sid,target,count(*) AS episode_rows FROM intervals GROUP BY sid,target), concepts AS ( SELECT itemid,label,abbreviation,category,unitname,  CASE  WHEN lower(label) ~ '(non[ -]?invasive|non invasive|\mnibp\M|\mnbp\M)'   AND lower(label) ~ '(mean|\mmap\M)' THEN 'cuff'  WHEN lower(label) ~ '(arterial.*(blood pressure|pressure)|\mart[ .-]*bp\M|\mabp\M)'   AND lower(label) ~ '(mean|\mmap\M)' THEN 'arterial'  ELSE NULL END AS target FROM mimiciv_icu.d_items WHERE linksto='chartevents'  AND lower(label) !~ '(alarm|target|iabp|pulmonary|\mpap\M|streef)' ), observations AS (             SELECT c.stay_id AS sid,d.target,                 count(*) AS n_mapped_source_rows,                 count(c.valuenum) AS n_source_numeric,                 count(*) FILTER (WHERE c.valuenum BETWEEN 20 AND 200) AS n_source_in_range,                 count(*) FILTER (WHERE nullif(btrim(c.value),'') IS NOT NULL) AS n_source_nonblank_text             FROM mimiciv_icu.chartevents c JOIN concepts d ON c.itemid=d.itemid             JOIN pt ON c.stay_id=pt.sid AND d.target=pt.target             WHERE d.target IS NOT NULL AND EXISTS (                 SELECT 1 FROM intervals i WHERE i.sid=c.stay_id AND i.target=d.target                 AND c.charttime>=timestamp '1970-01-01'+i.a*interval '1 second'                 AND c.charttime<timestamp '1970-01-01'+i.b*interval '1 second')             GROUP BY c.stay_id,d.target)             SELECT pt.sid,pt.target,pt.episode_rows,                 coalesce(o.n_mapped_source_rows,0) AS n_mapped_source_rows,                 coalesce(o.n_source_numeric,0) AS n_source_numeric,                 coalesce(o.n_source_in_range,0) AS n_source_in_range,                 coalesce(o.n_source_nonblank_text,0) AS n_source_nonblank_text             FROM pt LEFT JOIN observations o USING(sid,target) ORDER BY pt.target,pt.sid) TO '@CUFF_WORKSPACE@/phase0/_restricted/audit_missing_source_mimic_patient_targets.csv' CSV HEADER
COMMIT;
