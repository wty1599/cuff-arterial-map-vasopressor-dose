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
\copy audit_missing_sample FROM '@CUFF_WORKSPACE@/phase0/_restricted/missing_check_input_eicu.csv' CSV HEADER
ANALYZE audit_missing_sample;
BEGIN READ ONLY;
SET LOCAL statement_timeout=0;
\copy (WITH intervals AS (SELECT DISTINCT sid::bigint AS sid,target,episode_start::double precision AS a,                 episode_stop::double precision AS b FROM audit_missing_sample),             pt AS (SELECT sid,target,count(*) AS episode_rows FROM intervals GROUP BY sid,target), observations AS (SELECT v.patientunitstayid AS sid,'arterial' AS target,                 count(*) AS n_mapped_source_rows,count(v.systemicmean) AS n_source_numeric,                 count(*) FILTER (WHERE v.systemicmean BETWEEN 20 AND 200) AS n_source_in_range,                 count(v.systemicmean) AS n_source_nonblank_text                 FROM public.vitalperiodic v JOIN pt ON v.patientunitstayid=pt.sid AND pt.target='arterial'                 WHERE EXISTS (SELECT 1 FROM intervals i WHERE i.sid=v.patientunitstayid AND i.target='arterial'                     AND v.observationoffset::double precision*60>=i.a                     AND v.observationoffset::double precision*60<i.b)                 GROUP BY v.patientunitstayid UNION ALL SELECT v.patientunitstayid AS sid,'cuff' AS target,                 count(*) AS n_mapped_source_rows,count(v.noninvasivemean) AS n_source_numeric,                 count(*) FILTER (WHERE v.noninvasivemean BETWEEN 20 AND 200) AS n_source_in_range,                 count(v.noninvasivemean) AS n_source_nonblank_text                 FROM public.vitalaperiodic v JOIN pt ON v.patientunitstayid=pt.sid AND pt.target='cuff'                 WHERE EXISTS (SELECT 1 FROM intervals i WHERE i.sid=v.patientunitstayid AND i.target='cuff'                     AND v.observationoffset::double precision*60>=i.a                     AND v.observationoffset::double precision*60<i.b)                 GROUP BY v.patientunitstayid)             SELECT pt.sid,pt.target,pt.episode_rows,                 coalesce(o.n_mapped_source_rows,0) AS n_mapped_source_rows,                 coalesce(o.n_source_numeric,0) AS n_source_numeric,                 coalesce(o.n_source_in_range,0) AS n_source_in_range,                 coalesce(o.n_source_nonblank_text,0) AS n_source_nonblank_text             FROM pt LEFT JOIN observations o USING(sid,target) ORDER BY pt.target,pt.sid) TO '@CUFF_WORKSPACE@/phase0/_restricted/audit_missing_source_eicu_patient_targets.csv' CSV HEADER
COMMIT;
