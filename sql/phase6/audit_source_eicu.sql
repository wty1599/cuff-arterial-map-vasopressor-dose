-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
SET client_min_messages=warning;
SET TIME ZONE 'UTC';
CREATE TEMP TABLE audit_sample ("audit_id" text,"sid" text,"t" text,"episode_start" text,"episode_stop" text,"cuff_map" text,"art_map" text,"art_2" text,"art_same" text,"dose_stratum" text,"changed" text);
\copy audit_sample FROM '@CUFF_WORKSPACE@/phase6/_restricted/audit_source_sample_eicu.csv' CSV HEADER
BEGIN READ ONLY;
SET LOCAL statement_timeout=0;
\copy (SELECT a.audit_id,v.observationoffset::double precision*60 source_t,'arterial' modality,'arterial' family,'map' component,v.systemicmean::double precision value,'vitalperiodic.systemicmean' itemid,'mmHg' valueuom FROM audit_sample a JOIN public.vitalperiodic v ON v.patientunitstayid=a.sid::bigint WHERE v.observationoffset::double precision*60 BETWEEN a.t::double precision-300 AND a.t::double precision+300 UNION ALL SELECT a.audit_id,v.observationoffset::double precision*60 source_t,'arterial' modality,'arterial' family,'sbp' component,v.systemicsystolic::double precision value,'vitalperiodic.systemicsystolic' itemid,'mmHg' valueuom FROM audit_sample a JOIN public.vitalperiodic v ON v.patientunitstayid=a.sid::bigint WHERE v.observationoffset::double precision*60 BETWEEN a.t::double precision-300 AND a.t::double precision+300 UNION ALL SELECT a.audit_id,v.observationoffset::double precision*60 source_t,'arterial' modality,'arterial' family,'dbp' component,v.systemicdiastolic::double precision value,'vitalperiodic.systemicdiastolic' itemid,'mmHg' valueuom FROM audit_sample a JOIN public.vitalperiodic v ON v.patientunitstayid=a.sid::bigint WHERE v.observationoffset::double precision*60 BETWEEN a.t::double precision-300 AND a.t::double precision+300 UNION ALL SELECT a.audit_id,v.observationoffset::double precision*60 source_t,'cuff' modality,'cuff' family,'map' component,v.noninvasivemean::double precision value,'vitalaperiodic.noninvasivemean' itemid,'mmHg' valueuom FROM audit_sample a JOIN public.vitalaperiodic v ON v.patientunitstayid=a.sid::bigint WHERE v.observationoffset::double precision*60=a.t::double precision UNION ALL SELECT a.audit_id,v.observationoffset::double precision*60 source_t,'cuff' modality,'cuff' family,'sbp' component,v.noninvasivesystolic::double precision value,'vitalaperiodic.noninvasivesystolic' itemid,'mmHg' valueuom FROM audit_sample a JOIN public.vitalaperiodic v ON v.patientunitstayid=a.sid::bigint WHERE v.observationoffset::double precision*60=a.t::double precision UNION ALL SELECT a.audit_id,v.observationoffset::double precision*60 source_t,'cuff' modality,'cuff' family,'dbp' component,v.noninvasivediastolic::double precision value,'vitalaperiodic.noninvasivediastolic' itemid,'mmHg' valueuom FROM audit_sample a JOIN public.vitalaperiodic v ON v.patientunitstayid=a.sid::bigint WHERE v.observationoffset::double precision*60=a.t::double precision) TO '@CUFF_WORKSPACE@/phase6/_restricted/audit_source_rows_eicu.csv' CSV HEADER
COMMIT;
