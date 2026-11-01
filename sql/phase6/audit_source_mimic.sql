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
\copy audit_sample FROM '@CUFF_WORKSPACE@/phase6/_restricted/audit_source_sample_mimic.csv' CSV HEADER
CREATE TEMP TABLE audit_map ("itemid" text,"label" text,"abbreviation" text,"category" text,"unitname" text,"linksto" text,"classification" text,"component" text,"family" text,"selected" text,"rows" text,"human_review" text,"count_source" text,"count_status" text);
\copy audit_map FROM '@CUFF_WORKSPACE@/phase0/map_bp_mimic.csv' CSV HEADER
BEGIN READ ONLY;
SET LOCAL statement_timeout=0;
\copy (SELECT a.audit_id,extract(epoch FROM c.charttime) source_t,m.classification modality,           m.family,m.component,c.valuenum value,c.itemid,c.valueuom           FROM audit_sample a JOIN mimiciv_icu.chartevents c ON c.stay_id=a.sid::bigint           JOIN audit_map m ON m.itemid::bigint=c.itemid AND m.selected='1'           WHERE c.charttime BETWEEN timestamp '1970-01-01'+(a.t::double precision-300)*interval '1 second'              AND timestamp '1970-01-01'+(a.t::double precision+300)*interval '1 second'           AND (m.classification='arterial' OR (m.classification='cuff' AND             c.charttime=timestamp '1970-01-01'+a.t::double precision*interval '1 second'))) TO '@CUFF_WORKSPACE@/phase6/_restricted/audit_source_rows_mimic.csv' CSV HEADER
COMMIT;
