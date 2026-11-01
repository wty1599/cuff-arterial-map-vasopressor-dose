-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
SET client_min_messages=warning;
SET TIME ZONE 'UTC';
CREATE TEMP TABLE audit_input(cid text,sid text,t text);
\copy audit_input FROM '@CUFF_WORKSPACE@/phase2/_restricted/audit_mimic_all_pair_keys.csv' CSV HEADER
CREATE TEMP TABLE audit_targets AS SELECT cid::bigint cid,sid::bigint stay_id,timestamp '1970-01-01'+t::double precision*interval '1 second' anchor FROM audit_input;
CREATE INDEX ON audit_targets(stay_id,anchor);
BEGIN READ ONLY;
SET LOCAL statement_timeout=0;
SET LOCAL work_mem='16MB';
SET LOCAL max_parallel_workers_per_gather=0;
\copy (WITH waveform AS MATERIALIZED (       SELECT c.stay_id,c.charttime,c.value       FROM mimiciv_icu.chartevents c JOIN mimiciv_icu.d_items d ON c.itemid=d.itemid       JOIN (SELECT DISTINCT stay_id FROM audit_targets) s USING(stay_id)       WHERE d.label='Arterial line Waveform Appear'     )     SELECT a.cid,a.stay_id AS sid,extract(epoch FROM a.anchor) AS t,       count(w.charttime) AS quality_records,       count(w.charttime) FILTER(WHERE w.value='Positional') AS positional_records,       min(abs(extract(epoch FROM(w.charttime-a.anchor)))) FILTER(WHERE w.value='Positional') AS nearest_positional_diff_seconds     FROM audit_targets a LEFT JOIN waveform w ON w.stay_id=a.stay_id       AND w.charttime BETWEEN a.anchor-interval '5 minutes' AND a.anchor+interval '5 minutes'     GROUP BY a.cid,a.stay_id,a.anchor) TO '@CUFF_WORKSPACE@/phase2/_restricted/audit_mimic_waveform_all_pairs.csv' CSV HEADER
COMMIT;
