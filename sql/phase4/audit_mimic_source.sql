-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
SET TIME ZONE 'UTC';
CREATE TEMP TABLE sampled_anchors(cid text,sid bigint,t double precision,expected_site text,expected_reason text,art_map double precision,cuff_map double precision,nee double precision);
\copy sampled_anchors FROM '@CUFF_WORKSPACE@/phase4/_restricted/site_audit_mimic.csv' CSV HEADER
BEGIN READ ONLY;
\copy (SELECT a.cid,a.sid,a.t,p.itemid,d.label,p.orderid,extract(epoch FROM p.starttime) start_seconds,extract(epoch FROM p.endtime) end_seconds,p.location,p.locationcategory,p.statusdescription FROM sampled_anchors a LEFT JOIN mimiciv_icu.procedureevents p ON p.stay_id=a.sid AND p.itemid IN (SELECT itemid FROM mimiciv_icu.d_items WHERE label='Arterial Line' AND linksto='procedureevents') LEFT JOIN mimiciv_icu.d_items d ON p.itemid=d.itemid ORDER BY a.cid,p.starttime,p.orderid) TO '@CUFF_WORKSPACE@/phase4/_restricted/audit_source_mimic_intervals.csv' CSV HEADER
\copy (SELECT a.cid,c.stay_id sid,c.itemid,d.label,extract(epoch FROM c.charttime) bp_seconds,c.valuenum FROM sampled_anchors a JOIN mimiciv_icu.chartevents c ON c.stay_id=a.sid AND c.charttime BETWEEN timestamp '1970-01-01'+(a.t-300)*interval '1 second' AND timestamp '1970-01-01'+(a.t+300)*interval '1 second' JOIN mimiciv_icu.d_items d ON d.itemid=c.itemid WHERE d.label IN ('Arterial Blood Pressure mean','ART BP Mean') ORDER BY a.cid,c.charttime,c.itemid) TO '@CUFF_WORKSPACE@/phase4/_restricted/audit_source_mimic_bp.csv' CSV HEADER
ROLLBACK;
