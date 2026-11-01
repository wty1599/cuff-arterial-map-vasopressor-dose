-- Private-workspace guard: use analysis/run_private.py after preparation.
\if :{?CUFF_PREPARED_WORKSPACE}
\else
\echo 'Use the prepared-workspace launcher; extraction is disabled in the checkout.'
\quit 3
\endif

\set ON_ERROR_STOP on
SET timezone='UTC';
CREATE TEMP TABLE selected_stays(sid bigint);
\copy selected_stays FROM '@CUFF_WORKSPACE@/phase4/_restricted/site_stays_mimic.csv' CSV HEADER
BEGIN READ ONLY;
\copy (SELECT p.stay_id AS sid,extract(epoch FROM p.starttime) AS start,extract(epoch FROM p.endtime) AS stop,p.itemid,p.location,p.locationcategory,p.orderid,p.linkorderid,p.statusdescription FROM mimiciv_icu.procedureevents p JOIN selected_stays s ON p.stay_id=s.sid JOIN mimiciv_icu.d_items d USING(itemid) WHERE d.label='Arterial Line' AND d.linksto='procedureevents') TO '@CUFF_WORKSPACE@/phase4/_restricted/site_source_mimic.csv' CSV HEADER
\copy (SELECT location,locationcategory,count(*) AS records FROM mimiciv_icu.procedureevents p JOIN mimiciv_icu.d_items d USING(itemid) WHERE d.label='Arterial Line' AND d.linksto='procedureevents' GROUP BY 1,2 ORDER BY 1,2) TO '@CUFF_WORKSPACE@/phase4/site_location_dictionary_mimic.csv' CSV HEADER
\copy (SELECT d.itemid,d.label,c.value,count(*) AS records FROM mimiciv_icu.chartevents c JOIN selected_stays s ON c.stay_id=s.sid JOIN mimiciv_icu.d_items d USING(itemid) WHERE d.label ILIKE '%arterial%site%' GROUP BY 1,2,3 ORDER BY 1,3) TO '@CUFF_WORKSPACE@/phase4/site_chartevents_dictionary_mimic.csv' CSV HEADER
ROLLBACK;
