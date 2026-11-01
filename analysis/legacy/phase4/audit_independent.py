# Private-workspace guard: run only after analysis/prepare_workspace.py.
from pathlib import Path as _ResearchPath
_research_workspace = _ResearchPath('@CUFF_WORKSPACE@')
if '@' in str(_research_workspace) or not _research_workspace.is_absolute():
    raise RuntimeError('Prepare an absolute private workspace before running analysis templates.')
_research_workspace = _research_workspace.resolve()
for _candidate in _ResearchPath(__file__).resolve().parents:
    if (_candidate/'analysis/prepare_workspace.py').is_file() and (_candidate/'outputs/aggregate').is_dir():
        if _research_workspace == _candidate or _candidate in _research_workspace.parents:
            raise RuntimeError('Patient-data workspace must remain outside the repository checkout.')
# End private-workspace guard.

import os as _cuff_env
"""Bounded independent source and aggregate audit. Never imports production functions."""
from pathlib import Path
import os,sys,csv,json,re,subprocess,hashlib,datetime,math
ROOT=Path('@CUFF_WORKSPACE@');OUT=ROOT/'phase4';RES=OUT/'_restricted';CODE=ROOT/'sql/phase4'
import numpy as np,pandas as pd,duckdb
def objects(d):
 d=d.copy()
 for col in d:
  if isinstance(d[col].dtype,pd.StringDtype):d[col]=d[col].astype(object)
 return d
def outcsv(d,name):pd.DataFrame(d).to_csv(OUT/name,index=False,encoding='utf-8-sig')
def sql_sources():
 q=r"""\set ON_ERROR_STOP on
SET TIME ZONE 'UTC';
CREATE TEMP TABLE sampled_anchors(cid text,sid bigint,t double precision,expected_site text,expected_reason text,art_map double precision,cuff_map double precision,nee double precision);
\copy sampled_anchors FROM '@CUFF_WORKSPACE@/phase4/_restricted/site_audit_mimic.csv' CSV HEADER
BEGIN READ ONLY;
\copy (SELECT a.cid,a.sid,a.t,p.itemid,d.label,p.orderid,extract(epoch FROM p.starttime) start_seconds,extract(epoch FROM p.endtime) end_seconds,p.location,p.locationcategory,p.statusdescription FROM sampled_anchors a LEFT JOIN mimiciv_icu.procedureevents p ON p.stay_id=a.sid AND p.itemid IN (SELECT itemid FROM mimiciv_icu.d_items WHERE label='Arterial Line' AND linksto='procedureevents') LEFT JOIN mimiciv_icu.d_items d ON p.itemid=d.itemid ORDER BY a.cid,p.starttime,p.orderid) TO '@CUFF_WORKSPACE@/phase4/_restricted/audit_source_mimic_intervals.csv' CSV HEADER
\copy (SELECT a.cid,c.stay_id sid,c.itemid,d.label,extract(epoch FROM c.charttime) bp_seconds,c.valuenum FROM sampled_anchors a JOIN mimiciv_icu.chartevents c ON c.stay_id=a.sid AND c.charttime BETWEEN timestamp '1970-01-01'+(a.t-300)*interval '1 second' AND timestamp '1970-01-01'+(a.t+300)*interval '1 second' JOIN mimiciv_icu.d_items d ON d.itemid=c.itemid WHERE d.label IN ('Arterial Blood Pressure mean','ART BP Mean') ORDER BY a.cid,c.charttime,c.itemid) TO '@CUFF_WORKSPACE@/phase4/_restricted/audit_source_mimic_bp.csv' CSV HEADER
ROLLBACK;
"""
 p=CODE/'audit_mimic_source.sql';p.write_text(q,encoding='utf-8')
 cmd=[_cuff_env.environ.get('CUFF_PSQL','psql'),'-X','-v','ON_ERROR_STOP=1','-d','mimiciv','-f','@CUFF_WORKSPACE@/sql/phase4/'+p.name]
 p=subprocess.run(cmd,capture_output=True);(RES/'audit_mimic_source.log').write_bytes(p.stdout+b'\n'+p.stderr)
 if p.returncode:raise RuntimeError('Independent source query failed; see restricted log')
 # Decode the original Amsterdam file, retaining all process types for the sampled admissions.
 a=pd.read_csv(RES/'site_audit_amsterdam.csv',dtype={'sid':str,'cid':str});ids=set(a.sid);parts=[];scanned=0
 for chunk in pd.read_csv('@AMSTERDAM_ROOT@/processitems.csv',encoding='cp1252',chunksize=100000,dtype={'admissionid':str}):
  scanned+=len(chunk);parts.append(chunk[chunk.admissionid.isin(ids)])
 source=pd.concat(parts,ignore_index=True);c=duckdb.connect();c.register('raw_processitems',objects(source));c.register('sampled_anchors',objects(a))
 query="SELECT a.cid,a.sid,a.t,p.itemid,p.item AS label,p.start/1000.0 AS start_seconds,p.stop/1000.0 AS end_seconds FROM sampled_anchors a LEFT JOIN raw_processitems p ON a.sid=p.admissionid AND lower(p.item) LIKE 'arterielijn %' ORDER BY a.cid,p.start,p.itemid"
 got=c.execute(query).df();c.close();got.to_csv(RES/'audit_source_amsterdam_intervals.csv',index=False)
 (CODE/'audit_amsterdam_source.sql').write_text('-- Executed with DuckDB over a direct cp1252 read of original @AMSTERDAM_ROOT@/processitems.csv; all process types retained before SQL. No production interval file used.\n'+query+';\n',encoding='utf-8')
 return scanned
def vessel(label):
 x=str(label).strip().casefold()
 for rx,name in [(r'radial','radial'),(r'femor','femoral'),(r'brachi','brachial'),(r'axill','axillary'),(r'ulnar','ulnar'),(r'dorsalis','dorsalis_pedis')]:
  if re.search(rx,x):return name
 return 'unknown'
def independently_classify(rows,t,db):
 active=[]
 for x in rows:
  try:s=float(x['start_seconds']);e=float(x['end_seconds'])
  except (ValueError,TypeError):continue
  if not (math.isfinite(s) and math.isfinite(e) and s<=t<=e):continue
  if db=='mimic' and x.get('statusdescription')=='Rewritten':continue
  label=x.get('location','') if db=='mimic' else x.get('label','');site=vessel(label)
  if db=='amsterdam' and site=='unknown':continue # Non-site process labels cannot establish an anatomical interval.
  key=str(label).strip().casefold() if db=='mimic' else site
  active.append((key if site!='unknown' else 'unknown',site))
 known={k for k,v in active if v!='unknown'}
 if len(known)>1:return 'mixed','concurrent_distinct_locations',len(active)
 if any(v=='unknown' for _,v in active):return 'unknown','active_line_without_anatomical_site',len(active)
 if not known:return 'unknown','no_valid_interval_at_anchor',0
 return next(v for _,v in active),'one_documented_site',len(active)
def source_checks():
 summary=[];allrows=[]
 for db in ['mimic','amsterdam']:
  a=pd.read_csv(RES/f'site_audit_{db}.csv',dtype={'cid':str,'sid':str});s=pd.read_csv(RES/f'audit_source_{db}_intervals.csv',dtype={'cid':str,'sid':str},keep_default_na=False);by={key:g.to_dict('records') for key,g in s.groupby('cid')}
  checked=[]
  for x in a.to_dict('records'):
   site,reason,n=independently_classify(by.get(x['cid'],[]),float(x['t']),db)
   checked.append(dict(cid=x['cid'],sid=x['sid'],expected_site=x['site'],source_site=site,expected_reason=x['site_reason'],source_reason=reason,active_source_rows=n,site_match=site==x['site'],reason_match=reason==x['site_reason']))
  pd.DataFrame(checked).to_csv(RES/f'audit_site_comparison_{db}.csv',index=False)
  summary.append(dict(database=db,audit='independent_source_catheter_intervals',sampled_anchors=len(a),sampled_stays=a.sid.nunique(),site_mismatches=sum(not r['site_match'] for r in checked),reason_mismatches=sum(not r['reason_match'] for r in checked),source='original PostgreSQL procedureevents joined to d_items' if db=='mimic' else 'original cp1252 processitems.csv; independent DuckDB SQL',scope='stratified bounded sample; not a population error estimate'))
  for site,g in pd.DataFrame(checked).groupby('expected_site'):allrows.append(dict(database=db,site=site,anchors=len(g),site_mismatches=int((~g.site_match).sum())))
 # Classify BP source timestamps, not their nearby cuff anchors, for item/site cross-tab.
 s=pd.read_csv(RES/'audit_source_mimic_intervals.csv',dtype={'sid':str,'cid':str},keep_default_na=False);by={k:g.to_dict('records') for k,g in s.groupby('sid')}
 b=pd.read_csv(RES/'audit_source_mimic_bp.csv',dtype={'sid':str,'cid':str}).drop_duplicates(['sid','itemid','bp_seconds','valuenum'])
 b['source_time_site']=[independently_classify(by.get(r.sid,[]),r.bp_seconds,'mimic')[0] for r in b.itertuples()]
 cross=b.groupby(['itemid','label','source_time_site'],dropna=False).agg(observations=('sid','size'),stays=('sid','nunique')).reset_index()
 outcsv(cross,'phase4_abp_item_site_audit.csv');outcsv(allrows,'phase4_source_audit_strata.csv');outcsv(summary,'phase4_source_audit_summary.csv')
 return summary,cross.to_dict('records')
def aggregate_checks():
 checks=[];rows=[];quality=[];S4=pd.read_csv(OUT/'phase4_S4_bias_by_range.csv');S3=pd.read_csv(OUT/'phase4_S3_site.csv');compare=pd.read_csv(OUT/'phase4_S3_comparison.csv');baseline=pd.read_csv(OUT/'phase4_S3_baseline.csv');coverage=pd.read_csv(OUT/'phase4_site_availability.csv');distribution=pd.read_csv(OUT/'phase4_site_distribution.csv')
 for db in ['mimic','eicu']:
  p=(ROOT/f'phase1/_restricted/pairs_{db}.parquet').as_posix();c=duckdb.connect()
  for axis,expr in [('arterial_MAP','art_map'),('mean_of_cuff_and_arterial_MAP','(art_map+cuff_map)/2.0')]:
   sql=f"WITH d AS (SELECT *,CASE WHEN {expr}<55 THEN '<55' WHEN {expr}<65 THEN '55-64' WHEN {expr}<75 THEN '65-74' WHEN {expr}<85 THEN '75-84' WHEN {expr}<95 THEN '85-94' ELSE '>=95' END AS bin FROM read_parquet('{p}')) SELECT bin,count(*) n,count(DISTINCT patient) patients,avg(cuff_map-art_map) delta FROM d GROUP BY bin"
   got=c.execute(sql).df();exp=S4[S4.database.eq(db)&S4.axis.eq(axis)].set_index('pressure_range')
   for r in got.itertuples():
    e=exp.loc[r.bin];passed=int(e.denominator)==r.n and int(e.denominator_patients)==r.patients and abs(e.estimate-r.delta)<1e-10;rows.append(dict(database=db,axis=axis,pressure_range=r.bin,pairs=r.n,source_mean=r.delta,published_mean=e.estimate,absolute_difference=abs(e.estimate-r.delta),passed=passed))
  c.close()
  orig=pd.read_parquet(p);checks.append(dict(database=db,check='paired_MAP_values_finite',passed=np.isfinite(orig[['art_map','cuff_map']].to_numpy()).all()))
  a=S4[S4.database.eq(db)].reset_index(drop=True);rep=np.load(RES/f'bias_bootstrap_{db}.npz')['estimates'];ci=np.nanquantile(rep,[.025,.975],axis=0)
  checks.append(dict(database=db,check='S4_saved_replica_CI',passed=rep.shape==(2000,len(a)) and np.allclose(ci[0],a.ci_lower) and np.allclose(ci[1],a.ci_upper)))
 for db in ['mimic','amsterdam']:
  z=pd.read_parquet(RES/f'pairs_with_site_{db}.parquet');orig=pd.read_parquet(ROOT/f'phase1/_restricted/pairs_{db}.parquet');assert z.cid.is_unique
  cv=coverage[coverage.database.eq(db)].iloc[0];known=~z.site.isin(['unknown','mixed']);checks.append(dict(database=db,check='coverage_and_pair_count',passed=len(z)==len(orig)==cv.total_pairs and int(known.sum())==cv.assigned_pairs and z.loc[known,'patient'].nunique()==cv.assigned_patients and z.patient.nunique()==cv.total_patients))
  checks.append(dict(database=db,check='pair_columns_unchanged',passed=z.set_index('cid')[orig.columns.difference(['cid'])].sort_index().equals(orig.set_index('cid')[orig.columns.difference(['cid'])].sort_index())))
  dist=distribution[distribution.database.eq(db)];checks.append(dict(database=db,check='site_partition_sums_to_original',passed=int(dist.pairs.sum())==len(orig)))
  a=S3[S3.database.eq(db)].reset_index(drop=True);rep=np.load(RES/f'site_bootstrap_{db}.npz')['estimates'];ci=np.nanquantile(rep,[.025,.975],axis=0);checks.append(dict(database=db,check='S3_saved_replica_CI',passed=rep.shape==(2000,len(a)) and np.allclose(ci[0],a.ci_lower,equal_nan=True) and np.allclose(ci[1],a.ci_upper,equal_nan=True)))
  fi=a.index[a.site.eq('femoral')&a.metric.eq('missed_hypotension')&a.stratum.eq('all')][0];ri=a.index[a.site.eq('radial')&a.metric.eq('missed_hypotension')&a.stratum.eq('all')][0];ci=np.nanquantile(rep[:,fi]-rep[:,ri],[.025,.975]);e=compare[compare.database.eq(db)].iloc[0];checks.append(dict(database=db,check='joint_site_difference_CI',passed=np.allclose(ci,[e.ci_lower,e.ci_upper]) and np.isclose(a.loc[fi,'estimate']-a.loc[ri,'estimate'],e.estimate)))
  fem=z[z.site.eq('femoral')];sparse=len(fem)<1000 or fem.patient.nunique()<100;checks.append(dict(database=db,check='sparse_database_overall_only',passed=not sparse or a.stratum.eq('all').all()))
  checks.append(dict(database=db,check='sparse_femoral_precision_label',passed=not sparse or a.loc[a.site.eq('femoral'),'precision'].eq('insufficient_femoral_sample').all()))
  pt=pd.read_csv(ROOT/f'phase2/_restricted/table1_patients_{db}.csv')
  for row in a.itertuples():
   mask=z.site.eq(row.site)
   if row.stratum=='<=0.25':mask&=z.nee.le(.25)
   elif row.stratum=='>0.25-1.0':mask&=z.nee.gt(.25)&z.nee.le(1)
   elif row.stratum=='>1.0':mask&=z.nee.gt(1)
   if row.metric=='missed_hypotension':eligible=mask&z.art_map.lt(65);num=(eligible&z.cuff_map.ge(65)).sum()
   elif row.metric=='deep_missed_hypotension':eligible=mask&z.art_map.lt(55);num=(eligible&z.cuff_map.ge(65)).sum()
   elif row.metric=='art_lt65_given_cuff_ge65':eligible=mask&z.cuff_map.ge(65);num=(eligible&z.art_map.lt(65)).sum()
   else:eligible=mask&z.art_map.lt(55);num=(z.loc[eligible,'cuff_map']-z.loc[eligible,'art_map']).sum()
   passed=int(eligible.sum())==row.denominator and abs(num-row.numerator)<1e-8 and int(mask.sum())==row.layer_pairs and z.loc[mask,'patient'].nunique()==row.layer_patients and int((mask&z.art_map.lt(65)).sum())==row.layer_low_arterial_pairs
   checks.append(dict(database=db,check='S3_counts_'+row.site+'_'+row.metric+'_'+row.stratum,passed=passed))
  for row in baseline[baseline.database.eq(db)&baseline.variable.ne('age_group')].itertuples():
   patients=set(z.loc[z.site.eq(row.site),'patient']);x=pd.to_numeric(pt.loc[pt.patient.isin(patients),row.variable],errors='coerce');checks.append(dict(database=db,check='baseline_missingness_'+row.site+'_'+row.variable,passed=len(x)==row.patients and x.isna().sum()==row.missing and x.notna().sum()==row.available and np.isclose(x.median(),row.median,equal_nan=True) and np.isclose(x.quantile(.25),row.q1,equal_nan=True) and np.isclose(x.quantile(.75),row.q3,equal_nan=True)))
   quality.append(dict(database=db,site=row.site,variable=row.variable,patients=len(x),missing=int(x.isna().sum()),available=int(x.notna().sum()),missing_fraction=x.isna().mean(),minimum=x.min(),maximum=x.max(),median=x.median(),q1=x.quantile(.25),q3=x.quantile(.75),known_events=int(x.eq(1).sum()) if row.variable=='hospital_death' else np.nan))
  overlapping=len(set(z.loc[z.site.eq('radial'),'patient'])&set(z.loc[z.site.eq('femoral'),'patient']));checks.append(dict(database=db,check='overlapping_site_patients_preserved',passed=baseline.loc[baseline.database.eq(db),'radial_femoral_overlapping_patients'].eq(overlapping).all()))
  if db=='amsterdam':
   for site,g in baseline[baseline.database.eq(db)&baseline.variable.eq('age_group')].groupby('site'):
    ids=set(z.loc[z.site.eq(site),'patient']);counts=pt.loc[pt.patient.isin(ids),'age_group'].fillna('unknown').value_counts();reported=g.set_index('category')['n'];checks.append(dict(database=db,check='age_group_categories_'+site,passed=counts.to_dict()==reported.to_dict()))
 # Reused full-range REML numbers must agree exactly with delivered original numbers.
 old=pd.read_csv(ROOT/'phase1/phase1_secondary_bland_altman.csv');new=pd.read_csv(OUT/'phase4_S4_existing_overall_BA.csv')
 for r in new.itertuples():
  e=old[old.database.eq(r.database)&old.stratification.eq('overall')&old.stratum.eq('all')].iloc[0];checks.append(dict(database=r.database,check='old_REML_BA_reused_exactly',passed=all(getattr(r,k)==e[k] for k in ['bias','bias_ci_lower','bias_ci_upper','loa_lower','loa_upper','n_pairs','n_patients'])))
 outcsv(rows,'phase4_S4_independent_point_checks.csv');outcsv(checks,'phase4_independent_checks.csv');outcsv(quality,'phase4_independent_baseline_quality.csv')
 return rows,checks
def main():
 scanned=sql_sources();summary,cross=source_checks();s4,checks=aggregate_checks()
 (RES/'audit_work_receipt.json').write_text(json.dumps(dict(completed_at=datetime.datetime.now().astimezone().isoformat(),source_summary=summary,amsterdam_original_rows_scanned=scanned,all_aggregate_checks_pass=all(r['passed'] for r in checks),S4_points_match=all(r['passed'] for r in s4)),indent=2),encoding='utf-8')
 print(json.dumps(dict(source_summary=summary,abp_cross_tab=cross,aggregate_checks=len(checks),failed=[x for x in checks if not x['passed']],S4_cells=len(s4),S4_failures=sum(not x['passed'] for x in s4)),ensure_ascii=False))
if __name__=='__main__':main()
