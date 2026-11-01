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

from pathlib import Path
import os,sys,json,csv,hashlib,datetime,re
os.environ.setdefault('OPENBLAS_NUM_THREADS','2');os.environ.setdefault('OMP_NUM_THREADS','2')
ROOT=Path('@CUFF_WORKSPACE@');OUT=ROOT/'phase4';RES=OUT/'_restricted'

import numpy as np,pandas as pd,duckdb
P=json.loads((OUT/'phase4_params.yaml').read_text(encoding='utf-8'));PH=hashlib.sha256((OUT/'phase4_params.yaml').read_bytes()).hexdigest()
def save(rows,name):pd.DataFrame(rows).to_csv(OUT/name,index=False,encoding='utf-8-sig')
def load(db):
 d=pd.read_parquet(ROOT/f'phase1/_restricted/pairs_{db}.parquet');d.sid=d.sid.astype(str);assert d.cid.is_unique;return d
def intervals(db):
 s=pd.read_csv(RES/f'site_source_{db}.csv',keep_default_na=False)
 if db=='mimic':
  s['sid']=s.sid.astype(str);s['raw_site']=s.location.astype(str).str.strip();s=s[~s.statusdescription.eq('Rewritten')].copy()
  s['site']=s.raw_site.str.lower().map(lambda x:next((v for k,v in [('radial','radial'),('femoral','femoral'),('brachial','brachial'),('axilla','axillary'),('ulnar','ulnar'),('dorsalis','dorsalis_pedis')] if k in x),'unknown'))
  s['location_key']=np.where(s.site.eq('unknown'),'unknown',s.raw_site.str.lower());s['start']=pd.to_numeric(s.start,errors='coerce');s['stop']=pd.to_numeric(s.stop,errors='coerce')
 else:
  mp={9159:'radial',9160:'femoral',9186:'brachial',10619:'ulnar',10620:'dorsalis_pedis',13099:'axillary'}
  s=s[s.itemid.isin(mp)].copy();s['sid']=s.admissionid.astype(str);s['site']=s.itemid.map(mp);s['location_key']=s.site;s['raw_site']=s['item'];s['start']=pd.to_numeric(s.start,errors='coerce')/1000;s['stop']=pd.to_numeric(s.stop,errors='coerce')/1000
 valid=s.start.notna()&s.stop.notna()&s.stop.ge(s.start)
 save([dict(database=db,records=len(s),valid_intervals=int(valid.sum()),invalid_intervals=int((~valid).sum()),known_site_intervals=int((valid&~s.site.eq('unknown')).sum()))],f'site_interval_quality_{db}.csv')
 mapping=s[['itemid','raw_site','site','location_key']].drop_duplicates();save(mapping,f'site_mapping_{db}.csv')
 s=s[valid][['sid','start','stop','site','location_key']].drop_duplicates();s.to_parquet(RES/f'site_intervals_{db}.parquet',index=False);return s
def assign(d,i):
 c=duckdb.connect();a0=d[['cid','sid','t']].copy();i0=i.copy()
 for frame in [a0,i0]:
  for col in frame:
   if isinstance(frame[col].dtype,pd.StringDtype):frame[col]=frame[col].astype(object)
 c.register('anchors',a0);c.register('lines',i0)
 a=c.execute('SELECT a.cid,list_distinct(list(l.location_key)) locs,list_distinct(list(l.site)) sites FROM anchors a LEFT JOIN lines l ON a.sid=l.sid AND a.t>=l.start AND a.t<=l.stop GROUP BY a.cid').df();c.close()
 def cat(row):
  locs=[v for v in row.locs if v is not None];known=[v for v in locs if v!='unknown']
  if len(known)>1:return 'mixed','concurrent_distinct_locations'
  if 'unknown' in locs:return 'unknown','active_line_without_anatomical_site'
  if not known:return 'unknown','no_valid_interval_at_anchor'
  return next(v for v in row.sites if v is not None),'one_documented_site'
 labels=[cat(r) for r in a.itertuples()];a['site']=[x[0] for x in labels];a['site_reason']=[x[1] for x in labels]
 z=d.merge(a[['cid','site','site_reason']],on='cid',validate='one_to_one');assert len(z)==len(d);return z
def bootstrap(d,specs):
 pats=pd.Index(sorted(d.patient.unique()));codes=pats.get_indexer(d.patient);n=len(pats);vec=[];rows=[]
 for spec in specs:
  row,mask,event,kind=spec;mask=np.asarray(mask,bool);event=np.asarray(event);num=np.where(mask,event,0);cnt=np.bincount(codes,weights=mask,minlength=n);vec.extend([np.bincount(codes,weights=num,minlength=n),cnt]);nn=float(num.sum());den=int(mask.sum())
  row.update(numerator=nn if kind=='mean' else int(nn),denominator=den,denominator_patients=int((cnt>0).sum()),estimate=nn/den if den else np.nan,statistic='pair_weighted_mean_difference' if kind=='mean' else 'conditional_probability',parameter_sha256=PH);rows.append(row)
 a=np.array(vec).T;rng=np.random.default_rng(20261027);out=[]
 for b in range(0,2000,50):out.append(rng.multinomial(n,np.full(n,1/n),size=50)@a)
 cnt=np.concatenate(out);den=cnt[:,1::2];br=np.divide(cnt[:,0::2],den,out=np.full_like(den,np.nan),where=den>0)
 for j,row in enumerate(rows):
  good=np.isfinite(br[:,j]);lo,hi=np.quantile(br[good,j],[.025,.975]) if good.any() else (np.nan,np.nan);boundary=bool(row['statistic']=='conditional_probability' and row['denominator'] and row['numerator'] in [0,row['denominator']] and lo==hi)
  row.update(ci_lower=lo,ci_upper=hi,bootstrap_valid=int(good.sum()),ci_method='patient_cluster_percentile_2000',precision='not_estimable' if not row['denominator'] else 'boundary_not_informative' if boundary else 'insufficient' if row['denominator']<30 or (row['statistic']=='conditional_probability' and hi-lo>.30) else 'adequate')
  assert not row['denominator'] or (np.isfinite(row['estimate']) and lo<=hi)
  if row['statistic']=='conditional_probability':assert 0<=row['numerator']<=row['denominator']
 return rows,br
def main():
 coverage=[];distributions=[];site_results=[];bias_results=[];diffs=[];base=[];checks=[];decisions=[]
 for db in ['mimic','eicu','amsterdam','sicdb']:
  d=load(db);n=len(d);N=d.patient.nunique();checks.append(dict(database=db,check='unique_pair_id',passed=d.cid.is_unique))
  if db not in ['mimic','amsterdam']:
   coverage.append(dict(database=db,total_pairs=n,total_patients=N,assigned_pairs=0,assigned_patients=0,pair_coverage=0,patient_coverage=0,status='not_available',reason='No time-resolved anatomical catheter site; pulse/reflex sites are not catheter sites' if db=='eicu' else 'ArterialLine intervals present but the reviewed source Data payloads empty; dictionary and data_ref schema provide no anatomical site'))
   distributions.append(dict(database=db,site='unknown',pairs=n,patients=N,pair_fraction=1,patient_fraction=1))
  else:
   i=intervals(db);z=assign(d,i);z.to_parquet(RES/f'pairs_with_site_{db}.parquet',index=False);known=~z.site.isin(['unknown','mixed']);fem=z[z.site.eq('femoral')];sparse=len(fem)<1000 or fem.patient.nunique()<100
   coverage.append(dict(database=db,total_pairs=n,total_patients=N,assigned_pairs=int(known.sum()),assigned_patients=z.loc[known,'patient'].nunique(),pair_coverage=float(known.mean()),patient_coverage=z.loc[known,'patient'].nunique()/N,status='available',femoral_sparse=sparse))
   for s,g in z.groupby('site',sort=True):distributions.append(dict(database=db,site=s,pairs=len(g),patients=g.patient.nunique(),pair_fraction=len(g)/n,patient_fraction=g.patient.nunique()/N))
   others=[s for s,g in z[~z.site.isin(['radial','femoral','mixed','unknown'])].groupby('site') if len(g)>=1000 and g.patient.nunique()>=100];sites=['radial','femoral']+others
   art=z.art_map.to_numpy();cuff=z.cuff_map.to_numpy();lo=art<65;normal=cuff>=65;delta=cuff-art
   masks={'all':np.ones(n,bool),'<=0.25':z.nee.le(.25).to_numpy(),'>0.25-1.0':(z.nee.gt(.25)&z.nee.le(1)).to_numpy(),'>1.0':z.nee.gt(1).to_numpy()};spec=[]
   for s in sites:
    for metric,eligible,event,labs,kind in [('missed_hypotension',lo,normal,['all','<=0.25','>0.25-1.0','>1.0'],'rate'),('deep_missed_hypotension',art<55,normal,['all','>1.0'],'rate'),('art_lt65_given_cuff_ge65',normal,lo,['all','>1.0'],'rate'),('mean_difference_art_lt55',art<55,delta,['all'],'mean')]:
     for lab in (['all'] if sparse else labs):
      lay=z.site.eq(s).to_numpy()&masks[lab];row=dict(database=db,site=s,metric=metric,stratum=lab,layer_pairs=int(lay.sum()),layer_patients=z.loc[lay,'patient'].nunique(),layer_low_arterial_pairs=int((lay&lo).sum()),layer_low_arterial_prevalence=float((lay&lo).sum()/lay.sum()) if lay.any() else np.nan,femoral_sparse=sparse)
      spec.append((row,lay&eligible,event,kind))
   rr,br=bootstrap(z,spec)
   # The prespecified femoral sample threshold is separate from CI-width precision.
   for row in rr:
    if sparse and row['site']=='femoral':row['precision']='insufficient_femoral_sample'
   site_results.extend(rr);np.savez_compressed(RES/f'site_bootstrap_{db}.npz',estimates=br)
   f=next(j for j,r in enumerate(rr) if r['site']=='femoral' and r['metric']=='missed_hypotension' and r['stratum']=='all');r=next(j for j,x in enumerate(rr) if x['site']=='radial' and x['metric']=='missed_hypotension' and x['stratum']=='all');bb=br[:,f]-br[:,r];good=np.isfinite(bb);ci=np.quantile(bb[good],[.025,.975]) if good.any() else [np.nan,np.nan];effect=rr[f]['estimate']-rr[r]['estimate']
   rule='insufficient_precision_direction_only' if sparse else 'femoral_lower_over_10pp' if effect<-.10 else 'femoral_not_lower_over_10pp'
   diffs.append(dict(database=db,metric='femoral_minus_radial_missed_probability',estimate=effect,ci_lower=ci[0],ci_upper=ci[1],bootstrap_valid=int(good.sum()),femoral_pairs=len(fem),femoral_patients=fem.patient.nunique(),interpretation=rule))
   pt=pd.read_csv(ROOT/f'phase2/_restricted/table1_patients_{db}.csv');assert pt.patient.is_unique
   for s in sites:
    pop=pt[pt.patient.isin(z.loc[z.site.eq(s),'patient'])];both=len(set(z.loc[z.site.eq('radial'),'patient'])&set(z.loc[z.site.eq('femoral'),'patient']))
    for v in ['age','max_nee_observed','hospital_death']:
     if v=='age' and db=='amsterdam':
      for label,num in pop.age_group.fillna('unknown').value_counts().items():base.append(dict(database=db,site=s,variable='age_group',category=label,patients=len(pop),n=int(num),radial_femoral_overlapping_patients=both))
      continue
     x=pd.to_numeric(pop[v],errors='coerce');row=dict(database=db,site=s,variable=v,patients=len(pop),available=int(x.notna().sum()),missing=int(x.isna().sum()),missing_fraction=float(x.isna().mean()) if len(pop) else np.nan,median=x.median(),q1=x.quantile(.25),q3=x.quantile(.75),radial_femoral_overlapping_patients=both)
     if v=='hospital_death':row.update(events=int(x.eq(1).sum()),proportion=float(x.eq(1).sum()/x.notna().sum()) if x.notna().any() else np.nan)
     base.append(row)
   audit=pd.concat([g.sample(min(6,len(g)),random_state=20261027) for _,g in z.groupby('site')]).drop_duplicates('cid')
   if len(audit)>30:audit=audit.sample(30,random_state=20261027)
   if len(audit)<30:audit=pd.concat([audit,z[~z.cid.isin(audit.cid)].sample(min(30-len(audit),len(z)-len(audit)),random_state=20261027)])
   audit[['cid','sid','t','site','site_reason','art_map','cuff_map','nee']].to_csv(RES/f'site_audit_{db}.csv',index=False)
   checks.append(dict(database=db,check='site_pair_partition',passed=sum(x['pairs'] for x in distributions if x['database']==db)==n))
  if db in ['mimic','eicu']:
   art=d.art_map.to_numpy();delta=d.cuff_map.to_numpy()-art;spec=[]
   for axis,vals in [('arterial_MAP',art),('mean_of_cuff_and_arterial_MAP',(d.cuff_map.to_numpy()+art)/2)]:
    bins=np.digitize(vals,[55,65,75,85,95],right=False)
    for k,label in enumerate(P['S4']['labels']):spec.append((dict(database=db,axis=axis,pressure_range=label,n_pairs=int((bins==k).sum())),bins==k,delta,'mean'))
    checks.append(dict(database=db,check='S4_partition_'+axis,passed=sum(int((bins==k).sum()) for k in range(6))==n))
   rr,bb=bootstrap(d,spec);bias_results.extend(rr);np.savez_compressed(RES/f'bias_bootstrap_{db}.npz',estimates=bb)
   for row in rr:assert row['denominator']==row['n_pairs']
   frame=pd.DataFrame(rr);a=frame[frame.axis.eq('arterial_MAP')].set_index('pressure_range');m=frame[frame.axis.eq('mean_of_cuff_and_arterial_MAP')].set_index('pressure_range');low=['<55','55-64'];high=['85-94','>=95'];pattern_a=(a.loc[low,'estimate']>0).all() and (a.loc[high,'estimate']<0).all();pattern_m=(m.loc[low,'estimate']>0).all() and (m.loc[high,'estimate']<0).all();atten=bool((m.loc[low,'estimate']<a.loc[low,'estimate']).any())
   decisions.append(dict(database=db,arterial_axis_reversal=bool(pattern_a),mean_axis_reversal=bool(pattern_m),mean_axis_low_attenuation=atten,interpretation='descriptive_pressure_pattern_with_attenuation' if pattern_a and pattern_m and atten else 'pressure_dependent_direction_on_both_axes' if pattern_a and pattern_m else 'only_arterial_pattern_descriptive' if pattern_a else 'no_prespecified_complete_pattern',parameter_sha256=PH))
  print(db,'supplementary calculations completed',flush=True)
 assert all(r['passed'] for r in checks)
 save(coverage,'phase4_site_availability.csv');save(distributions,'phase4_site_distribution.csv');save(site_results,'phase4_S3_site.csv');save(diffs,'phase4_S3_comparison.csv');save(base,'phase4_S3_baseline.csv');save(bias_results,'phase4_S4_bias_by_range.csv');save(decisions,'phase4_S4_interpretation.csv');save(checks,'phase4_invariants.csv')
 ba=pd.read_csv(ROOT/'phase1/phase1_secondary_bland_altman.csv');save(ba[ba.database.isin(['mimic','eicu'])&ba.stratification.eq('overall')].assign(source='phase1/phase1_secondary_bland_altman.csv; reused unchanged'),'phase4_S4_existing_overall_BA.csv')
 (OUT/'phase4_analysis_receipt.json').write_text(json.dumps(dict(completed_at=datetime.datetime.now().astimezone().isoformat(),parameter_sha256=PH,source_audit='pending',s3_rows=len(site_results),s4_rows=len(bias_results)),indent=2),encoding='utf-8')
 print('ANALYSIS_COMPLETE_SOURCE_AUDIT_PENDING',flush=True)
if __name__=='__main__':main()
