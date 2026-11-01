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
from pathlib import Path
import sys,json,hashlib,datetime
ROOT=Path('@CUFF_WORKSPACE@');A=ROOT/'phase5';RES=A/'_restricted'
import numpy as np,pandas as pd
PROB_TOL=1e-7;RATIO_TOL=1e-6
(A/'phase5_independent_comparison_rules.json').write_text(json.dumps({'recorded_before_comparison':datetime.datetime.now().astimezone().isoformat(),'probability_absolute_tolerance':PROB_TOL,'ratio_absolute_tolerance':RATIO_TOL,'count_tolerance':0,'interval_check_tolerance':1e-12,'scope':'Independent R point estimates and one explicit patient bootstrap reconstruction per database; remaining saved replicates checked without refitting.'},indent=2),encoding='utf-8')
primary=pd.read_csv(A/'phase5_S5_counterfactual.csv');ind=pd.read_csv(A/'phase5_independent_R_points.csv');models=pd.read_csv(A/'phase5_independent_R_models.csv')
labels=['<=0.1','>0.1-0.25','>0.25-0.5','>0.5-1.0','>1.0'];metrics=['missed_hypotension','art_lt65_given_cuff_ge65'];checks=[];boot_receipt=[]
def check(db,scope,key,field,reference,audit,tol):
 reference=float(reference);audit=float(audit);difference=audit-reference
 checks.append(dict(database=db,scope=scope,key=key,field=field,production=reference,independent=audit,absolute_difference=abs(difference),tolerance=tol,passed=bool(abs(difference)<=tol or np.isnan(reference) and np.isnan(audit))))
for db in ['mimic','eicu']:
 data=np.load(RES/f'bootstrap_{db}.npz');logs=json.loads((RES/f'fit_log_{db}.json').read_text(encoding='utf-8'))
 for _,r in primary[primary.database.eq(db)].iterrows():
  matches=ind[ind.database.eq(db)&ind.replicate.eq(0)&ind.analysis.eq(r.analysis)&ind.metric.eq(r.metric)&ind.stratum.eq(r.stratum)]
  assert len(matches)==1;s=matches.iloc[0]
  fields=['observed_increase','expected_increase','explained_ratio'] if r.metric=='explained_ratio' else ['observed','expected','difference']
  key='|'.join([r.analysis,r.metric,r.stratum])
  for f in fields:check(db,'independent_R_point',key,f,r[f],s[f],RATIO_TOL if f=='explained_ratio' else PROB_TOL)
  array=data['ratio' if r.analysis=='primary' else 'empirical_ratio'] if r.metric=='explained_ratio' else data['primary' if r.analysis=='primary' else 'empirical'][:,labels.index(r.stratum),metrics.index(r.metric),:]
  sr=ind[ind.database.eq(db)&ind.replicate.eq(1)&ind.analysis.eq(r.analysis)&ind.metric.eq(r.metric)&ind.stratum.eq(r.stratum)].iloc[0]
  for i,f in enumerate(fields):
   check(db,'independent_R_first_patient_bootstrap',key,f,array[0,i],sr[f],RATIO_TOL if f=='explained_ratio' else PROB_TOL)
   finite=array[:,i][np.isfinite(array[:,i])];lo,hi=np.quantile(finite,[.025,.975])
   check(db,'saved_replicate_interval',key,f+'_ci_lower',r[f+'_ci_lower'],lo,1e-12)
   check(db,'saved_replicate_interval',key,f+'_ci_upper',r[f+'_ci_upper'],hi,1e-12)
   check(db,'saved_replicate_interval',key,f+'_bootstrap_valid',r[f+'_bootstrap_valid'],len(finite),0)
 for method,ratname in [('primary','ratio'),('empirical','empirical_ratio')]:
  arr=data[method];rat=data[ratname]
  check(db,'replicate_arithmetic',method,'max_abs_difference_identity',0,np.nanmax(abs(arr[:,:,:,2]-(arr[:,:,:,0]-arr[:,:,:,1]))),1e-12)
  check(db,'replicate_arithmetic',method,'max_abs_observed_increase_identity',0,np.nanmax(abs(rat[:,0]-(arr[:,-1,1,0]-arr[:,0,1,0]))),1e-12)
  check(db,'replicate_arithmetic',method,'max_abs_expected_increase_identity',0,np.nanmax(abs(rat[:,1]-(arr[:,-1,1,1]-arr[:,0,1,1]))),1e-12)
  valid=abs(rat[:,0])>=1e-12
  check(db,'replicate_arithmetic',method,'max_abs_unclipped_ratio_identity',0,np.nanmax(abs(rat[valid,2]-rat[valid,1]/rat[valid,0])),1e-12)
  boot_receipt.append(dict(database=db,method=method,replicates=len(rat),finite_ratio=int(np.isfinite(rat[:,2]).sum()),ratio_below_zero=int((rat[:,2]<0).sum()),ratio_above_one=int((rat[:,2]>1).sum()),ratio_min=float(np.nanmin(rat[:,2])),ratio_max=float(np.nanmax(rat[:,2])),observed_contrast_negative=int((rat[:,0]<0).sum())))
 assert len(logs)==2000 and [r['replicate'] for r in logs]==list(range(1,2001))
 assert all(r['converged'] for r in logs)
 boot_receipt.append(dict(database=db,method='fit_log',replicates=len(logs),unique_knot_vectors=len({tuple(r['knots']) for r in logs}),max_scaled_score=max(r['scaled_score'] for r in logs),min_iterations=min(r['iterations'] for r in logs),max_iterations=max(r['iterations'] for r in logs)))
 for replicate in [0,1]:
  k=models[models.database.eq(db)&models.replicate.eq(replicate)].iloc[0].knots
  independent=[float(t) for t in k.split(';')]
  expected=logs[0]['knots'] if replicate else next(x for x in json.loads((A/'phase5_model_receipt.json').read_text(encoding='utf-8')) if x['database']==db)['knots']
  for i,(x,y) in enumerate(zip(expected,independent)):check(db,'knots','replicate_'+str(replicate),'knot_'+str(i+1),x,y,1e-10)
 # This links an unchanged input to applicable prior source evidence; not a new source query.
 earlier=pd.read_csv(ROOT/'phase1/phase1_source_audit_summary.csv');pop=pd.read_csv(A/'phase5_independent_population.csv')
 assert earlier[earlier.database.eq(db)].iloc[0].production_file_sha256==pop[pop.database.eq(db)].iloc[0].original_pair_sha256
counts=pd.read_csv(A/'phase5_independent_dose_counts.csv');q=pd.read_csv(A/'phase5_quality_checks.csv')
for _,c in counts.iterrows():
 r=q[q.database.eq(c.database)&q.stratum.eq(c.stratum)].iloc[0]
 for target,source in [('stratum_pairs','pairs'),('stratum_patients','patients'),('arterial_low_pairs','low_arterial_pairs'),('cuff_normal_pairs','normal_cuff_pairs'),('joint_low_normal_pairs','joint_low_normal')]:
  check(c.database,'independent_count',c.stratum,target,r[target],c[source],0)
old_b=pd.read_csv(ROOT/'phase2/phase2_B_cuff_perspective.csv')
for db in ['mimic','eicu']:
 old_a=pd.read_csv(ROOT/f'phase1/phase1_primary_{db}.csv')
 for _,r in primary[primary.database.eq(db)&primary.analysis.eq('primary')&~primary.metric.eq('explained_ratio')].iterrows():
  if r.metric=='missed_hypotension':
   old=old_a[old_a.scenario.eq('primary_5min')&old_a.endpoint.eq('missed')&old_a.stratum.eq(r.stratum)]
  else:
   old=old_b[old_b.database.eq(db)&old_b.scenario.eq('primary_5min')&old_b.metric.eq('art_lt65_given_cuff_ge65')&old_b.stratum.eq(r.stratum)&old_b.category.eq('all_normal_cuff')]
  assert len(old)==1
  for source,target in [('estimate','observed'),('numerator','observed_numerator'),('denominator','observed_denominator')]:
   check(db,'preserved_prior_point',r.metric+'|'+r.stratum,target,old.iloc[0][source],r[target],1e-12 if source=='estimate' else 0)
protected=[]
for r in json.loads((A/'protected_inputs_before.json').read_text(encoding='utf-8')):
 p=Path(r['path']);actual=hashlib.sha256(p.read_bytes()).hexdigest()
 protected.append(dict(path=str(p),expected_sha256=r['sha256'],actual_sha256=actual,unchanged=actual==r['sha256']))
pd.DataFrame(protected).to_csv(A/'phase5_independent_protected_inputs.csv',index=False,encoding='utf-8-sig')
df=pd.DataFrame(checks);df.to_csv(A/'phase5_independent_comparisons.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(boot_receipt).to_csv(A/'phase5_independent_bootstrap_checks.csv',index=False,encoding='utf-8-sig')
print(json.dumps({'checks':len(df),'failed':int((~df.passed).sum()),'point_max_abs_diff':df[df.scope.eq('independent_R_point')].absolute_difference.max(),'one_bootstrap_max_abs_diff':df[df.scope.eq('independent_R_first_patient_bootstrap')].absolute_difference.max(),'all_passed':bool(df.passed.all())}))
