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
import sys,json,hashlib,datetime
R=Path('@CUFF_WORKSPACE@');O=R/'phase3/icm_revision4';A=O/'analysis'
import numpy as np,pandas as pd
from scipy.sparse import csr_matrix
P=json.loads((A/'icm_rev4_params.yaml').read_text(encoding='utf-8'));checks=[];a1=[];a2=[];counts=[];flow=[];interpret=[]
def check(db,name,ok,detail='',critical=True):
 checks.append(dict(database=db,check=name,passed=bool(ok),detail=str(detail),critical=critical))
 if critical and not ok:raise AssertionError((db,name,detail))
def interval(v):
 v=v[np.isfinite(v)];return (*np.quantile(v,[.025,.975]),len(v)) if len(v) else (np.nan,np.nan,0)
for db in P['databases']:
 d0=pd.read_parquet(R/'phase1/_restricted'/f'pairs_{db}.parquet');d=d0[d0.nee.notna()].copy()
 check(db,'unique_cuff_anchors',d0.cid.is_unique)
 check(db,'known_dose_stratum',set(d.dose_stratum.astype(str))==set(P['dose_strata']))
 check(db,'complete_MAP',d[['cuff_map','art_map']].notna().all().all())
 for lab,q in [('all',d0),('known',d),('unknown',d0[d0.nee.isna()])]:flow.append(dict(database=db,population=lab,pairs=len(q),patients=q.pid.nunique(),arterial_lt65_pairs=int((q.art_map<65).sum()),arterial_lt65_prevalence=(q.art_map<65).mean(),dose_missing_fraction=q.nee.isna().mean(),cuff_missing=int(q.cuff_map.isna().sum()),arterial_missing=int(q.art_map.isna().sum())))
 ids,patients=pd.factorize(d.pid,sort=True);n=len(patients);c=d.cuff_map.to_numpy();art=d.art_map.to_numpy();dose=d.nee.to_numpy();groups=[];masks=[]
 for x in P['dose_strata']:groups.append(('five_strata',x));masks.append(d.dose_stratum.astype(str).eq(x).to_numpy())
 for x,m in [('<=0.25',dose<=.25),('>0.25-1.0',(dose>.25)&(dose<=1)),('>1.0',dose>1)]:groups.append(('three_groups',x));masks.append(m)
 check(db,'five_strata_count_conservation',sum(m.sum() for m in masks[:5])==len(d))
 check(db,'three_groups_count_conservation',sum(m.sum() for m in masks[5:])==len(d))
 cols=[];meta=[]
 for (scheme,st),m in zip(groups,masks):
  defs={'TP':m&(c<65)&(art<65),'FN':m&(c>=65)&(art<65),'FP':m&(c<65)&(art>=65),'TN':m&(c>=65)&(art>=65),'A1_num':m&(c>=70)&(art<65),'A1_den':m&(c>=70)}
  z={k:int(v.sum()) for k,v in defs.items()};z.update(database=db,grouping=scheme,stratum=st,all_pairs=int(m.sum()),all_patients=int(np.unique(ids[m]).size),all_normal_cuff=int((m&(c>=65)).sum()),cuff_65_69=int((m&(c>=65)&(c<70)).sum()),cuff_65_69_low=int((m&(c>=65)&(c<70)&(art<65)).sum()),positive_patients=int(np.unique(ids[m&(art<65)]).size),negative_patients=int(np.unique(ids[m&(art>=65)]).size),A1_den_patients=int(np.unique(ids[defs['A1_den']]).size))
  check(db,'2x2 conservation '+str((scheme,st)),z['TP']+z['FN']+z['FP']+z['TN']==z['all_pairs'])
  check(db,'A1 nested counts '+str((scheme,st)),z['A1_den']==z['all_normal_cuff']-z['cuff_65_69'] and z['A1_num']==z['FN']-z['cuff_65_69_low'])
  counts.append(z);meta.append(z)
  for k in ['TP','FN','FP','TN','A1_num','A1_den']:cols.append(np.bincount(ids[defs[k]],minlength=n))
 C=csr_matrix(np.column_stack(cols).astype(float));B=P['bootstrap']['replicates'];boot_counts=np.empty((B,len(groups),6));rng=np.random.default_rng(P['bootstrap']['seed'])
 for start in range(0,B,100):
  w=rng.multinomial(n,np.full(n,1/n),size=min(100,B-start));boot_counts[start:start+len(w)]=np.asarray(C.T@w.T).T.reshape(len(w),len(groups),6)
 with np.errstate(divide='ignore',invalid='ignore'):
  b1=boot_counts[:,:,4]/boot_counts[:,:,5]
  bminus=(boot_counts[:,:,1]/(boot_counts[:,:,0]+boot_counts[:,:,1]))/(boot_counts[:,:,3]/(boot_counts[:,:,2]+boot_counts[:,:,3]))
  bplus=(boot_counts[:,:,0]/(boot_counts[:,:,0]+boot_counts[:,:,1]))/(boot_counts[:,:,2]/(boot_counts[:,:,2]+boot_counts[:,:,3]))
 np.savez_compressed(A/'_restricted'/f'bootstrap_{db}.npz',counts=boot_counts,A1=b1,LR_negative=bminus,LR_positive=bplus,patient_counts=C.toarray())
 for j,z in enumerate(meta):
  base={k:z[k] for k in ['database','grouping','stratum','all_pairs','all_patients','positive_patients','negative_patients']}
  value=z['A1_num']/z['A1_den'];lo,hi,k=interval(b1[:,j]);risk65=z['FN']/z['all_normal_cuff']
  a1.append(dict(**base,numerator=z['A1_num'],denominator=z['A1_den'],denominator_patients=z['A1_den_patients'],estimate=value,ci_lower=lo,ci_upper=hi,bootstrap_valid=k,bootstrap_requested=B,boundary=z['A1_num'] in [0,z['A1_den']],precision='insufficient' if z['A1_den']<30 or z['A1_den_patients']<30 or hi-lo>.30 else 'adequate',cuff_ge65_risk=risk65))
  check(db,'A1 empirical monotonic check '+str((z['grouping'],z['stratum'])),value<=risk65,critical=False)
  for metric,v,b in [('LR_negative',(z['FN']/(z['TP']+z['FN']))/(z['TN']/(z['TN']+z['FP'])),bminus),('LR_positive',(z['TP']/(z['TP']+z['FN']))/(z['FP']/(z['TN']+z['FP'])),bplus)]:
   lo,hi,k=interval(b[:,j]);a2.append(dict(**base,metric=metric,estimate=v,ci_lower=lo,ci_upper=hi,bootstrap_valid=k,bootstrap_requested=B,undefined_replicates=int(np.isnan(b[:,j]).sum()),infinite_replicates=int(np.isinf(b[:,j]).sum()),TP=z['TP'],FN=z['FN'],FP=z['FP'],TN=z['TN'],missed_rate=z['FN']/(z['TP']+z['FN']),specificity=z['TN']/(z['TN']+z['FP']),precision='insufficient' if min(z['positive_patients'],z['negative_patients'])<30 else 'adequate'))
 # Local expected-count comparisons are omitted from this release.
 for metric in ['LR_negative','LR_positive']:
  rows=[x for x in a2 if x['database']==db and x['grouping']=='five_strata' and x['metric']==metric]
  overlap=max(x['ci_lower'] for x in rows)<=min(x['ci_upper'] for x in rows);low,high=rows[0],rows[-1]
  worse=high['estimate']>low['estimate'] and high['ci_lower']>low['ci_upper'] if metric=='LR_negative' else high['estimate']<low['estimate'] and high['ci_upper']<low['ci_lower']
  interpret.append(dict(database=db,metric=metric,all_five_CI_overlap=overlap,highest_clearly_worse_than_lowest=worse,highest_estimate=high['estimate'],lowest_estimate=low['estimate'],rule='descriptive_similarity_only' if overlap else 'omit_similarity_high_worse' if worse else 'neither_branch_report_direction_only'))
  if metric=='LR_negative':
   for x,expected in zip(rows,P['quality']['A2_expected_rounded'][db]):check(db,'user_LR_negative_'+x['stratum'],round(x['estimate'],2)==expected)
 print(json.dumps({'database':db,'known_pairs':len(d),'known_patients':n,'A1_high':a1[-1],'interpretation':interpret[-2:]},ensure_ascii=False,default=lambda x:x.item()),flush=True)
for data,name in [(a1,'A1_cuff_ge70.csv'),(a2,'A2_likelihood_ratios.csv'),(counts,'A_contingency_counts.csv'),(flow,'A_flow.csv'),(checks,'A_checks.csv'),(interpret,'A2_interpretation.csv')]:pd.DataFrame(data).to_csv(A/name,index=False,encoding='utf-8-sig')
pd.DataFrame(checks).to_csv(A/'A_checks.csv',index=False,encoding='utf-8-sig')
(A/'analysis_receipt.json').write_text(json.dumps(dict(completed_at=datetime.datetime.now().astimezone().isoformat(),A1_rows=len(a1),A2_rows=len(a2),checks=len(checks),all_passed=all(x['passed'] for x in checks),parameters_sha256=hashlib.sha256((A/'icm_rev4_params.yaml').read_bytes()).hexdigest()),indent=2),encoding='utf-8')
