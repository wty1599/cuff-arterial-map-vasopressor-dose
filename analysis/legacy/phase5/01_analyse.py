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
import os,sys,json,hashlib,datetime,time
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
ROOT=Path('@CUFF_WORKSPACE@');A=ROOT/'phase5';RES=A/'_restricted'

import numpy as np,pandas as pd
from scipy import sparse
from scipy.special import expit
P=json.loads((A/'phase5_params.yaml').read_text());PH=hashlib.sha256((A/'phase5_params.yaml').read_bytes()).hexdigest()
LABELS=P['dose_strata'];METRICS=['missed_hypotension','art_lt65_given_cuff_ge65']

def wquantile(x,n):
 h=(n.sum()-1)*np.array(P['primary']['knot_quantiles']);lo=np.floor(h).astype(int);hi=np.ceil(h).astype(int);cum=np.cumsum(n)
 return x[np.searchsorted(cum,lo,side='right')]*(1-h+lo)+x[np.searchsorted(cum,hi,side='right')]*(h-lo)

def fit_spline(x,n,y):
 k=wquantile(x,n)
 if len(np.unique(k))!=5:raise RuntimeError('Duplicate spline knots')
 t1=np.maximum(x-k[-2],0)**3;t2=np.maximum(x-k[-1],0)**3
 basis=[np.ones(len(x)),x]
 for kj in k[:-2]:basis.append((np.maximum(x-kj,0)**3-t1*(k[-1]-kj)/(k[-1]-k[-2])+t2*(k[-2]-kj)/(k[-1]-k[-2]))/(k[-1]-k[0])**2)
 X=np.column_stack(basis);mu=np.average(X[:,1:],weights=n,axis=0);sd=np.sqrt(np.average((X[:,1:]-mu)**2,weights=n,axis=0));X[:,1:]=(X[:,1:]-mu)/sd
 b=np.zeros(5);b[0]=np.log(y.sum()/(n.sum()-y.sum()))
 for it in range(100):
  eta=X@b;pr=expit(eta);score=X.T@(y-n*pr)
  if max(abs(score))/n.sum()<1e-10:return pr,k,it+1,max(abs(score))/n.sum()
  H=(X.T*(n*pr*(1-pr)))@X
  step=np.linalg.solve(H,score);old=np.sum(y*eta-n*np.logaddexp(0,eta));rate=1.
  for back in range(30):
   trial=b+rate*step;et=X@trial;ll=np.sum(y*et-n*np.logaddexp(0,et))
   if np.isfinite(ll) and ll>=old-1e-8:break
   rate*=.5
  else:raise RuntimeError('Line search failed')
  b=trial
  if max(abs(rate*step))<1e-8:return expit(X@b),k,it+1,max(abs(X.T@(y-n*expit(X@b))))/n.sum()
 raise RuntimeError('Model did not converge in 100 iterations')

def empirical(x,n,y):
 bins=np.where(x<40,0,np.where(x>120,82,np.where(x==120,81,np.floor(x).astype(int)-39)))
 den=np.bincount(bins,weights=n,minlength=83);num=np.bincount(bins,weights=y,minlength=83)
 p=np.divide(num,den,out=np.zeros_like(num),where=den>0)
 return p[bins]

def calculate(x,n,y,p):
 low=x<65;num=y[:,low].sum(axis=1);nlow=n[:,low].sum(axis=1);yn=y.sum(axis=1)
 obs=np.column_stack([num/nlow,num/yn]);en=(n[:,low]*p[low]).sum(axis=1)
 exp=np.column_stack([en/nlow,en/(n*p).sum(axis=1)])
 od=obs[-1,1]-obs[0,1];ed=exp[-1,1]-exp[0,1];ratio=ed/od if abs(od)>=1e-12 else np.nan
 return np.stack([obs,exp,obs-exp],axis=-1),np.array([od,ed,ratio])

def interval(v):
 good=np.isfinite(v);q=np.quantile(v[good],[.025,.975]) if good.any() else [np.nan,np.nan]
 return float(q[0]),float(q[1]),int(good.sum())

allrows=[];quality=[];model=[];interpret=[]
started=datetime.datetime.now().astimezone().isoformat()
for db in ['mimic','eicu']:
 d0=pd.read_parquet(ROOT/f'phase1/_restricted/pairs_{db}.parquet');assert d0.cid.is_unique
 d=d0[d0.nee.notna()].copy()
 assert d[['art_map','cuff_map','patient','nee']].notna().all().all()
 d['dose']=np.searchsorted([.1,.25,.5,1.0],d.nee.to_numpy(),side='left')
 pats,pc=np.unique(d.patient.to_numpy(),return_inverse=True);x,xc=np.unique(d.art_map.to_numpy(),return_inverse=True);g=d.dose.to_numpy();y0=d.cuff_map.ge(65).to_numpy().astype(float)
 N=len(pats);K=len(x);col=g*K+xc
 ns=sparse.csr_matrix((np.ones(len(d)),(pc,col)),shape=(N,5*K));ys=sparse.csr_matrix((y0,(pc,col)),shape=(N,5*K))
 n=np.asarray(ns.sum(axis=0)).reshape(5,K);y=np.asarray(ys.sum(axis=0)).reshape(5,K)
 nall=n.sum(axis=0);yall=y.sum(axis=0)
 method='restricted_cubic_spline_5_knots';fail_reason=''
 try:prob,knots,it,score=fit_spline(x,nall,yall)
 except (RuntimeError,np.linalg.LinAlgError) as e:
  method='empirical_1mmHg_fallback';fail_reason=str(e);prob=empirical(x,nall,yall);knots=np.array([]);it=0;score=np.nan
 point,ratio=calculate(x,n,y,prob);ep=empirical(x,nall,yall);sens,sensratio=calculate(x,n,y,ep)
 # Aggregated sufficient counts are retained privately for independent verification.
 pd.DataFrame([dict(database=db,arterial_MAP=a,stratum=LABELS[j],pairs=n[j,i],cuff_ge65=y[j,i]) for j in range(5) for i,a in enumerate(x) if n[j,i]>0]).to_csv(RES/f'point_counts_{db}.csv',index=False)
 br=np.full((2000,5,2,3),np.nan);bs=np.full_like(br,np.nan);rr=np.full((2000,3),np.nan);sr=np.full_like(rr,np.nan);fits=[];rng=np.random.default_rng(20261027);t=time.time()
 for b in range(2000):
  w=rng.multinomial(N,np.full(N,1/N));nb=np.asarray(w@ns).reshape(5,K);yb=np.asarray(w@ys).reshape(5,K);na=nb.sum(axis=0);ya=yb.sum(axis=0)
  try:
   if method=='restricted_cubic_spline_5_knots':pb,kb,ib,scb=fit_spline(x,na,ya)
   else:pb=empirical(x,na,ya);kb=[];ib=0;scb=np.nan
   br[b],rr[b]=calculate(x,nb,yb,pb);fits.append({'replicate':b+1,'converged':True,'iterations':ib,'scaled_score':scb,'knots':list(kb)})
  except (RuntimeError,np.linalg.LinAlgError) as e:fits.append({'replicate':b+1,'converged':False,'reason':str(e)})
  bs[b],sr[b]=calculate(x,nb,yb,empirical(x,na,ya))
  if (b+1)%500==0:print(json.dumps({'database':db,'replicates':b+1,'elapsed_seconds':round(time.time()-t,1)}),flush=True)
 np.savez_compressed(RES/f'bootstrap_{db}.npz',primary=br,ratio=rr,empirical=bs,empirical_ratio=sr)
 (RES/f'fit_log_{db}.json').write_text(json.dumps(fits),encoding='utf-8')
 for typ,pm,bm,rat,bra,idxs in [('primary',point,br,ratio,rr,range(5)),('empirical_1mmHg',sens,bs,sensratio,sr,[4])]:
  for j in idxs:
   mask=g==j
   for m,metric in enumerate(METRICS):
    row=dict(database=db,analysis=typ,method=method if typ=='primary' else 'empirical_1mmHg',metric=metric,stratum=LABELS[j],n_pairs=int(n[j].sum()),n_patients=len(np.unique(pc[mask])),arterial_hypotension_prevalence=n[j,x<65].sum()/n[j].sum(),observed_numerator=int(y[j,x<65].sum()),observed_denominator=int(n[j,x<65].sum() if m==0 else y[j].sum()),precision='insufficient' if n[j].sum()<500 else 'adequate',parameter_sha256=PH)
    for z,key in enumerate(['observed','expected','difference']):
     lo,hi,valid=interval(bm[:,j,m,z]);row[key]=pm[j,m,z];row[key+'_ci_lower']=lo;row[key+'_ci_upper']=hi;row[key+'_bootstrap_valid']=valid
    assert 0<=row['observed']<=1 and 0<=row['expected']<=1
    allrows.append(row)
  row=dict(database=db,analysis=typ,method=method if typ=='primary' else 'empirical_1mmHg',metric='explained_ratio',stratum='>1.0_vs_<=0.1',parameter_sha256=PH)
  for i,key in enumerate(['observed_increase','expected_increase','explained_ratio']):
   lo,hi,valid=interval(bra[:,i]);row[key]=rat[i];row[key+'_ci_lower']=lo;row[key+'_ci_upper']=hi;row[key+'_bootstrap_valid']=valid
  row['undefined_ratio_replicates']=int((~np.isfinite(bra[:,2])).sum());allrows.append(row)
 hirow=next(r for r in allrows if r['database']==db and r['analysis']=='primary' and r['metric']==METRICS[1] and r['stratum']=='>1.0')
 include=hirow['difference_ci_lower']<=0<=hirow['difference_ci_upper'];f=ratio[2]
 branch='dominant_compatible' if f>=.75 and include else 'partial' if .25<=f<.75 else 'low_explanation' if f<.25 else 'large_explanation_with_residual' if np.isfinite(f) else 'not_estimable'
 interpret.append(dict(database=db,branch=branch,explained_ratio=f,high_difference=hirow['difference'],high_difference_ci_lower=hirow['difference_ci_lower'],high_difference_ci_upper=hirow['difference_ci_upper'],expected_miss_high=point[4,0,1],expected_miss_low=point[0,0,1],expected_miss_nonincreasing=bool(np.all(np.diff(point[:,0,1])<=0))))
 model.append(dict(database=db,method=method,knots=list(knots),iterations=it,scaled_score=score,original_fit_fallback_reason=fail_reason,bootstrap_converged=sum(r['converged'] for r in fits),bootstrap_failed=sum(not r['converged'] for r in fits),patients=N,pairs=len(d),distinct_arterial_values=K))
 for j in range(5):quality.append(dict(database=db,stratum=LABELS[j],all_pairs=len(d0),known_pairs=len(d),known_fraction=len(d)/len(d0),stratum_pairs=int(n[j].sum()),stratum_patients=len(np.unique(pc[g==j])),arterial_low_pairs=int(n[j,x<65].sum()),cuff_normal_pairs=int(y[j].sum()),joint_low_normal_pairs=int(y[j,x<65].sum()),pair_fields_missing=int(d[['art_map','cuff_map','nee','patient']].isna().sum().sum()),pair_ids_unique=d.cid.is_unique))
 assert sum(r['stratum_pairs'] for r in quality if r['database']==db)==len(d)
pd.DataFrame(allrows).to_csv(A/'phase5_S5_counterfactual.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(quality).to_csv(A/'phase5_quality_checks.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(interpret).to_csv(A/'phase5_interpretation.csv',index=False,encoding='utf-8-sig')
(A/'phase5_model_receipt.json').write_text(json.dumps(model,indent=2),encoding='utf-8')
(A/'phase5_analysis_receipt.json').write_text(json.dumps({'started':started,'completed':datetime.datetime.now().astimezone().isoformat(),'parameter_sha256':PH,'seed':20261027,'replicates':2000,'numpy':np.__version__,'pandas':pd.__version__,'source_level_verification':'Unchanged phase1 pairs; reuse phase0-4 independent source checks, no new raw-source query claimed.','independent_computational_audit':'pending'},indent=2),encoding='utf-8')
print(json.dumps({'rows':len(allrows),'interpretation':interpret,'models':model},indent=2))
