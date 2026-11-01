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
R=Path('@CUFF_WORKSPACE@'); A=R/'phase6'; I=R/'phase1/_restricted'
import numpy as np,pandas as pd
from scipy.sparse import csr_matrix
P=json.loads((A/'phase6_params.yaml').read_text(encoding='utf-8'))
STRATA=P['additional_strata']+P['dose_strata']; WINDOWS=['5min','2min','same_minute']
METRICS=list(P['metrics']); checks=[]; rates=[]; effects=[]; dose_contrasts=[]; flows=[]; profiles=[]; transitions=[]
start=datetime.datetime.now().astimezone().isoformat()
def check(db,name,ok,detail=''):
    checks.append(dict(database=db,check=name,passed=bool(ok),detail=str(detail)))
    if not ok: raise AssertionError((db,name,detail))
def save(rows,name):pd.DataFrame(rows).to_csv(A/name,index=False,encoding='utf-8-sig')
def interval(x):
    x=x[np.isfinite(x)];return (*np.quantile(x,[.025,.975]),len(x)) if len(x) else (np.nan,np.nan,0)
def summaries(db,frame,d):
    people=d.drop_duplicates('pid')
    for scope,sub in [('patient',people),('pair',d)]:
        for col in ['age','weight','height','bmi','nee','af_model']:
            x=pd.to_numeric(sub[col],errors='coerce');z=x.dropna()
            profiles.append(dict(database=db,frame=frame,unit=scope,variable=col,n=len(x),missing=int(x.isna().sum()),missing_fraction=x.isna().mean(),median=z.median(),q25=z.quantile(.25),q75=z.quantile(.75)))
    for s in STRATA:
        m=np.ones(len(d),bool) if s=='all' else d.dose_stratum.ne('unknown').to_numpy() if s=='known' else d.dose_stratum.eq(s).to_numpy()
        q=d.loc[m]
        flows.append(dict(database=db,frame=frame,stratum=s,pairs=len(q),patients=q.pid.nunique(),known_dose_pairs=int(q.nee.notna().sum())))

for db in P['databases']:
    d=pd.read_parquet(I/f'pairs_{db}.parquet').rename(columns={'art_map':'art_5min'})
    check(db,'unique_original_cuff_ids',d.cid.is_unique)
    d['dose_stratum']=d.dose_stratum.astype(str)
    for stem,w in [('window2','2min'),('same_minute','same_minute')]:
        s=pd.read_parquet(I/f'{stem}_{db}.parquet')
        check(db,f'{w}_unique_cuff_ids',s.cid.is_unique)
        d=d.merge(s[['cid','art_map']].rename(columns={'art_map':f'art_{w}'}),on='cid',how='left',validate='one_to_one')
    eligible={w:d[f'art_{w}'].notna().to_numpy() for w in WINDOWS}
    common=np.logical_and.reduce(list(eligible.values()))
    check(db,'same_minute_subset_of_2min',np.all(~eligible['same_minute']|eligible['2min']))
    check(db,'original_5min_complete',eligible['5min'].all())
    check(db,'common_equals_same_minute',np.array_equal(common,eligible['same_minute']))
    d['common_anchor']=common
    d.to_parquet(A/'_restricted'/f'aligned_{db}.parquet',index=False)
    for f,m in [('native_5min',eligible['5min']),('native_2min',eligible['2min']),('common',common),('excluded_from_common',~common)]:summaries(db,f,d.loc[m])
    pid,unique=pd.factorize(d.pid,sort=True);n=len(unique);c=d.cuff_map.to_numpy(float)
    cols=[];meta=[];lookup={}
    for frame in ['native','common']:
        for w in WINDOWS:
            a=d[f'art_{w}'].to_numpy(float);valid=eligible[w] if frame=='native' else common
            for s in STRATA:
                sm=np.ones(len(d),bool) if s=='all' else d.dose_stratum.ne('unknown').to_numpy() if s=='known' else d.dose_stratum.eq(s).to_numpy()
                m=valid&sm
                defs={
                  'missed':(m&(a<65)&(c>=65),m&(a<65)),
                  'deep_missed':(m&(a<55)&(c>=65),m&(a<55)),
                  'cuff_conditioned':(m&(a<65)&(c>=65),m&(c>=65)),
                  'cuff_65_69':(m&(a<65)&(c>=65)&(c<70),m&(c>=65)&(c<70)),
                  'cuff_65_69_deep':(m&(a<60)&(c>=65)&(c<70),m&(c>=65)&(c<70)),
                  'hypotension_prevalence':(m&(a<65),m)}
                for metric,(num,den) in defs.items():
                    nc=np.bincount(pid[num],minlength=n);dc=np.bincount(pid[den],minlength=n)
                    check(db,f'count_{frame}_{w}_{s}_{metric}',np.all(nc<=dc))
                    index=len(meta);lookup[(frame,w,s,metric)]=index
                    N,D=int(num.sum()),int(den.sum())
                    meta.append(dict(database=db,frame=frame,window=w,stratum=s,metric=metric,numerator=N,denominator=D,denominator_patients=int((dc>0).sum()),all_pairs=int(m.sum()),all_patients=int(np.unique(pid[m]).size),estimate=N/D if D else np.nan,bootstrap_column=index))
                    cols.extend([nc,dc])
    count=csr_matrix(np.column_stack(cols).astype(float));point=np.array([x['estimate'] for x in meta])
    rng=np.random.default_rng(P['bootstrap']['seed']);B=P['bootstrap']['replicates'];boot=np.full((B,len(meta)),np.nan)
    for b in range(0,B,100):
        weights=rng.multinomial(n,np.full(n,1/n),size=min(100,B-b))
        sums=(count.T@weights.T).T
        with np.errstate(divide='ignore',invalid='ignore'):boot[b:b+len(weights)]=sums[:,::2]/sums[:,1::2]
    for i,row in enumerate(meta):
        lo,hi,k=interval(boot[:,i]);row.update(ci_lower=lo,ci_upper=hi,bootstrap_valid=k,boundary=bool(row['denominator'] and row['numerator'] in [0,row['denominator']]),precision='insufficient' if row['denominator']<30 or row['denominator_patients']<30 or (np.isfinite(hi-lo) and hi-lo>.30) else 'adequate')
    rates.extend(meta)
    np.savez_compressed(A/'_restricted'/f'joint_bootstrap_{db}.npz',estimates=boot,counts=count.toarray())
    (A/'_restricted'/f'bootstrap_columns_{db}.json').write_text(json.dumps(meta,allow_nan=True),encoding='utf-8')
    for w in ['2min','same_minute']:
        for s in STRATA:
            for metric in METRICS:
                ix=[lookup[(f,win,s,metric)] for f,win in [('native','5min'),('common','5min'),('common',w)]]
                baseline,selected,target=point[ix];b0,b1,b2=boot[:,ix].T
                pieces={'total':(target-baseline,b2-b0),'selection':(selected-baseline,b1-b0),'within_anchor_window':(target-selected,b2-b1)}
                if all(np.isfinite([baseline,selected,target])):
                    check(db,f'decomposition_{w}_{s}_{metric}',abs(pieces['total'][0]-pieces['selection'][0]-pieces['within_anchor_window'][0])<1e-12)
                for effect,(val,vec) in pieces.items():
                    lo,hi,k=interval(vec)
                    effects.append(dict(database=db,target_window=w,stratum=s,metric=metric,component=effect,estimate=val,ci_lower=lo,ci_upper=hi,bootstrap_valid=k,baseline_native5=baseline,common5=selected,common_target=target))
                finite=np.isfinite(b0)&np.isfinite(b1)&np.isfinite(b2)
                check(db,f'bootstrap_identity_{w}_{s}_{metric}',np.all(np.abs((b2[finite]-b0[finite])-(b1[finite]-b0[finite])-(b2[finite]-b1[finite]))<1e-12))
    for frame in ['native','common']:
        for w in WINDOWS:
            for metric in METRICS:
                loidx=lookup[(frame,w,STRATA[3],metric)];hiidx=lookup[(frame,w,STRATA[-1],metric)]
                lo,hi,k=interval(boot[:,hiidx]-boot[:,loidx])
                dose_contrasts.append(dict(database=db,frame=frame,window=w,metric=metric,contrast='>1.0 minus <=0.1',estimate=point[hiidx]-point[loidx],ci_lower=lo,ci_upper=hi,bootstrap_valid=k))
    # Counts are mutually exclusive at fixed cuff anchors; they expose threshold reclassification.
    for w in ['2min','same_minute']:
        for s in STRATA:
            sm=np.ones(len(d),bool) if s=='all' else d.dose_stratum.ne('unknown').to_numpy() if s=='known' else d.dose_stratum.eq(s).to_numpy()
            for threshold in [55,60,65]:
                a0=d.art_5min.to_numpy()<threshold;a1=d[f'art_{w}'].to_numpy()<threshold
                for condition,cm in [('all_cuff',np.ones(len(d),bool)),('cuff_ge65',c>=65),('cuff_65_69',(c>=65)&(c<70))]:
                    base=common&sm&cm;total=0
                    for old in [False,True]:
                        for new in [False,True]:
                            m=base&(a0==old)&(a1==new);N=int(m.sum());total+=N
                            transitions.append(dict(database=db,target_window=w,stratum=s,arterial_threshold=threshold,cuff_condition=condition,low_at_5min=old,low_at_target=new,pairs=N,patients=int(np.unique(pid[m]).size)))
                    check(db,f'transition_sum_{w}_{s}_{threshold}_{condition}',total==int(base.sum()))
    # Reproduce published native point estimates without replacing their original confidence intervals.
    for filename,w in [(f'phase1_primary_{db}.csv','5min'),('phase1_sensitivity_2min.csv','2min'),('phase1_sensitivity_same_minute.csv','same_minute')]:
        old=pd.read_csv(R/'phase1'/filename)
        old=old[(old.database==db)&(old.endpoint=='missed')]
        for row in old.itertuples():
            x=meta[lookup[('native',w,row.stratum,'missed')]]
            check(db,f'reproduce_{w}_{row.stratum}',x['numerator']==row.numerator and x['denominator']==row.denominator and abs(x['estimate']-row.estimate)<1e-12)
    print(json.dumps({'database':db,'original_pairs':len(d),'common_pairs':int(common.sum()),'patients':n,'rates':len(meta),'bootstrap':B}),flush=True)
save(rates,'phase6_window_rates.csv');save(effects,'phase6_decomposition.csv');save(dose_contrasts,'phase6_dose_contrasts.csv')
save(flows,'phase6_flow.csv');save(profiles,'phase6_profile.csv');save(transitions,'phase6_threshold_transitions.csv');save(checks,'phase6_quality_checks.csv')
(A/'analysis_receipt.json').write_text(json.dumps({'start':start,'end':datetime.datetime.now().astimezone().isoformat(),'params_sha256':hashlib.sha256((A/'phase6_params.yaml').read_bytes()).hexdigest(),'checks':len(checks),'all_passed':all(x['passed'] for x in checks),'rates':len(rates),'contrasts':len(effects),'protected_files_unchanged':None},indent=2),encoding='utf-8')
print('COMPLETE',len(rates),'rates',len(effects),'decomposition rows',len(checks),'checks',flush=True)
