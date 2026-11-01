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

"""Independent patient multiplicity bootstrap; no main bootstrap arrays read."""
from pathlib import Path
import sys,json,datetime
import numpy as np
import pandas as pd
ROOT=Path('@CUFF_WORKSPACE@'); BASE=ROOT/'phase3/icm_revision4'; QA=BASE/'qa'
STRATA=[('five','<=0.1',0,.1),('five','>0.1-0.25',.1,.25),('five','>0.25-0.5',.25,.5),('five','>0.5-1.0',.5,1),('five','>1.0',1,np.inf),('three','<=0.25',0,.25),('three','>0.25-1.0',.25,1),('three','>1.0',1,np.inf)]
outputs=[]
for db in ['mimic','eicu']:
    d=pd.read_parquet(ROOT/f'phase1/_restricted/pairs_{db}.parquet')
    d=d[d.nee.notna()].copy()
    pi,patients=pd.factorize(d.pid.astype(str),sort=True); n=len(patients)
    counts=np.zeros((n,48),dtype=np.int64)
    for j,(_,_,lo,hi) in enumerate(STRATA):
        g=d.nee.gt(lo)&d.nee.le(hi); pos=d.art_map.lt(65); cp=d.cuff_map.lt(65)
        masks=[g&pos&cp,g&pos&~cp,g&~pos&cp,g&~pos&~cp,g&pos&d.cuff_map.ge(70),g&d.cuff_map.ge(70)]
        for k,m in enumerate(masks): counts[:,j*6+k]=np.bincount(pi,weights=m.astype(np.int64),minlength=n)
    rng=np.random.default_rng(20261027)
    bs=[]
    for start in range(0,2000,100):
        w=rng.multinomial(n,np.full(n,1/n),size=100)
        bs.append(w@counts)
    sums=np.concatenate(bs).reshape(2000,8,6)
    with np.errstate(divide='ignore',invalid='ignore'):
        neg=(sums[:,:,1]/(sums[:,:,0]+sums[:,:,1]))/(sums[:,:,3]/(sums[:,:,3]+sums[:,:,2]))
        pos=(sums[:,:,0]/(sums[:,:,0]+sums[:,:,1]))/(sums[:,:,2]/(sums[:,:,3]+sums[:,:,2]))
        a1=sums[:,:,4]/sums[:,:,5]
    for j,(grouping,stratum,lo,hi) in enumerate(STRATA):
        for metric,vals in [('A1',a1[:,j]),('LR_negative',neg[:,j]),('LR_positive',pos[:,j])]:
            finite=vals[np.isfinite(vals)];ci=np.quantile(finite,[.025,.975])
            outputs.append(dict(database=db,grouping=grouping,stratum=stratum,metric=metric,
                ci_lower=ci[0],ci_upper=ci[1],bootstrap_valid=len(finite),bootstrap_requested=2000,
                bootstrap_frame_patients=n,undefined_replicates=int(np.isnan(vals).sum()),infinite_replicates=int(np.isinf(vals).sum())))
    print(json.dumps({'database':db,'independent_bootstrap_complete':True,'frame_patients':n,'replicates':2000}))
out=pd.DataFrame(outputs);out.to_csv(QA/'audit_bootstrap_intervals.csv',index=False)
# Main aggregate CSV intervals are read only after all independent draws finish.
maina=pd.read_csv(BASE/'analysis/A1_cuff_ge70.csv');mainb=pd.read_csv(BASE/'analysis/A2_likelihood_ratios.csv')
checks=[]
for row in outputs:
    main=maina if row['metric']=='A1' else mainb[mainb.metric.eq(row['metric'])]
    main=main[(main.database==row['database'])&(main.stratum==row['stratum'])]
    if len(main)>1: main=main[main.grouping.str.contains(row['grouping'],case=False)]
    assert len(main)==1
    main=main.iloc[0]
    for key in ['ci_lower','ci_upper','bootstrap_valid']:
        value=float(main[key]);own=float(row[key])
        checks.append(dict(database=row['database'],grouping=row['grouping'],stratum=row['stratum'],metric=row['metric'],field=key,independent=own,main=value,absolute_difference=abs(own-value),match_1e12=bool(np.isclose(own,value,rtol=1e-12,atol=1e-12))))
pd.DataFrame(checks).to_csv(QA/'audit_bootstrap_comparison.csv',index=False)
print(json.dumps({'independent_intervals':len(outputs),'comparison_cells':len(checks),'nonidentical_cells':sum(not c['match_1e12'] for c in checks),'largest_interval_difference':max(c['absolute_difference'] for c in checks if c['field']!='bootstrap_valid')}))
