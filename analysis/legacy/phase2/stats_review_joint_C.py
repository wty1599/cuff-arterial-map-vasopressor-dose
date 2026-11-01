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
"""Independently group patient/bin counts and replay the declared child streams for C review."""
from pathlib import Path
import sys,json
ROOT=Path('@CUFF_WORKSPACE@');
import pandas as pd,numpy as np
OUT=ROOT/'phase2';RES=OUT/'_restricted';P1=ROOT/'phase1';levels=['<=0.1','>0.1-0.25','>0.25-0.5','>0.5-1.0','>1.0'];bins=['60-64','55-59','<55']
children=np.random.SeedSequence(20261027).spawn(2);totals=[];tables={};raw={};checks=[]
for db,child in zip(['mimic','eicu'],children):
    d=pd.read_parquet(P1/'_restricted'/f'pairs_{db}.parquet');patients=pd.Index(sorted(d.patient.unique()));N=len(patients)
    low=d[d.art_map.lt(65)].copy();low['mapbin']=np.select([low.art_map.ge(60),low.art_map.ge(55)],['60-64','55-59'],default='<55');low['missed']=low.cuff_map.ge(65).astype(int)
    pop=low.groupby(['patient','mapbin'],observed=True).size().unstack(fill_value=0).reindex(index=patients,columns=bins,fill_value=0).to_numpy()
    matrices=[pop];keys=[]
    scenarios=[('primary_5min',low,levels)]
    if db=='eicu':
        rec=pd.read_csv(P1/'_restricted/eicu_recovered_pairs.csv',usecols=['cid','recovered_nee']);z=low.merge(rec,on='cid',validate='one_to_one')
        ok=z.nee_raw.isna()&z.recovered_nee.notna()&~z.dose_weight_flag&z.recovered_nee.between(0,5)
        z.loc[ok,'nee']=z.loc[ok,'recovered_nee'];z['dose_stratum']=pd.cut(z.nee,[-np.inf,.1,.25,.5,1,np.inf],labels=levels).astype('string').fillna('unknown')
        scenarios.append(('recovered_dose',z,['all','>1.0']))
    offset=3
    for scenario,z,labs in scenarios:
        for lab in labs:
            u=z if lab=='all' else z[z.dose_stratum.eq(lab)]
            counts=u.groupby(['patient','mapbin'],observed=True).agg(num=('missed','sum'),den=('missed','size'))
            mat=counts.unstack(fill_value=0).reindex(index=patients,columns=pd.MultiIndex.from_product([['num','den'],bins]),fill_value=0).to_numpy()
            matrices.append(mat);keys.append((scenario,lab,offset));offset+=6
    matrix=np.concatenate(matrices,axis=1);rng=np.random.default_rng(child);bb=[]
    for _ in range(40):bb.append(rng.multinomial(N,np.repeat(1/N,N),size=50)@matrix)
    out=np.concatenate(bb);totals.append(out[:,:3]);tables[db]=(out,keys);raw[db]=matrix.sum(axis=0)
weight=sum(totals).astype(float);weight/=weight.sum(axis=1)[:,None]
expected=pd.read_csv(OUT/'phase2_C_standardized.csv')
saved=np.load(RES/'C_bootstrap_standard_weights.npz')['weights']
checks.append(dict(database='both',scenario='joint',stratum='standard',field='all_2000_standard_weights',difference=float(np.max(np.abs(weight-saved))),passed=bool(np.allclose(weight,saved,atol=1e-14))))
for db,(b,keys) in tables.items():
    for scenario,lab,off in keys:
        rate=np.divide(b[:,off:off+3],b[:,off+3:off+6],out=np.full((2000,3),np.nan),where=b[:,off+3:off+6]>0)
        ans=(weight*rate).sum(axis=1);valid=np.isfinite(ans);ci=np.quantile(ans[valid],[.025,.975]);row=expected[expected.database.eq(db)&expected.scenario.eq(scenario)&expected.stratum.eq(lab)].iloc[0]
        for field,value in [('ci_lower',ci[0]),('ci_upper',ci[1]),('bootstrap_valid',valid.sum())]:
            diff=abs(row[field]-value);checks.append(dict(database=db,scenario=scenario,stratum=lab,field=field,difference=float(diff),passed=bool(diff<1e-12)))
sm=pd.read_csv(OUT/'phase2_sensitivity_same_minute.csv')
for db in ['mimic','eicu']:
    d=pd.read_parquet(P1/'_restricted'/f'pairs_{db}.parquet');s=pd.read_parquet(P1/'_restricted'/f'same_minute_{db}.parquet');d=d.drop(columns=['art_map','art_values','nearest_diff']).merge(s,on='cid',validate='one_to_one')
    for lab in ['all','>1.0']:
        z=d if lab=='all' else d[d.dose_stratum.eq(lab)];rows=sm[sm.database.eq(db)&sm.domain.eq('B')&sm.stratum.eq(lab)]
        for row in rows.itertuples():
            for field,value in [('layer_low_arterial_pairs',z.art_map.lt(65).sum()),('layer_pairs',len(z)),('layer_low_arterial_prevalence',z.art_map.lt(65).mean())]:
                diff=abs(getattr(row,field)-value);checks.append(dict(database=db,scenario='same_minute',stratum=lab,field=field,difference=float(diff),passed=bool(diff<1e-12)))
out=pd.DataFrame(checks);out.to_csv(OUT/'phase2_statistical_review_repaired_paths.csv',index=False,encoding='utf-8-sig')
print(json.dumps({'checks':len(out),'failed':int((~out.passed).sum()),'maximum_difference':float(out.difference.max())}))
