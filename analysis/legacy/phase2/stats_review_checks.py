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
"""Bounded independent aggregate review; does not call production functions or raw source audits."""
from pathlib import Path
import sys,json
ROOT=Path('@CUFF_WORKSPACE@')

import pandas as pd,numpy as np
OUT=ROOT/'phase2';RES=OUT/'_restricted';P1=ROOT/'phase1'
checks=[]
def check(name,actual,expected,db='',key='',tol=1e-12):
    passed=(pd.isna(actual) and pd.isna(expected)) or (not pd.isna(actual) and not pd.isna(expected) and abs(float(actual)-float(expected))<=tol)
    checks.append(dict(check=name,database=db,key=key,actual=actual,expected=expected,passed=bool(passed)))
def select(d,lab):
    if lab=='all':return np.ones(len(d),bool)
    if lab=='known':return d.nee.notna().to_numpy()
    if lab=='unknown':return d.nee.isna().to_numpy()
    if lab=='dosegroup:<=0.25':return d.nee.le(.25).to_numpy()
    if lab=='dosegroup:>0.25-1.0':return (d.nee.gt(.25)&d.nee.le(1)).to_numpy()
    if lab=='dosegroup:>1.0':return d.nee.gt(1).to_numpy()
    return d.dose_stratum.eq(lab).to_numpy()
dist=pd.read_csv(OUT/'phase2_A_miss_distribution.csv');A=pd.read_csv(OUT/'phase2_A_deep_miss.csv');B=pd.read_csv(OUT/'phase2_B_cuff_perspective.csv');C=pd.read_csv(OUT/'phase2_C_cross_table.csv')
standard=pd.read_csv(OUT/'phase2_C_standard_weights.csv');std=pd.read_csv(OUT/'phase2_C_standardized.csv');interp=pd.read_csv(OUT/'phase2_C_interpretation.csv')
data={};total_std=np.zeros(3);bins=['60-64','55-59','<55']
for db in ['mimic','eicu','amsterdam','sicdb']:
    d=pd.read_parquet(P1/'_restricted'/f'pairs_{db}.parquet');data[db]=d
    for row in A[A.database.eq(db)].itertuples():
        mask=select(d,row.stratum);threshold=55 if row.metric.endswith('55') else 60;den=mask&d.art_map.lt(threshold).to_numpy();ev=den&d.cuff_map.ge(65).to_numpy()
        key=f'{row.metric};{row.stratum}';check('A_denominator',row.denominator,int(den.sum()),db,key);check('A_numerator',row.numerator,int(ev.sum()),db,key)
    for lab,z in dist[dist.database.eq(db)].groupby('stratum'):
        check('A_distribution_numerator_conservation',z.numerator.sum(),z.denominator.iloc[0],db,lab)
        check('A_distribution_denominator_constant',z.denominator.nunique(),1,db,lab)
    if db not in ['mimic','eicu']:continue
    a=d.art_map.to_numpy();cf=d.cuff_map.to_numpy();low=a<65
    abins={'60-64':(a>=60)&(a<65),'55-59':(a>=55)&(a<60),'<55':a<55}
    total_std+=np.array([abins[x].sum() for x in bins])
    for row in B[B.database.eq(db)].itertuples():
        mask=select(d,row.stratum);den=mask.copy();outcome=low
        if row.metric=='art_lt65_given_cuff_ge65':den&=cf>=65
        elif row.metric=='art_lt65_given_cuff_bin':
            ca={'65-69':(cf>=65)&(cf<70),'70-74':(cf>=70)&(cf<75),'75-79':(cf>=75)&(cf<80),'>=80':cf>=80}[row.category];den&=ca
        elif row.metric=='art_lt60_given_cuff65_69':den&=(cf>=65)&(cf<70);outcome=a<60
        elif row.metric=='art_lt65_given_cuff_lt65':den&=cf<65
        key=f'{row.metric};{row.stratum};{row.category}'
        check('B_denominator',row.denominator,int(den.sum()),db,key);check('B_numerator',row.numerator,int((den&outcome).sum()),db,key)
        check('B_layer_low',row.layer_low_arterial_pairs,int((mask&low).sum()),db,key);check('B_layer_total',row.layer_pairs,int(mask.sum()),db,key)
        check('B_condition_prevalence',row.condition_low_arterial_prevalence,(den&low).sum()/den.sum() if den.sum() else np.nan,db,key)
    for lab,z in B[B.database.eq(db)].groupby('stratum'):
        n=z[z.metric.eq('art_lt65_given_cuff_ge65')].iloc[0];q=z[z.metric.eq('art_lt65_given_cuff_bin')]
        check('B_cuff_bin_denominator_conservation',q.denominator.sum(),n.denominator,db,lab);check('B_cuff_bin_numerator_conservation',q.numerator.sum(),n.numerator,db,lab)
    for row in C[C.database.eq(db)].itertuples():
        mask=select(d,row.stratum);den=mask&abins[row.category];key=f'{row.stratum};{row.category}'
        check('C_cross_denominator',row.denominator,int(den.sum()),db,key);check('C_cross_numerator',row.numerator,int((den&(cf>=65)).sum()),db,key)
    for lab,z in C[C.database.eq(db)].groupby('stratum'):
        check('C_bin_denominator_conservation',z.denominator.sum(),int((select(d,lab)&low).sum()),db,lab)
weights=total_std/total_std.sum()
for i,b in enumerate(bins):
    r=standard[standard.arterial_bin.eq(b)].iloc[0];check('C_standard_population_count',r.n_pairs,total_std[i],key=b);check('C_standard_population_weight',r.weight,weights[i],key=b)
for row in std[std.scenario.eq('primary_5min')].itertuples():
    q=C[C.database.eq(row.database)&C.stratum.eq(row.stratum)].set_index('category').loc[bins]
    p=np.sum(weights*q.numerator.to_numpy()/q.denominator.to_numpy());check('C_standardized_point',row.estimate,p,row.database,row.stratum)
for row in interp.itertuples():
    z=std[std.database.eq(row.database)&std.scenario.eq('primary_5min')]
    span=z.estimate.max()-z.estimate.min();raw=z.raw_missed_rate.max()-z.raw_missed_rate.min()
    check('C_standardized_span',row.standardized_span,span,row.database);check('C_5pp_rule',int(row.rule_satisfied),int(span<=.05 and span<raw),row.database)
# Verify percentile intervals directly from saved aggregate replicates, not by rerunning model/data extraction.
for db in data:
    specs=pd.read_csv(RES/f'ratio_specs_{db}_primary_5min.csv');boot=np.load(RES/f'bootstrap_{db}_primary_5min.npz')
    ne=boot['numerator'];de=boot['denominator'];br=np.divide(ne,de,out=np.full_like(ne,np.nan),where=de>0)
    for row in specs.itertuples():
        i=int(row.bootstrap_column);good=np.isfinite(br[:,i]);ll,uu=np.quantile(br[good,i],[.025,.975]) if good.any() else (np.nan,np.nan)
        check('percentile_CI_lower',row.ci_lower,ll,db,str(i));check('percentile_CI_upper',row.ci_upper,uu,db,str(i));check('bootstrap_replicates',ne.shape[0],2000,db)
result=pd.DataFrame(checks);result.to_csv(OUT/'phase2_statistical_review_checks.csv',index=False,encoding='utf-8-sig')
print(json.dumps({'checks':len(result),'failed':int((~result.passed).sum()),'standard_total':int(total_std.sum())}))
# This descriptive correlation does not prove independent RNG streams.
bs=[]
for db in ['mimic','eicu']:
    s=pd.read_csv(RES/f'ratio_specs_{db}_primary_5min.csv');z=s[s.domain.eq('C_cross')&s.stratum.eq('all')].set_index('category').loc[bins]
    b=np.load(RES/f'bootstrap_{db}_primary_5min.npz');bs.append(b['denominator'][:,z.bootstrap_column.to_numpy(int)])
corr=np.corrcoef(np.concatenate(bs,axis=1).T)[:3,3:]
pd.DataFrame(corr,index=['mimic_'+b for b in bins],columns=['eicu_'+b for b in bins]).to_csv(OUT/'phase2_statistical_review_cross_stream_correlation.csv',encoding='utf-8-sig')
print(json.dumps({'cross_bin_replicate_correlations_max_abs':float(np.abs(corr).max())}))
