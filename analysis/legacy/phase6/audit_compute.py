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
"""Independent DuckDB reconstruction of all rates; SQL patient counts for CI audit."""
from pathlib import Path
import sys,json,csv,hashlib,datetime
R=Path('@CUFF_WORKSPACE@');A=R/'phase6';I=R/'phase1/_restricted'

import duckdb,numpy as np,pandas as pd
P=json.loads((A/'phase6_params.yaml').read_text(encoding='utf-8'))
all_checks=[];all_points=[];cis=[]
def ck(db,name,ok,detail=''):
    all_checks.append(dict(database=db,check=name,passed=bool(ok),detail=str(detail)))
def esc(x):return "'"+str(x).replace("'","''")+"'"
for db in ['mimic','eicu']:
    c=duckdb.connect();c.execute("SET threads=2; SET memory_limit='1GB'")
    c.execute(f"CREATE TABLE x AS SELECT b.*,b.art_map a5,w.art_map a2,s.art_map a0 FROM read_parquet({esc((I/f'pairs_{db}.parquet').as_posix())}) b LEFT JOIN read_parquet({esc((I/f'window2_{db}.parquet').as_posix())}) w USING(cid) LEFT JOIN read_parquet({esc((I/f'same_minute_{db}.parquet').as_posix())}) s USING(cid)")
    ck(db,'unique_cuff_ids',c.execute('SELECT count(*)=count(DISTINCT cid) FROM x').fetchone()[0])
    ck(db,'same_minute_subset_2min',c.execute('SELECT count(*) FROM x WHERE a0 IS NOT NULL AND a2 IS NULL').fetchone()[0]==0)
    ck(db,'all_phase1_rows_have_5min_reference',c.execute('SELECT count(*) FROM x WHERE a5 IS NULL').fetchone()[0]==0)
    strata={'all':'true','known':"dose_stratum <> 'unknown'",'unknown':"dose_stratum = 'unknown'"}
    strata.update({v:f'dose_stratum={esc(v)}' for v in P['dose_strata']})
    points=[];ci_sql=[];ci_meta=[]
    for frame in ['native','common']:
      for win,a in [('5min','a5'),('2min','a2'),('same_minute','a0')]:
        elig=f'{a} IS NOT NULL' if frame=='native' else 'a5 IS NOT NULL AND a2 IS NOT NULL AND a0 IS NOT NULL'
        for st,cond in strata.items():
          selected=f'({elig}) AND ({cond})'
          defs={
            'missed':(f'{a}<65 AND cuff_map>=65',f'{a}<65'),
            'deep_missed':(f'{a}<55 AND cuff_map>=65',f'{a}<55'),
            'cuff_conditioned':(f'{a}<65 AND cuff_map>=65','cuff_map>=65'),
            'cuff_65_69':(f'{a}<65 AND cuff_map>=65 AND cuff_map<70','cuff_map>=65 AND cuff_map<70'),
            'cuff_65_69_deep':(f'{a}<60 AND cuff_map>=65 AND cuff_map<70','cuff_map>=65 AND cuff_map<70'),
            'hypotension_prevalence':(f'{a}<65','true')}
          expressions=['count(*) n','count(DISTINCT pid) p']
          for metric,(num,den) in defs.items():expressions.extend([f'count(*) FILTER(WHERE {num}) {metric}_num',f'count(*) FILTER(WHERE {den}) {metric}_den',f'count(DISTINCT pid) FILTER(WHERE {den}) {metric}_patients'])
          z=c.execute(f'SELECT '+','.join(expressions)+f' FROM x WHERE {selected}').df().iloc[0]
          for metric,(num,den) in defs.items():
            N,D=int(z[metric+'_num']),int(z[metric+'_den'])
            points.append(dict(database=db,frame=frame,window=win,stratum=st,metric=metric,numerator=N,denominator=D,denominator_patients=int(z[metric+'_patients']),all_pairs=int(z.n),all_patients=int(z.p),estimate=N/D if D else np.nan))
            # Every principal endpoint overall and highest dose, all three common windows.
            if frame=='common' and st in ['all','>1.0']:
              ix=len(ci_meta);ci_meta.append((frame,win,st,metric))
              ci_sql.extend([f'count(*) FILTER(WHERE ({selected}) AND ({num})) n{ix}',f'count(*) FILTER(WHERE ({selected}) AND ({den})) d{ix}'])
    actual=pd.read_csv(A/'phase6_window_rates.csv');actual=actual[actual.database.eq(db)]
    keys=['database','frame','window','stratum','metric']
    independent=pd.DataFrame(points);merge=independent.merge(actual,on=keys,suffixes=('_ind','_prod'),validate='one_to_one')
    ck(db,'all_expected_rate_cells',len(merge)==288 and len(actual)==288)
    for field in ['numerator','denominator','denominator_patients','all_pairs','all_patients','estimate']:
        x,y=merge[field+'_ind'],merge[field+'_prod'];delta=(x-y).abs();ok=(delta<1e-12)|(x.isna()&y.isna())
        ck(db,f'all_{field}_match',ok.all(),f'max_difference={delta.max()}, cells={len(merge)}')
    all_points.extend(points)
    # Independently assembled patient sufficient counts; reproduce all 2,000 multiplicities.
    q='SELECT pid,'+','.join(ci_sql)+' FROM x GROUP BY pid ORDER BY pid'
    counts=c.execute(q).df();n=len(counts);matrix=counts.iloc[:,1:].to_numpy(float)
    rng=np.random.default_rng(20261027);rep=[]
    for k in range(20):
        w=rng.multinomial(n,np.repeat(1/n,n),size=100);s=w@matrix
        with np.errstate(divide='ignore',invalid='ignore'):rep.append(s[:,::2]/s[:,1::2])
    rep=np.vstack(rep)
    prod_boot=np.load(A/'_restricted'/f'joint_bootstrap_{db}.npz')['estimates']
    for j,(frame,win,st,metric) in enumerate(ci_meta):
        z=actual[(actual.frame==frame)&(actual.window==win)&(actual.stratum==st)&(actual.metric==metric)].iloc[0]
        vals=rep[:,j];valid=np.isfinite(vals);l,u=np.quantile(vals[valid],[.025,.975]) if valid.any() else [np.nan,np.nan]
        v=prod_boot[:,int(z.bootstrap_column)];same=np.allclose(vals,v,atol=1e-12,rtol=0,equal_nan=True)
        ci_ok=np.allclose([l,u],[z.ci_lower,z.ci_upper],atol=1e-12,rtol=0,equal_nan=True)
        ck(db,f'bootstrap_{frame}_{win}_{st}_{metric}',same and ci_ok and valid.sum()==z.bootstrap_valid)
        cis.append(dict(database=db,frame=frame,window=win,stratum=st,metric=metric,ind_ci_lower=l,ind_ci_upper=u,bootstrap_valid=int(valid.sum()),all_2000_replicates_match=bool(same)))
    # Point algebra independently derives each three-part path from independently reconstructed rates.
    ix={(z['frame'],z['window'],z['stratum'],z['metric']):z['estimate'] for z in points}
    decom=pd.read_csv(A/'phase6_decomposition.csv');decom=decom[decom.database.eq(db)]
    for row in decom.itertuples():
        b=ix[('native','5min',row.stratum,row.metric)];s=ix[('common','5min',row.stratum,row.metric)];t=ix[('common',row.target_window,row.stratum,row.metric)]
        val={'total':t-b,'selection':s-b,'within_anchor_window':t-s}[row.component]
        ck(db,f'decomposition_{row.target_window}_{row.stratum}_{row.metric}_{row.component}',(np.isnan(val) and np.isnan(row.estimate)) or abs(val-row.estimate)<1e-12)
    c.close()
pd.DataFrame(all_points).to_csv(A/'phase6_independent_points.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(cis).to_csv(A/'phase6_independent_bootstrap.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(all_checks).to_csv(A/'phase6_independent_checks.csv',index=False,encoding='utf-8-sig')
(A/'phase6_independent_compute_receipt.json').write_text(json.dumps({'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':'original phase1 Parquet inputs; independent SQL expressions, no production module import','checks':len(all_checks),'passed':all(x['passed'] for x in all_checks),'points':len(all_points),'bootstrap_cells':len(cis),'all_2000_replicates_per_cell':True},indent=2),encoding='utf-8')
print(json.dumps({'checks':len(all_checks),'passed':all(x['passed'] for x in all_checks),'points':len(all_points),'bootstrap_cells':len(cis)}))
