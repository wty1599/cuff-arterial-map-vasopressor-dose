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
"""Independent ICM rev4 point estimates and fresh high-dose source sample.

No production analysis code is imported. Patient-level evidence stays restricted.
"""
from pathlib import Path
import sys, json, hashlib, datetime

import pandas as pd
import numpy as np
import duckdb
ROOT=Path('@CUFF_WORKSPACE@')
BASE=ROOT/'phase3/icm_revision4'
QA=BASE/'qa'; R=BASE/'analysis/_restricted'; SQL=ROOT/'sql/icm_rev4'
for p in [QA,R,SQL]: p.mkdir(parents=True,exist_ok=True)
STRATA=[('five','<=0.1',0,.1),('five','>0.1-0.25',.1,.25),
        ('five','>0.25-0.5',.25,.5),('five','>0.5-1.0',.5,1),
        ('five','>1.0',1,np.inf),('three','<=0.25',0,.25),
        ('three','>0.25-1.0',.25,1),('three','>1.0',1,np.inf)]
allrows=[]; qc=[]; samples=[]
for db in ['mimic','eicu']:
    p=ROOT/f'phase1/_restricted/pairs_{db}.parquet'
    d=pd.read_parquet(p)
    assert d.win.eq(5).all(), 'Expected primary 5-minute pairs only'
    assert d[['sid','t']].duplicated().sum()==0
    known=d.nee.notna()
    q={'database':db,'pairs':len(d),'patients':d.pid.nunique(),
       'known_pairs':int(known.sum()),'unknown_pairs':int((~known).sum()),
       'known_patients':d.loc[known,'pid'].nunique(),
       'art_lt65_pairs':int(d.art_map.lt(65).sum()),
       'art_lt65_rate':float(d.art_map.lt(65).mean()),
       'cuff_or_art_missing':int(d[['cuff_map','art_map']].isna().any(axis=1).sum()),
       'duplicates_sid_t':int(d[['sid','t']].duplicated().sum())}
    q['known_art_lt65_rate']=float(d.loc[known,'art_map'].lt(65).mean())
    q['known_unknown_conserved']=q['known_pairs']+q['unknown_pairs']==q['pairs']
    five_n=0
    for grouping,stratum,lo,hi in STRATA:
        s=d[known & d.nee.gt(lo) & d.nee.le(hi)]
        pos=s.art_map.lt(65); cuffpos=s.cuff_map.lt(65)
        TP=int((pos&cuffpos).sum()); FN=int((pos&~cuffpos).sum())
        FP=int((~pos&cuffpos).sum()); TN=int((~pos&~cuffpos).sum())
        a1den=int(s.cuff_map.ge(70).sum()); a1num=int((pos&s.cuff_map.ge(70)).sum())
        den65=int(s.cuff_map.ge(65).sum()); denbin=int((s.cuff_map.ge(65)&s.cuff_map.lt(70)).sum())
        lrminus=(FN/(TP+FN))/(TN/(TN+FP))
        lrplus=(TP/(TP+FN))/(FP/(TN+FP))
        allrows.append(dict(database=db,grouping=grouping,stratum=stratum,
            pairs=len(s),patients=s.pid.nunique(),TP=TP,FN=FN,FP=FP,TN=TN,
            numerator=a1num,denominator=a1den,estimate=a1num/a1den,
            LR_negative=lrminus,LR_positive=lrplus,
            denominator_ge65=den65,denominator_65to69=denbin,
            denominator_partition_pass=a1den==den65-denbin,
            confusion_conservation_pass=TP+FN+FP+TN==len(s)))
        if grouping=='five': five_n+=len(s)
    q['five_strata_pairs']=five_n
    q['five_strata_conserved']=five_n==int(known.sum())
    qc.append(q)
    high=d[d.nee.gt(1)&d.nee.le(5)].copy()
    # Exclude prior independent source samples; selecting a new deterministic 50.
    prior=set()
    for old in [ROOT/f'phase1/_restricted/audit_phase1_{db}.csv',
                ROOT/f'phase2/_restricted/audit_phase2_{db}.csv',
                ROOT/f'phase6/_restricted/audit_sample_{db}.csv']:
        if old.exists():
            x=pd.read_csv(old,dtype={'sid':str})
            if {'sid','t'}<=set(x): prior.update(zip(x.sid.astype(str),x.t.astype(float)))
    high=high[[ (str(s),float(t)) not in prior for s,t in zip(high.sid,high.t) ]].copy()
    high['audit_hash']=[hashlib.sha256(f'ICM-rev4-new-source-50-20261031|{db}|{s}|{float(t):.6f}'.encode()).hexdigest() for s,t in zip(high.sid,high.t)]
    sample=high.sort_values('audit_hash').head(50).copy()
    assert len(sample)==50
    con=duckdb.connect(str(ROOT/f'phase0/_restricted/{db}.duckdb'),read_only=True)
    ep=con.execute('SELECT cid,episode_start,episode_stop FROM cuff60').df();con.close()
    sample=sample.merge(ep,on='cid',how='left',validate='one_to_one')
    assert sample.episode_start.notna().all()
    sample['audit_id']=[f'{db}_{i+1:03d}' for i in range(len(sample))]
    sample['expected_phase0_weight']=sample.weight_raw
    sample['expected_phase0_height']=sample.height_raw
    sample.to_csv(R/f'audit_sample_{db}.csv',index=False)
    samples.append({'database':db,'sample_pairs':len(sample),'sample_patients':sample.pid.nunique(),
      'previous_sample_overlap':0,'sample_sha256':hashlib.sha256((R/f'audit_sample_{db}.csv').read_bytes()).hexdigest()})
pd.DataFrame(allrows).to_csv(QA/'audit_point_estimates.csv',index=False)
pd.DataFrame(qc).to_csv(QA/'audit_cohort_qc.csv',index=False)
pd.DataFrame(samples).to_csv(QA/'audit_sample_receipt.csv',index=False)
# Reuse only the independently authored source auditor, redirect every output.
# The production extraction/analysis modules remain unimported.
src=(ROOT/'sql/phase1/audit_pg_source.py').read_text(encoding='utf-8')
src=src.replace("SQLDIR = ROOT / 'sql' / 'phase1'", "SQLDIR = ROOT / 'sql' / 'icm_rev4'")
src=src.replace("OUT = ROOT / 'phase1'", "OUT = ROOT / 'phase3' / 'icm_revision4' / 'analysis'")
src=src.replace('from audit_quality import P as QUALITY, anthro, mimic_anthro',
    "import sys\nsys.path.insert(0,str(ROOT/'sql/phase1'))\nfrom audit_quality import P as QUALITY, anthro, mimic_anthro")
src=src.replace('audit_phase1_', 'audit_sample_')
src=src.replace("SQLDIR/'phase1_parameters.json'", "ROOT/'phase3/icm_revision4/analysis/icm_rev4_params.yaml'")
src=src.replace("phase1/_restricted/audit_failure.log", "phase3/icm_revision4/analysis/_restricted/audit_failure.log")
src=src.replace('阶段0结果之后；新50对样本固定后执行。阶段1主分层结果相对时序见 audit_execution_timeline.csv',
                'ICM rev4主分析结果已知之后；新高剂量50对样本写定后；独立重查真实源数据库。')
(SQL/'audit_source_local.py').write_text(src,encoding='utf-8')
# Only now compare the independently constructed point table with main outputs.
a=pd.DataFrame(allrows)
maina=pd.read_csv(BASE/'analysis/A1_cuff_ge70.csv')
mainb=pd.read_csv(BASE/'analysis/A2_likelihood_ratios.csv')
checks=[]
for _,row in a.iterrows():
    # Main grouping labels are selected from their actual strata membership.
    ma=maina[(maina.database==row.database)&(maina.stratum==row.stratum)]
    mb=mainb[(mainb.database==row.database)&(mainb.stratum==row.stratum)]
    if len(ma)>1:
        ma=ma[ma.grouping.str.contains('five' if row.grouping=='five' else 'three',case=False)]
        mb=mb[mb.grouping.str.contains('five' if row.grouping=='five' else 'three',case=False)]
    assert len(ma)==1, 'Main A1 comparator key unresolved'
    ma=ma.iloc[0]
    for key in ['numerator','denominator','estimate']:
        checks.append(dict(database=row.database,grouping=row.grouping,stratum=row.stratum,metric='A1_'+key,
            independent=float(row[key]),main=float(ma[key]),match=bool(np.isclose(row[key],ma[key],rtol=1e-12,atol=1e-12))))
    for metric,own in [('LR_negative',row.LR_negative),('LR_positive',row.LR_positive)]:
        mm=mb[mb.metric.eq(metric)]
        if mm.empty: mm=mb[mb.metric.str.contains('negative' if metric=='LR_negative' else 'positive',case=False)]
        assert len(mm)==1, 'Main A2 comparator key unresolved'
        value=float(mm.iloc[0].estimate)
        checks.append(dict(database=row.database,grouping=row.grouping,stratum=row.stratum,metric=metric,
            independent=float(own),main=value,match=bool(np.isclose(own,value,rtol=1e-12,atol=1e-12))))
pd.DataFrame(checks).to_csv(QA/'audit_point_comparison.csv',index=False)
receipt={'created_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'point_checks':len(checks),'point_mismatches':sum(not x['match'] for x in checks),
    'source_samples':samples,'source_queries_executed':False,
    'independent_algorithm':'Direct four-cell counts and threshold denominators; no production imports.',
    'timing':'Parent results had already been reported; independent counts computed before reading comparator CSV.',
    'parameter_sha256':hashlib.sha256((BASE/'analysis/icm_rev4_params.yaml').read_bytes()).hexdigest()}
(QA/'audit_prepare_receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(receipt,ensure_ascii=False))
