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
"""Aggregate-only independent receipts for ICM rev4, preserving previous work."""
from pathlib import Path
import sys,json,hashlib,datetime,re

import pandas as pd
import numpy as np
ROOT=Path('@CUFF_WORKSPACE@'); BASE=ROOT/'phase3/icm_revision4'; QA=BASE/'qa'; R=BASE/'analysis/_restricted'
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
params=json.loads((BASE/'analysis/_restricted/validation_expected_checks.json').read_text(encoding='utf-8'))
points=pd.read_csv(QA/'audit_point_estimates.csv'); comps=pd.read_csv(QA/'audit_point_comparison.csv')
boots=pd.read_csv(QA/'audit_bootstrap_comparison.csv'); cohort=pd.read_csv(QA/'audit_cohort_qc.csv')
expect=[]
for db in ['mimic','eicu']:
    for gr,st,key in [('three','<=0.25',db+'_low'),('five','>1.0',db+'_high')]:
        row=points[(points.database==db)&(points.grouping==gr)&(points.stratum==st)].iloc[0]
        num,den=params['quality']['A1_expected_counts'][key]
        expect.append(dict(database=db,check='user_A1_'+st,pass_check=int(row.numerator)==num and int(row.denominator)==den,independent=f'{int(row.numerator)}/{int(row.denominator)}',expected=f'{num}/{den}'))
    five=points[(points.database==db)&(points.grouping=='five')]
    for (_,row),ex in zip(five.iterrows(),params['quality']['A2_expected_rounded'][db]):
        expect.append(dict(database=db,check='user_LR_negative_'+row.stratum,pass_check=round(row.LR_negative,2)==ex,independent=round(row.LR_negative,2),expected=ex))
pd.DataFrame(expect).to_csv(QA/'audit_user_expected_checks.csv',index=False)
# Reused baseline QC is explicit, not presented as new extraction.
t1path=ROOT/'phase2/phase2_table1.csv'; t1=pd.read_csv(t1path)
t1=t1[t1.database.isin(['mimic','eicu'])].drop_duplicates(['database','variable'])
t1['qc_status']='复用阶段2原主队列Table1；本轮未重新提取基线'
t1['source_sha256']=sha(t1path)
t1.to_csv(QA/'audit_table1_missingness_reused.csv',index=False)
source=[]
for db in ['mimic','eicu']:
    path=R/f'audit_summary_{db}.json'
    if not path.exists(): raise RuntimeError(f'{db} actual source query has not completed')
    s=json.loads(path.read_text(encoding='utf-8'))
    detail=pd.read_csv(R/f'audit_comparison_{db}.csv')
    clean_match=np.isclose(detail.source_clean_nee,detail.nee,rtol=1e-6,atol=1e-6,equal_nan=True)
    count=lambda name: len(pd.read_csv(R/f'audit_source_{db}_{name}.csv'))
    source.append(dict(database=db,sample_pairs=s['sample_pairs'],sample_patients=detail.sid.nunique(),
      source_bp_records=count('bp'),source_drug_records=count('drugs'),
      cuff_art_time_episode_source_discrepancy_pairs=s['discrepancy_pairs'],
      quality_clean_nee_discrepancy_pairs=int((~clean_match).sum()),
      known_reconstructed_dose_pairs=int(detail.source_clean_nee.notna().sum()),
      actual_postgres_query_executed=s['source_queries_executed'],
      source_query_started_at=s['source_query_started_at'],source_query_completed_at=s['source_query_completed_at'],
      sample_sha256=s['sample_sha256'],source_query_sql_sha256=sha(ROOT/f'sql/icm_rev4/audit_source_{db}.sql'),
      pass_check=bool(s['sample_pairs']==50 and s['discrepancy_pairs']==0 and clean_match.all())))
pd.DataFrame(source).to_csv(QA/'audit_source_summary.csv',index=False)
