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

from common import *
checks=[];artifacts=[]
def ck(db,name,passed,detail=''):
 checks.append(dict(database=db,check=name,passed=bool(passed),details=detail))
A=pd.read_csv(OUT/'phase2_A_deep_miss.csv');AD=pd.read_csv(OUT/'phase2_A_miss_distribution.csv');B=pd.read_csv(OUT/'phase2_B_cuff_perspective.csv');C=pd.read_csv(OUT/'phase2_C_cross_table.csv');S=pd.read_csv(OUT/'phase2_C_standardized.csv');W=pd.read_csv(OUT/'phase2_C_standard_weights.csv').set_index('arterial_bin').weight
for db in DBS:
 d=read_pairs(db);p1=pd.read_csv(P1/f'phase1_primary_{db}.csv').set_index('stratum')
 ck(db,'pair_count_unchanged',len(d)==int(p1.loc['all','all_pairs']))
 for lab,mask in masks(d).items():
  z=d[mask];miss=z.art_map.lt(65)&z.cuff_map.ge(65)
  q=AD[AD.database.eq(db)&AD.stratum.eq(lab)]
  ck(db,'A_missed_bins_sum_'+lab,int(q.numerator.sum())==int(miss.sum()) and q.denominator.nunique()==1 and int(q.denominator.iloc[0])==int(miss.sum()))
  for th in [55,60]:
   r=A[A.database.eq(db)&A.stratum.eq(lab)&A.metric.eq(f'deep_miss_art_lt{th}')].iloc[0]
   ck(db,f'A_deep_counts_{lab}_{th}',r.denominator==z.art_map.lt(th).sum() and r.numerator==(z.art_map.lt(th)&z.cuff_map.ge(65)).sum())
  if db in MAIN:
   q=B[B.database.eq(db)&B.stratum.eq(lab)&B.metric.eq('art_lt65_given_cuff_bin')]
   r=B[B.database.eq(db)&B.stratum.eq(lab)&B.metric.eq('art_lt65_given_cuff_ge65')].iloc[0]
   ck(db,'B_cuff_bins_sum_'+lab,q.denominator.sum()==r.denominator and q.numerator.sum()==r.numerator)
   pr=B[B.database.eq(db)&B.stratum.eq(lab)&B.metric.eq('arterial_lt65_prevalence')].iloc[0]
   ck(db,'B_prevalence_'+lab,pr.denominator==len(z) and pr.numerator==z.art_map.lt(65).sum())
   q=C[C.database.eq(db)&C.stratum.eq(lab)]
   ck(db,'C_arterial_bins_sum_'+lab,q.denominator.sum()==z.art_map.lt(65).sum() and q.numerator.sum()==miss.sum())
   if lab in LABELS:
    calc=sum(W.loc[r.category]*r.numerator/r.denominator for r in q.itertuples()) if q.denominator.gt(0).all() else np.nan
    rr=S[S.database.eq(db)&S.scenario.eq('primary_5min')&S.stratum.eq(lab)].iloc[0]
    ck(db,'C_direct_standardization_'+lab,np.isclose(calc,rr.estimate,atol=1e-12,equal_nan=True))
  ck(db,'phase1_missed_count_'+lab,p1.loc[lab,'numerator']==miss.sum())
 for metric in A.metric.unique():
  q=A[A.database.eq(db)&A.metric.eq(metric)].set_index('stratum');ss=q.loc[[*LABELS,'unknown']]
  ck(db,'A_dose_conservation_'+metric,ss.numerator.sum()==q.loc['all','numerator'] and ss.denominator.sum()==q.loc['all','denominator'])
for path in OUT.glob('phase2_*.csv'):
 q=pd.read_csv(path)
 if {'numerator','denominator'}.issubset(q.columns):
  valid=q[['numerator','denominator']].notna().all(axis=1)
  ck(path.name,'numerator_le_denominator',(q.loc[valid,'numerator']<=q.loc[valid,'denominator']).all())
 if {'estimate','ci_lower','ci_upper'}.issubset(q.columns):
  good=q.estimate.notna();ck(path.name,'probabilities_0_to_1',q.loc[good,'estimate'].between(0,1).all())
  valid=q[['ci_lower','ci_upper']].notna().all(axis=1);ck(path.name,'CI_ordered',(q.loc[valid,'ci_lower']<=q.loc[valid,'ci_upper']).all())
 if {'boundary','precision'}.issubset(q.columns):
  z=q[q.boundary.eq(True)];ck(path.name,'boundary_label',z.precision.eq('boundary_not_informative').all())
save(pd.DataFrame(checks),'phase2_validation.csv')
if not all(x['passed'] for x in checks):raise RuntimeError('impossible_count_or_consistency_result')
print(json.dumps({'checks':len(checks),'all_passed':True}))
