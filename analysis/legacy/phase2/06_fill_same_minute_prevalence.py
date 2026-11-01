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
q=pd.read_csv(OUT/'phase2_sensitivity_same_minute.csv')
for db in MAIN:
    d=read_pairs(db);sm=pd.read_parquet(INP/f'same_minute_{db}.parquet');d=d.drop(columns=['art_map','art_values','nearest_diff']).merge(sm,on='cid',validate='one_to_one')
    for lab in ['all','>1.0']:
        mask=np.ones(len(d),bool) if lab=='all' else d.dose_stratum.eq(lab);den=int(mask.sum());num=int((mask&d.art_map.lt(65)).sum());ix=q.database.eq(db)&q.stratum.eq(lab)
        q.loc[ix,'layer_low_arterial_pairs']=num;q.loc[ix,'layer_pairs']=den;q.loc[ix,'layer_low_arterial_prevalence']=num/den if den else np.nan
        ix=ix&q.domain.eq('B');q.loc[ix,'condition_low_arterial_pairs']=q.loc[ix,'numerator'];q.loc[ix,'condition_low_arterial_prevalence']=q.loc[ix,'estimate']
save(q,'phase2_sensitivity_same_minute.csv');print('SAME_MINUTE_PREVALENCE_COMPLETE_NO_RATE_CHANGES')
