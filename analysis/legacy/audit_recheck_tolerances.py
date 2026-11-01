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
from pathlib import Path
import pandas as pd
O=Path('@CUFF_WORKSPACE@/phase0');R=O/'_restricted'
db='amsterdam';p=R/f'audit_file_comparison_{db}.csv';x=pd.read_csv(p)
checks={'cuff_map':'source_cuff','art_map':'source_art','nee':'source_nee','nearest_diff':'source_diff','episode_start':'source_start','episode_stop':'source_stop'}
for old,new in checks.items():
    a=pd.to_numeric(x[old],errors='coerce');b=pd.to_numeric(x[new],errors='coerce');tol=1e-4 if old in ['nearest_diff','episode_start','episode_stop'] else 1e-6*(1+a.abs());x[old+'_match']=(a.isna()&b.isna())|((a-b).abs()<=tol)
x['all_match']=x[[k+'_match' for k in checks]].all(axis=1);x.to_csv(p,index=False)
rows=[{'metric':'sampled_pairs','value':len(x)},{'metric':'any_discrepancy_pairs','value':int((~x.all_match).sum())}]+[{'metric':k+'_discrepancies','value':int((~x[k+'_match']).sum())} for k in checks]
pd.DataFrame(rows).to_csv(O/f'phase0_source_audit_{db}.csv',index=False);print('Absolute-time tolerance audit:',int((~x.all_match).sum()),'discrepancies')
