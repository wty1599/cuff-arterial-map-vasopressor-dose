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
import sys,json,csv,hashlib
R=Path('@CUFF_WORKSPACE@');A=R/'phase6'

import pandas as pd,numpy as np
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
checks=pd.read_csv(A/'phase6_independent_checks.csv');source=pd.read_csv(A/'phase6_source_audit_summary.csv')
rates=pd.read_csv(A/'phase6_window_rates.csv');flow=pd.read_csv(A/'phase6_flow.csv');profile=pd.read_csv(A/'phase6_profile.csv')
notes=[]
for db in ['mimic','eicu']:
    d=pd.read_parquet(A/'_restricted'/f'aligned_{db}.parquet');cm=d[d.common_anchor]
    notes.append(dict(database=db,common_pairs=len(cm),same_value_2min_same_minute=int((cm.art_2min==cm.art_same_minute).sum()),changed_value_5min_same_minute=int((cm.art_5min!=cm.art_same_minute).sum()),common_dose_known=int(cm.nee.notna().sum()),native_pairs=len(d),common_retention=len(cm)/len(d)))
pd.DataFrame(notes).to_csv(A/'phase6_independent_timing_checks.csv',index=False,encoding='utf-8-sig')
