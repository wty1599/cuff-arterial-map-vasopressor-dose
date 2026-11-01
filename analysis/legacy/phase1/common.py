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

from pathlib import Path
import sys, json
ROOT=Path('@CUFF_WORKSPACE@')
SRC=ROOT/'phase0'/'_restricted'
OUT=ROOT/'phase1'
RES=OUT/'_restricted'
LOG=OUT/'_logs'
for p in (OUT,RES,LOG,OUT/'figures'): p.mkdir(parents=True,exist_ok=True)


import duckdb, numpy as np, pandas as pd
P=json.loads((Path(__file__).parent/'phase1_parameters.json').read_text(encoding='utf-8'))
DBS=['mimic','eicu','amsterdam','sicdb']
LABELS=['<=0.1','>0.1-0.25','>0.25-0.5','>0.5-1.0','>1.0']
def con(db):
    c=duckdb.connect(str(SRC/f'{db}.duckdb'),read_only=True)
    c.execute("SET memory_limit='700MB'; SET threads=1")
    return c
def save(d,name):
    d.to_csv(OUT/name,index=False,encoding='utf-8-sig')
def strata(x):
    return pd.cut(x,[-np.inf,*P['dose_edges'],np.inf],labels=LABELS).astype('string').fillna('unknown')
