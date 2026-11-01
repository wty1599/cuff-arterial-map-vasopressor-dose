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
import sys,json,os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
ROOT=Path('@CUFF_WORKSPACE@')
P1=ROOT/'phase1';P0=ROOT/'phase0';SRC=P0/'_restricted';INP=P1/'_restricted'
OUT=ROOT/'phase2';RES=OUT/'_restricted';LOG=OUT/'_logs'
for p in (OUT,RES,LOG,OUT/'figures'):p.mkdir(parents=True,exist_ok=True)


import duckdb,numpy as np,pandas as pd
P=json.loads(Path(__file__).with_name('phase2_parameters.json').read_text(encoding='utf-8'))
DBS=['mimic','eicu','amsterdam','sicdb'];MAIN=['mimic','eicu']
LABELS=['<=0.1','>0.1-0.25','>0.25-0.5','>0.5-1.0','>1.0']
ALLSTRATA=['all','known','unknown',*LABELS]
NAMES={'mimic':'MIMIC-IV','eicu':'eICU','amsterdam':'AmsterdamUMCdb','sicdb':'SICdb'}
def save(d,name):d.to_csv(OUT/name,index=False,encoding='utf-8-sig')
def con(db):
 c=duckdb.connect(str(SRC/f'{db}.duckdb'),read_only=True);c.execute("SET threads=1;SET memory_limit='700MB'");return c
def dose_strata(x):return pd.cut(x,[-np.inf,.1,.25,.5,1,np.inf],labels=LABELS).astype('string').fillna('unknown')
def read_pairs(db):return pd.read_parquet(INP/f'pairs_{db}.parquet')
def masks(d,groups=False):
 out={'all':np.ones(len(d),bool),'known':d.nee.notna().to_numpy(),'unknown':d.nee.isna().to_numpy()}
 out.update({s:d.dose_stratum.eq(s).to_numpy() for s in LABELS})
 if groups:out.update({'dosegroup:<=0.25':d.nee.le(.25).to_numpy(),'dosegroup:>0.25-1.0':(d.nee.gt(.25)&d.nee.le(1)).to_numpy(),'dosegroup:>1.0':d.nee.gt(1).to_numpy()})
 return out
