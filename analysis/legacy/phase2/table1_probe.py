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
import sys,json,gzip
ROOT=Path('@CUFF_WORKSPACE@');O=ROOT/'phase2';R=O/'_restricted';R.mkdir(exist_ok=True,parents=True)

import duckdb,pandas as pd
for db in ['mimic','eicu','amsterdam','sicdb']:
 c=duckdb.connect(str(ROOT/f'phase0/_restricted/{db}.duckdb'),read_only=True);c.execute("SET memory_limit='700MB';SET threads=1")
 p=c.execute(f"SELECT * FROM read_parquet('{(ROOT/f'phase1/_restricted/pairs_{db}.parquet').as_posix()}')").df();pt=p.drop_duplicates('patient');pt[['sid']].to_csv(R/f'table1_ids_{db}.csv',index=False)
 print(json.dumps({'db':db,'pair_columns':list(p.columns),'stays_columns':c.execute('DESCRIBE stays').df().column_name.tolist(),'patients':len(pt),'sids':pt.sid.nunique(),'pairs':len(p),'duplicate_patient_multiple_sid':int(p.groupby('patient').sid.nunique().gt(1).sum())}))
 print(c.execute('SELECT dose_reason,count(*) n FROM drugs GROUP BY 1').df().to_json(orient='records'))
 c.close()
src=Path('@SICDB_ROOT@')
s=pd.read_csv(src/'cases.csv.gz',dtype=str,keep_default_na=False);ref=pd.read_csv(src/'d_references.csv.gz',dtype=str,keep_default_na=False)
print(json.dumps({'sicdb_case_columns':list(s.columns),'sicdb_ref_columns':list(ref.columns)}))
rr=ref[ref.ReferenceName.str.contains('sex|gender|admi|death|discharg',case=False,regex=True)]
print(rr.to_json(orient='records'))
