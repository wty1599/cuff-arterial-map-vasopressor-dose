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

"""Extract raw hourly BP blocks with a full source item census; no patient output."""
from pathlib import Path
import csv,gzip,json,time
from collections import Counter
import pandas as pd
ROOT=Path('@CUFF_WORKSPACE@');OUT=ROOT/'phase0';R=OUT/'_restricted';SRC=Path('@SICDB_ROOT@')
with open(OUT/'map_bp_sicdb.csv',encoding='utf-8-sig',newline='') as f:ids={int(x['itemid']) for x in csv.DictReader(f) if x['selected']=='1'}
n=0;k=0;counts=Counter();t0=time.monotonic()
with gzip.open(R/'sicdb_bp_hourly.csv.gz','wt',encoding='utf-8',newline='',compresslevel=1) as f:
    for c in pd.read_csv(SRC/'data_float_h.csv.gz',dtype={'rawdata':str},keep_default_na=False,chunksize=300000):
        counts.update(c['DataID'].value_counts().to_dict());n+=len(c);s=c[c['DataID'].isin(ids)]
        s.to_csv(f,index=False,header=(k==0));k+=len(s)
        if n%3000000==0:(OUT/'_logs/sicdb_scan_progress.json').write_text(json.dumps({'rows_scanned':n,'bp_hourly_rows':k,'elapsed_seconds':round(time.monotonic()-t0)}))
with open(OUT/'sicdb_numeric_item_census.csv','w',encoding='utf-8-sig',newline='') as f:
    w=csv.writer(f);w.writerow(['itemid','rows']);w.writerows(sorted(counts.items()))
(OUT/'_logs/sicdb_scan_complete.json').write_text(json.dumps({'rows_scanned':n,'bp_hourly_rows':k,'elapsed_seconds':round(time.monotonic()-t0),'gzip_to_eof':True}))
print(json.dumps({'source':'sicdb_numeric','rows_scanned':n,'bp_hourly_rows':k}))
