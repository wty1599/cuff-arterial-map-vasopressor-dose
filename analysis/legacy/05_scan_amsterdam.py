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

"""One bounded-memory pass over original numeric rows, with item census."""
import csv,json,zipfile,gzip,time,re
from pathlib import Path
from collections import Counter
import pandas as pd
ROOT=Path('@CUFF_WORKSPACE@');OUT=ROOT/'phase0';R=OUT/'_restricted'
with open(OUT/'map_bp_amsterdam.csv',encoding='utf-8-sig',newline='') as f:
    bpids={int(x['itemid']) for x in csv.DictReader(f) if x['selected']=='1'}
dictionary=pd.read_csv(OUT/'amsterdam_dictionary_all.csv',encoding='utf-8-sig',dtype=str).fillna('')
extra=dictionary[(dictionary['item']+' '+dictionary['item_en']).str.contains(r'gewicht|weight|lengte|height|creatin|kreatinin',case=False,regex=True)]
extra[['itemid','item','item_en','unit','category']].drop_duplicates().to_csv(OUT/'amsterdam_covariate_lab_candidates.csv',index=False)
ids=bpids|set(extra['itemid'].astype(int))
counts=Counter();n=0;kept=0;t0=time.monotonic()
with zipfile.ZipFile(R/'amsterdam_numeric_inner.zip') as z,gzip.open(R/'amsterdam_numeric_selected.csv.gz','wt',encoding='utf-8',newline='',compresslevel=1) as out:
    assert z.namelist()==['numericitems.csv']
    with z.open('numericitems.csv') as f:
        for chunk in pd.read_csv(f,encoding='cp1252',usecols=['admissionid','itemid','item','value','unit','measuredat'],dtype={'admissionid':'int64','itemid':'int64','item':'str','value':'str','unit':'str','measuredat':'int64'},keep_default_na=False,chunksize=500000):
            n+=len(chunk)
            counts.update({tuple(k):int(v) for k,v in chunk.groupby(['itemid','item','unit'],dropna=False).size().items()})
            s=chunk[chunk['itemid'].isin(ids)]
            s.to_csv(out,index=False,header=(kept==0));kept+=len(s)
            if n%10000000==0:
                (OUT/'_logs/amsterdam_scan_progress.json').write_text(json.dumps({'rows_scanned':n,'rows_retained':kept,'elapsed_seconds':round(time.monotonic()-t0)}))
    # Reading to EOF also verifies the ZIP member CRC.
with open(OUT/'amsterdam_numeric_item_census.csv','w',encoding='utf-8-sig',newline='') as f:
    w=csv.writer(f);w.writerow(['itemid','label','unitname','rows']);w.writerows((*k,v) for k,v in sorted(counts.items()))
(OUT/'_logs/amsterdam_scan_complete.json').write_text(json.dumps({'rows_scanned':n,'rows_retained':kept,'elapsed_seconds':round(time.monotonic()-t0),'crc_to_eof':True}))
print(json.dumps({'source':'amsterdam_numeric','rows_scanned':n,'rows_retained':kept}))
