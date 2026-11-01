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

import sys,csv,gzip,json,re
from pathlib import Path
from collections import Counter

import duckdb,pandas as pd
ROOT=Path('@CUFF_WORKSPACE@');OUT=ROOT/'phase0';R=OUT/'_restricted';SRC=Path('@SICDB_ROOT@')
con=duckdb.connect();con.execute("SET memory_limit='1GB'; SET threads=2;")
def aggregate(name,q):
    d=con.execute(q).df();d.to_csv(OUT/(name+'.csv'),index=False);print(json.dumps({'aggregate':name,'rows':len(d)}))
aggregate('mimic_drug_source_units',f"SELECT itemid,rateuom,ordercategorydescription,statusdescription,count(*) n_rows,count(rate) rate_nonmissing FROM read_csv_auto('{R.as_posix()}/mimic_drugs.csv',sample_size=-1) GROUP BY ALL")
aggregate('sicdb_drug_source_units',f"SELECT DrugID,GivenState,IsSingleDose,count(*) n_rows,count(AmountPerMinute) rate_nonmissing FROM read_csv_auto('{SRC.as_posix()}/medication.csv.gz') GROUP BY ALL")
for db in ['amsterdam']:
    maps=pd.read_csv(OUT/f'map_vaso_{db}.csv',dtype=str)
    ids=set(maps['itemid'].astype(int))
    counter=Counter();kept=0
    with gzip.open(R/'amsterdam_drugs.csv.gz','wt',encoding='utf-8',newline='',compresslevel=1) as out:
        for c in pd.read_csv('@AMSTERDAM_ROOT@/drugitems.csv',encoding='cp1252',keep_default_na=False,chunksize=300000):
            counter.update({tuple(k):int(v) for k,v in c.groupby(['itemid','item','ordercategory']).size().items()})
            s=c[c['itemid'].isin(ids)];s.to_csv(out,index=False,header=kept==0);kept+=len(s)
    with open(OUT/'amsterdam_drug_item_census.csv','w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(['itemid','label','ordercategory','rows']);w.writerows((*k,v) for k,v in sorted(counter.items()))
aggregate('amsterdam_drug_source_units',f"SELECT itemid,ordercategory,iscontinuous,doseunit,doserateunit,doserateperkg,rateunit,ratetimeunitid,count(*) n_rows FROM read_csv_auto('{R.as_posix()}/amsterdam_drugs.csv.gz',all_varchar=true) GROUP BY ALL")
a=pd.read_csv('@AMSTERDAM_ROOT@/admissions.csv',encoding='cp1252',keep_default_na=False);a.to_csv(R/'amsterdam_stays.csv',index=False)
for col in ['agegroup','weightgroup','heightgroup','location']:
    a[col].value_counts(dropna=False).rename_axis(col).to_csv(OUT/f'amsterdam_categories_{col}.csv')
d=pd.read_csv(OUT/'amsterdam_dictionary_all.csv',dtype=str).fillna('');rid=set(d.loc[d['item'].str.fullmatch('Hartritme',case=False),'itemid'].astype(int))
kept=0;cnt=Counter()
with gzip.open(R/'amsterdam_rhythm.csv.gz','wt',encoding='utf-8',newline='',compresslevel=1) as out:
    for c in pd.read_csv('@AMSTERDAM_ROOT@/listitems.csv',encoding='cp1252',keep_default_na=False,chunksize=400000):
        cnt.update({tuple(k):int(v) for k,v in c.groupby(['itemid','item']).size().items()})
        s=c[c['itemid'].isin(rid)];s.to_csv(out,index=False,header=kept==0);kept+=len(s)
with open(OUT/'amsterdam_list_item_census.csv','w',encoding='utf-8-sig',newline='') as f:
    w=csv.writer(f);w.writerow(['itemid','label','rows']);w.writerows((*k,v) for k,v in sorted(cnt.items()))
aggregate('amsterdam_rhythm_values',f"SELECT value,count(*) n_rows FROM read_csv_auto('{R.as_posix()}/amsterdam_rhythm.csv.gz',all_varchar=true) GROUP BY ALL")
print('Source-unit probes and supplementary source extracts complete.')
