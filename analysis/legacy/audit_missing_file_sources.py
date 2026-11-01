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
"""Bounded negative-record source checks, independent of production BP caches."""
from pathlib import Path
import sys,csv,json,math,struct,re

import duckdb,pandas as pd
ROOT=Path('@CUFF_WORKSPACE@');O=ROOT/'phase0';R=O/'_restricted';S=Path('@SICDB_ROOT@');db=sys.argv[1]
c=duckdb.connect(str(R/f'audit_missing_file_{db}.duckdb'));c.execute("SET memory_limit='512MB';SET threads=2;")
def rd(p,latin=False):return f"read_csv_auto('{Path(p).as_posix()}',all_varchar=true"+(",encoding='latin-1'" if latin else '')+")"
c.execute(f"CREATE OR REPLACE TABLE input AS SELECT sid,target,CAST(episode_start AS DOUBLE) lo,CAST(episode_stop AS DOUBLE) hi FROM {rd(R/f'missing_check_input_{db}.csv')}")
if db=='amsterdam':
    q=f"""CREATE OR REPLACE TABLE found AS WITH dictionary AS (
      SELECT DISTINCT itemid,CASE WHEN lower(item_en) LIKE 'non-invasive%' OR lower(item) LIKE 'niet invasieve%' THEN 'cuff' ELSE 'arterial' END AS target
      FROM {rd(O/'amsterdam_dictionary_all.csv')} WHERE lower(item_en) LIKE '%mean arterial pressure%' OR lower(item) LIKE 'abp gemiddeld%' OR lower(item)='picco apm')
      SELECT DISTINCT i.sid,i.target FROM {rd(R/'amsterdam_numeric_source.csv',True)} x JOIN dictionary d ON x.itemid=d.itemid JOIN input i ON x.admissionid=i.sid AND d.target=i.target AND TRY_CAST(x.measuredat AS DOUBLE)/1000>=i.lo AND TRY_CAST(x.measuredat AS DOUBLE)/1000<i.hi WHERE TRY_CAST(x.value AS DOUBLE) IS NOT NULL"""
    c.execute(q);unknown={}
else:
    q=f"""WITH dictionary AS (
      SELECT DISTINCT ReferenceGlobalID AS itemid,CASE WHEN ReferenceValue LIKE 'BloodPressureNI%' THEN 'cuff' ELSE 'arterial' END AS target
      FROM {rd(S/'d_references.csv.gz')} WHERE ReferenceName='SignalFloat' AND regexp_matches(ReferenceValue,'^BloodPressure(Arterial|NI)MAP$'))
      SELECT i.sid,i.target,i.lo,i.hi,x.Offset,x.rawdata,x.cnt FROM {rd(S/'data_float_h.csv.gz')} x JOIN dictionary d ON x.DataID=d.itemid JOIN input i ON x.CaseID=i.sid AND d.target=i.target AND TRY_CAST(x.Offset AS DOUBLE)<i.hi AND TRY_CAST(x.Offset AS DOUBLE)+3600>i.lo"""
    blocks=c.execute(q).fetchall();found=set();unknown={}
    for sid,target,lo,hi,offset,h,cnt in blocks:
        if not h or len(h)!=482:
            unknown[(sid,target)]=unknown.get((sid,target),0)+1;continue
        data=bytes.fromhex(h[2:]);start=float(offset)
        for j in range(60):
            b=data[j*4:j*4+4]
            if b!=b'\0\0\0\0' and lo<=start+j*60<hi and math.isfinite(struct.unpack('<f',b)[0]):found.add((sid,target));break
    c.execute('CREATE OR REPLACE TABLE found(sid VARCHAR,target VARCHAR)')
    if found:c.executemany('INSERT INTO found VALUES(?,?)',sorted(found))
res=c.execute('SELECT a.sid,a.target,f.sid IS NOT NULL has_source_value FROM (SELECT DISTINCT sid,target FROM input) a LEFT JOIN found f USING(sid,target)').df();res['overlapping_hourly_only_blocks']=[unknown.get((x.sid,x.target),0) for x in res.itertuples()];res.to_csv(R/f'audit_missing_file_comparison_{db}.csv',index=False)
agg=res.groupby('target').agg(sample_patient_target_combinations=('sid','size'),with_source_values=('has_source_value','sum'),with_hourly_only_blocks=('overlapping_hourly_only_blocks',lambda x:int((x>0).sum()))).reset_index();agg['without_source_values']=agg.sample_patient_target_combinations-agg.with_source_values;agg['timing']='after preliminary missingness and R3';agg.to_csv(O/f'audit_missing_source_{db}_aggregate.csv',index=False)
(ROOT/'sql'/f'audit_missing_{db}_source_query.sql').write_text(q+';\n',encoding='utf-8');print(json.dumps({'db':db,'combinations':len(res),'with_source_values':int(res.has_source_value.sum())}));c.close()
