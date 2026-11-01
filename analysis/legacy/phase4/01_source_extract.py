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
import sys,csv,json,gzip,re,subprocess,datetime
ROOT=Path('@CUFF_WORKSPACE@');OUT=ROOT/'phase4';RES=OUT/'_restricted';CODE=ROOT/'sql/phase4'

import pandas as pd
import duckdb
def save(d,n):d.to_csv(OUT/n,index=False,encoding='utf-8-sig')
def run(db,sql):
 p=CODE/f'source_{db}.sql';p.write_text(sql,encoding='utf-8');cmd=['psql','-X','-v','ON_ERROR_STOP=1','-d','mimiciv' if db=='mimic' else 'eicu','-f','@CUFF_WORKSPACE@/sql/phase4/'+p.name]
 q=subprocess.run(cmd,capture_output=True);(RES/f'source_{db}.log').write_bytes(q.stdout+b'\n'+q.stderr)
 if q.returncode:raise RuntimeError(f'{db}: source query failed, see restricted log')
 print(db,'SOURCE_EXTRACT_COMPLETE',flush=True)
def main():
 assert (OUT/'phase4_params.yaml').exists()
 sql="\\set ON_ERROR_STOP on\nSET timezone='UTC';\nCREATE TEMP TABLE selected_stays(sid bigint);\n\\copy selected_stays FROM '@CUFF_WORKSPACE@/phase4/_restricted/site_stays_mimic.csv' CSV HEADER\nBEGIN READ ONLY;\n"
 sql+="\\copy (SELECT p.stay_id AS sid,extract(epoch FROM p.starttime) AS start,extract(epoch FROM p.endtime) AS stop,p.itemid,p.location,p.locationcategory,p.orderid,p.linkorderid,p.statusdescription FROM mimiciv_icu.procedureevents p JOIN selected_stays s ON p.stay_id=s.sid JOIN mimiciv_icu.d_items d USING(itemid) WHERE d.label='Arterial Line' AND d.linksto='procedureevents') TO '@CUFF_WORKSPACE@/phase4/_restricted/site_source_mimic.csv' CSV HEADER\n"
 sql+="\\copy (SELECT location,locationcategory,count(*) AS records FROM mimiciv_icu.procedureevents p JOIN mimiciv_icu.d_items d USING(itemid) WHERE d.label='Arterial Line' AND d.linksto='procedureevents' GROUP BY 1,2 ORDER BY 1,2) TO '@CUFF_WORKSPACE@/phase4/site_location_dictionary_mimic.csv' CSV HEADER\n"
 sql+="\\copy (SELECT d.itemid,d.label,c.value,count(*) AS records FROM mimiciv_icu.chartevents c JOIN selected_stays s ON c.stay_id=s.sid JOIN mimiciv_icu.d_items d USING(itemid) WHERE d.label ILIKE '%arterial%site%' GROUP BY 1,2,3 ORDER BY 1,3) TO '@CUFF_WORKSPACE@/phase4/site_chartevents_dictionary_mimic.csv' CSV HEADER\nROLLBACK;\n"
 if '--files-only' not in sys.argv:run('mimic',sql)
 sql="\\set ON_ERROR_STOP on\nBEGIN READ ONLY;\n"
 for table in ['nursecare','careplangeneral']:
  fields=['cellattributepath','celllabel','cellattribute'] if table=='nursecare' else ['cplgroup','cplitemvalue']
  concat="concat_ws(' ',"+','.join(fields)+")";query=f"SELECT {','.join(fields)},count(*) AS records FROM public.{table} WHERE {concat} ~* 'arteri|radial|femoral|brachial|axillar|catheter|line' GROUP BY "+','.join(str(i+1) for i in range(len(fields)))
  sql+=f"\\copy ({query}) TO '@CUFF_WORKSPACE@/phase4/site_{table}_dictionary_eicu.csv' CSV HEADER\n"
 sql+="ROLLBACK;\n"
 if '--files-only' not in sys.argv:run('eicu',sql)
 # Amsterdam processitems already includes source-relative millisecond start/stop.
 ams=Path('@AMSTERDAM_ROOT@/processitems.csv');ids=set(pd.read_csv(RES/'site_stays_amsterdam.csv').sid.astype(str));parts=[];census=[]
 for d in pd.read_csv(ams,encoding='cp1252',chunksize=100000):
  z=d[d.item.str.contains('Arterielijn',case=False,na=False)].copy();census.append(z.groupby(['itemid','item']).size().reset_index(name='records'))
  z=z[z.admissionid.astype(str).isin(ids)];parts.append(z)
 pd.concat(parts,ignore_index=True).to_csv(RES/'site_source_amsterdam.csv',index=False)
 save(pd.concat(census).groupby(['itemid','item'],as_index=False).records.sum(),'site_process_dictionary_amsterdam.csv')
 # Check actual SICdb payload for anatomical information without exposing free-text values.
 src=Path('@SICDB_ROOT@');ids=set(pd.read_csv(RES/'site_stays_sicdb.csv').sid.astype(str));parts=[];keys={};vals={};total=0
 for d in pd.read_csv(src/'data_range.csv.gz',chunksize=200000,keep_default_na=False):
  a=d[d.DataID.eq(3045)];total+=len(a)
  for value in a.Data.astype(str).unique():
   kind='empty' if value.strip().lower() in ['','null','none'] else 'has_anatomical_keyword' if re.search(r'radial|femor|brach|axill|dorsalis|ulnar',value,re.I) else 'nonanatomical_payload'
   vals[kind]=vals.get(kind,0)+int(a.Data.astype(str).eq(value).sum())
   try:
    obj=json.loads(value)
    if isinstance(obj,dict):
     for key in obj:keys[key]=keys.get(key,0)+1
   except (ValueError,TypeError):pass
  parts.append(a[a.CaseID.astype(str).isin(ids)])
 pd.concat(parts,ignore_index=True).to_csv(RES/'site_source_sicdb.csv',index=False)
 (OUT/'site_payload_sicdb.json').write_text(json.dumps(dict(all_arterial_line_records=total,payload_classes=vals,json_keys=keys),indent=2),encoding='utf-8')
 print('amsterdam and sicdb source inspection complete; no S3/S4 endpoints computed',flush=True)
if __name__=='__main__':main()
