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
"""Bounded, independent source verification of preselected missing-MAP cases.

Reads only the presampled missing-check input and real database source tables.
Patient-target evidence remains restricted; only target-level aggregate CSVs
are written to phase0. No production pair, blood-pressure, or drug table is read.
"""
import argparse
import collections
import csv
import datetime as dt
import hashlib
import json
import re
import subprocess
from pathlib import Path

from audit_source_recheck import ROOT,SQLDIR,OUT,RESTRICTED,read_csv,number,sid,sql_literal,wsl_path,temp_import


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# This is an independent interpretation of all source d_items labels, without
# importing the production item map or embedding any item IDs.
CONCEPTS = r"""
SELECT itemid,label,abbreviation,category,unitname,
 CASE
 WHEN lower(label) ~ '(non[ -]?invasive|non invasive|\mnibp\M|\mnbp\M)'
  AND lower(label) ~ '(mean|\mmap\M)' THEN 'cuff'
 WHEN lower(label) ~ '(arterial.*(blood pressure|pressure)|\mart[ .-]*bp\M|\mabp\M)'
  AND lower(label) ~ '(mean|\mmap\M)' THEN 'arterial'
 ELSE NULL END AS target
FROM mimiciv_icu.d_items
WHERE linksto='chartevents'
 AND lower(label) !~ '(alarm|target|iabp|pulmonary|\mpap\M|streef)'
"""


def build_sql(db):
    inp=RESTRICTED/f'missing_check_input_{db}.csv'
    txt='\\set ON_ERROR_STOP on\nSET client_min_messages=warning;\nSET TIME ZONE \'UTC\';\n'
    txt+=temp_import('audit_missing_sample',inp,['sid','target','episode_start','episode_stop'])
    txt+='ANALYZE audit_missing_sample;\nBEGIN READ ONLY;\nSET LOCAL statement_timeout=0;\n'
    output=RESTRICTED/f'audit_missing_source_{db}_patient_targets.csv'
    base='''intervals AS (SELECT DISTINCT sid::bigint AS sid,target,episode_start::double precision AS a,
                episode_stop::double precision AS b FROM audit_missing_sample),
            pt AS (SELECT sid,target,count(*) AS episode_rows FROM intervals GROUP BY sid,target)'''
    if db=='mimic':
        dictionary=RESTRICTED/'audit_missing_source_mimic_dictionary.csv'
        q=f'SELECT * FROM ({CONCEPTS}) d WHERE target IS NOT NULL ORDER BY target,itemid'
        txt+=f"\\copy ({q.replace(chr(10),' ')}) TO {sql_literal(wsl_path(dictionary))} CSV HEADER\n"
        q=f'''WITH {base}, concepts AS ({CONCEPTS}), observations AS (
            SELECT c.stay_id AS sid,d.target,
                count(*) AS n_mapped_source_rows,
                count(c.valuenum) AS n_source_numeric,
                count(*) FILTER (WHERE c.valuenum BETWEEN 20 AND 200) AS n_source_in_range,
                count(*) FILTER (WHERE nullif(btrim(c.value),'') IS NOT NULL) AS n_source_nonblank_text
            FROM mimiciv_icu.chartevents c JOIN concepts d ON c.itemid=d.itemid
            JOIN pt ON c.stay_id=pt.sid AND d.target=pt.target
            WHERE d.target IS NOT NULL AND EXISTS (
                SELECT 1 FROM intervals i WHERE i.sid=c.stay_id AND i.target=d.target
                AND c.charttime>=timestamp '1970-01-01'+i.a*interval '1 second'
                AND c.charttime<timestamp '1970-01-01'+i.b*interval '1 second')
            GROUP BY c.stay_id,d.target)
            SELECT pt.sid,pt.target,pt.episode_rows,
                coalesce(o.n_mapped_source_rows,0) AS n_mapped_source_rows,
                coalesce(o.n_source_numeric,0) AS n_source_numeric,
                coalesce(o.n_source_in_range,0) AS n_source_in_range,
                coalesce(o.n_source_nonblank_text,0) AS n_source_nonblank_text
            FROM pt LEFT JOIN observations o USING(sid,target) ORDER BY pt.target,pt.sid'''
    else:
        pieces=[]
        for target,table,column in [('arterial','vitalperiodic','systemicmean'),('cuff','vitalaperiodic','noninvasivemean')]:
            pieces.append(f'''SELECT v.patientunitstayid AS sid,'{target}' AS target,
                count(*) AS n_mapped_source_rows,count(v.{column}) AS n_source_numeric,
                count(*) FILTER (WHERE v.{column} BETWEEN 20 AND 200) AS n_source_in_range,
                count(v.{column}) AS n_source_nonblank_text
                FROM public.{table} v JOIN pt ON v.patientunitstayid=pt.sid AND pt.target='{target}'
                WHERE EXISTS (SELECT 1 FROM intervals i WHERE i.sid=v.patientunitstayid AND i.target='{target}'
                    AND v.observationoffset::double precision*60>=i.a
                    AND v.observationoffset::double precision*60<i.b)
                GROUP BY v.patientunitstayid''')
        q=f'''WITH {base}, observations AS ({' UNION ALL '.join(pieces)})
            SELECT pt.sid,pt.target,pt.episode_rows,
                coalesce(o.n_mapped_source_rows,0) AS n_mapped_source_rows,
                coalesce(o.n_source_numeric,0) AS n_source_numeric,
                coalesce(o.n_source_in_range,0) AS n_source_in_range,
                coalesce(o.n_source_nonblank_text,0) AS n_source_nonblank_text
            FROM pt LEFT JOIN observations o USING(sid,target) ORDER BY pt.target,pt.sid'''
    txt+=f"\\copy ({q.replace(chr(10),' ')}) TO {sql_literal(wsl_path(output))} CSV HEADER\nCOMMIT;\n"
    path=SQLDIR/f'audit_missing_source_{db}.sql'
    path.write_text(txt,encoding='utf-8')
    return path,inp,output


def execute(db):
    sql,inp,source_output=build_sql(db)
    sample=read_csv(inp); seen=collections.defaultdict(set)
    for r in sample:
        if r['target'] not in {'arterial','cuff'}:
            raise ValueError('Unexpected target')
        if number(r['episode_start']) is None or number(r['episode_stop']) is None or float(r['episode_stop'])<=float(r['episode_start']):
            raise ValueError('Invalid preselected episode interval')
        seen[r['target']].add(sid(r['sid']))
    if set(seen)!={'arterial','cuff'} or any(len(x)!=50 for x in seen.values()):
        raise ValueError('Input must contain exactly 50 preselected patient-target combinations per target')
    input_hash=digest(inp); start=dt.datetime.now(dt.timezone.utc).isoformat()
    database='mimiciv' if db=='mimic' else 'eicu'
    result=subprocess.run([_cuff_env.environ.get('CUFF_PSQL','psql'),'-X',
                           '-v','ON_ERROR_STOP=1','-d',database,'-f',wsl_path(sql)],
                           capture_output=True,text=True,encoding='utf-8',errors='replace')
    (RESTRICTED/f'audit_missing_source_{db}.log').write_text(result.stdout+'\n'+result.stderr,encoding='utf-8')
    if result.returncode:
        raise RuntimeError('Source SQL failed; diagnostics retained in restricted log')
    if input_hash!=digest(inp):
        raise RuntimeError('Input changed during source query')
    rows=read_csv(source_output)
    actual={(sid(r['sid']),r['target']) for r in rows}
    expected={(s,t) for t,ss in seen.items() for s in ss}
    if len(rows)!=len(expected) or actual!=expected:
        raise RuntimeError('Source-returned patient-target combinations do not conserve the presample')
    ended=dt.datetime.now(dt.timezone.utc).isoformat(); aggregates=[]
    for target in ['arterial','cuff']:
        selected=[r for r in rows if r['target']==target]
        with_numeric=sum(int(r['n_source_numeric'])>0 for r in selected)
        in_range=sum(int(r['n_source_in_range'])>0 for r in selected)
        nonblank=sum(int(r['n_source_nonblank_text'])>0 for r in selected)
        aggregates.append({'database':db,'target':target,'patient_target_combinations':len(selected),
            'with_source_numeric_map':with_numeric,'without_source_numeric_map':len(selected)-with_numeric,
            'with_source_map_20_to_200':in_range,'with_nonblank_source_value':nonblank,
            'episode_rows_checked':sum(int(r['episode_rows']) for r in selected),
            'source_query_started_at':start,'source_query_completed_at':ended,
            'timing':'初步 R3 结果产生后，对既有缺失红旗进行定向源回查',
            'scope':'预先抽取的缺失患者-目标组合；未重新定义队列、治疗期或主血压来源',
            'interpretation':'存在源数值不自动证明可配对或通过全部伪迹规则；无源值限于既定治疗区间与主来源'})
    dest=OUT/f'audit_missing_source_{db}_aggregate.csv'
    with dest.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(aggregates[0]));w.writeheader();w.writerows(aggregates)
    manifest={'database':db,'input_sha256':input_hash,'sql_sha256':digest(sql),
              'audit_script_sha256':digest(Path(__file__)),'source_output_sha256':digest(source_output),
              'aggregate_sha256':digest(dest),'source_query_started_at':start,'source_query_completed_at':ended,
              'patient_target_combinations_conserved':True,
              'timing':'初步 R3 结果之后；不依据核查结果改变预设定义'}
    (RESTRICTED/f'audit_missing_source_{db}_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    return {'database':db,'aggregate_file':str(dest),'results':aggregates,
            'patient_target_combinations_conserved':True,'source_queries_executed':True}


def main():
    p=argparse.ArgumentParser();p.add_argument('--db',choices=['mimic','eicu'],required=True)
    p.add_argument('--prepare-only',action='store_true');a=p.parse_args()
    try:
        if a.prepare_only:
            path,_,_=build_sql(a.db)
            output={'database':a.db,'sql_path':str(path),'source_queries_executed':False}
        else:
            output=execute(a.db)
        print(json.dumps(output,ensure_ascii=False))
    except Exception as e:
        import traceback
        (RESTRICTED/f'audit_missing_source_{a.db}_failure.log').write_text(traceback.format_exc(),encoding='utf-8')
        print(json.dumps({'database':a.db,'status':'failed','error_type':type(e).__name__,
                          'diagnostics':'Restricted audit log; no patient rows emitted'},ensure_ascii=False))
        raise SystemExit(1)


if __name__=='__main__':
    main()
