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
"""Independent source reconstruction: no production analysis imports or writes.
Samples and source records remain under phase6/_restricted. Reconstructs each
window from actual PostgreSQL blood-pressure rows, retaining inherited episode
bounds because treatment eligibility is unchanged in this revision.
"""
from pathlib import Path
import csv,sys,json,hashlib,datetime,subprocess,math,statistics,collections
R=Path('@CUFF_WORKSPACE@'); A=R/'phase6'; Q=A/'_restricted'

import pandas as pd
import duckdb

def write(p,rows):
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def num(x):
    try:
        x=float(x);return x if math.isfinite(x) else None
    except (ValueError,TypeError):return None
def med(x):
    x=[v for v in x if v is not None];return statistics.median(x) if x else None
def wsl(p):
    return Path(p).as_posix()
def q(x):return "'"+str(x).replace("'","''")+"'"
def prepare(db):
    p=Q/f'audit_source_sample_{db}.csv'
    if p.exists():return read(p)
    b=pd.read_parquet(R/'phase1/_restricted'/f'pairs_{db}.parquet')
    w=pd.read_parquet(R/'phase1/_restricted'/f'window2_{db}.parquet')[['cid','art_map']].rename(columns={'art_map':'art_2'})
    s=pd.read_parquet(R/'phase1/_restricted'/f'same_minute_{db}.parquet')[['cid','art_map']].rename(columns={'art_map':'art_same'})
    b=b.merge(w,on='cid',validate='one_to_one').merge(s,on='cid',validate='one_to_one')
    c=duckdb.connect(str(R/'phase0/_restricted'/f'{db}.duckdb'),read_only=True)
    ep=c.execute('SELECT cid,episode_start,episode_stop FROM cuff60').df();c.close()
    b=b.merge(ep,on='cid',validate='one_to_one')
    b['key']=[hashlib.sha256(f'20261027|phase6_source|{db}|{v}'.encode()).hexdigest() for v in b.cid]
    b['changed']=(b.art_map!=b.art_same)|(b.art_map!=b.art_2)
    chosen=[]
    for st in ['<=0.1','>0.1-0.25','>0.25-0.5','>0.5-1.0','>1.0']:
        part=b[b.dose_stratum.eq(st)].sort_values('key')
        both=pd.concat([part[part.changed].head(3),part[~part.changed].head(3)])
        if len(both)<6:both=pd.concat([both,part[~part.cid.isin(both.cid)].head(6-len(both))])
        chosen.append(both)
    z=pd.concat(chosen).reset_index(drop=True)
    z.insert(0,'audit_id',[f'{db}_phase6_{i:02d}' for i in range(1,len(z)+1)])
    rows=z[['audit_id','sid','t','episode_start','episode_stop','cuff_map','art_map','art_2','art_same','dose_stratum','changed']].to_dict('records')
    write(p,rows);return rows

def query(db,sample):
    path=R/'sql/phase6'/f'audit_source_{db}.sql'
    sample_path=Q/f'audit_source_sample_{db}.csv'
    cols=list(sample[0])
    sql='\\set ON_ERROR_STOP on\nSET client_min_messages=warning;\nSET TIME ZONE \'UTC\';\n'
    sql+='CREATE TEMP TABLE audit_sample ('+','.join(f'"{x}" text' for x in cols)+');\n'
    sql+=f'\\copy audit_sample FROM {q(wsl(sample_path))} CSV HEADER\n'
    if db=='mimic':
        mp=R/'phase0/map_bp_mimic.csv';header=list(read(mp)[0])
        sql+='CREATE TEMP TABLE audit_map ('+','.join(f'"{x}" text' for x in header)+');\n'
        sql+=f'\\copy audit_map FROM {q(wsl(mp))} CSV HEADER\n'
        query="""SELECT a.audit_id,extract(epoch FROM c.charttime) source_t,m.classification modality,
          m.family,m.component,c.valuenum value,c.itemid,c.valueuom
          FROM audit_sample a JOIN mimiciv_icu.chartevents c ON c.stay_id=a.sid::bigint
          JOIN audit_map m ON m.itemid::bigint=c.itemid AND m.selected='1'
          WHERE c.charttime BETWEEN timestamp '1970-01-01'+(a.t::double precision-300)*interval '1 second'
             AND timestamp '1970-01-01'+(a.t::double precision+300)*interval '1 second'
          AND (m.classification='arterial' OR (m.classification='cuff' AND
            c.charttime=timestamp '1970-01-01'+a.t::double precision*interval '1 second'))"""
    else:
        queries=[]
        for table,mod,prefix in [('vitalperiodic','arterial','systemic'),('vitalaperiodic','cuff','noninvasive')]:
            filt='v.observationoffset::double precision*60 BETWEEN a.t::double precision-300 AND a.t::double precision+300' if mod=='arterial' else 'v.observationoffset::double precision*60=a.t::double precision'
            for comp,col in [('map','mean'),('sbp','systolic'),('dbp','diastolic')]:
                queries.append(f"SELECT a.audit_id,v.observationoffset::double precision*60 source_t,'{mod}' modality,'{mod}' family,'{comp}' component,v.{prefix+col}::double precision value,'{table}.{prefix+col}' itemid,'mmHg' valueuom FROM audit_sample a JOIN public.{table} v ON v.patientunitstayid=a.sid::bigint WHERE {filt}")
        query=' UNION ALL '.join(queries)
    sql+='BEGIN READ ONLY;\nSET LOCAL statement_timeout=0;\n'
    sql+=f'\\copy ({query.replace(chr(10)," ")}) TO {q(wsl(Q/f"audit_source_rows_{db}.csv"))} CSV HEADER\nCOMMIT;\n'
    path.write_text(sql,encoding='utf-8')
    source=Q/f'audit_source_rows_{db}.csv'
    if not source.exists():
        res=subprocess.run([_cuff_env.environ.get('CUFF_PSQL','psql'),'-X','-v','ON_ERROR_STOP=1','-d','mimiciv' if db=='mimic' else 'eicu','-f',wsl(path)],capture_output=True,text=True,encoding='utf-8',errors='replace')
        (Q/f'audit_source_query_{db}.log').write_text(res.stdout+'\n'+res.stderr,encoding='utf-8')
        if res.returncode:raise RuntimeError(f'{db} source query failed, details restricted')
    return read(source),path

def compare(db,sample,raw):
    component=collections.defaultdict(list)
    for x in raw:component[(x['audit_id'],float(x['source_t']),x['modality'],x['family'],x['component'])].append(num(x['value']))
    family=collections.defaultdict(dict)
    for (*k,comp),values in component.items():family[tuple(k)][comp]=med(values)
    valid=collections.defaultdict(list);bad=0
    for (aid,t,mod,fam),v in family.items():
        m,s,d=(v.get(x) for x in ['map','sbp','dbp'])
        ok=m is not None and 20<=m<=200 and (s is None or s<=300) and (s is None or d is None or (s>d and s-d>=10))
        if ok:valid[(aid,t,mod)].append(m)
        else:bad+=1
    times=collections.defaultdict(dict)
    for (aid,t,mod),v in valid.items():times[(aid,mod)][t]=med(v)
    out=[]
    for x in sample:
        aid=x['audit_id'];t=float(x['t']);lo=float(x['episode_start']);hi=float(x['episode_stop'])
        cf=times[(aid,'cuff')].get(t)
        art={s:m for s,m in times[(aid,'arterial')].items() if lo<=s<hi}
        for name,key,pred in [('5min','art_map',lambda s:abs(s-t)<=300),('2min','art_2',lambda s:abs(s-t)<=120),('same_minute','art_same',lambda s:math.floor(s/60)==math.floor(t/60))]:
            vals=[m for s,m in art.items() if pred(s)];am=med(vals);expected=float(x[key])
            out.append(dict(audit_id=aid,window=name,expected_cuff=float(x['cuff_map']),source_cuff=cf,expected_art=expected,source_art=am,source_values=len(vals),cuff_matches=cf is not None and abs(cf-float(x['cuff_map']))<1e-9,arterial_matches=am is not None and abs(am-expected)<1e-9))
    write(Q/f'audit_source_comparison_{db}.csv',out)
    return dict(database=db,sampled_cuff_anchors=len(sample),source_records=len(raw),window_comparisons=len(out),cuff_mismatches=sum(not x['cuff_matches'] for x in out),arterial_mismatches=sum(not x['arterial_matches'] for x in out),source_family_rows_rejected=bad,passed=all(x['cuff_matches'] and x['arterial_matches'] for x in out),scope='Actual PostgreSQL BP source rows; retained original treatment-episode bounds; no waveform adjudication')

if __name__=='__main__':
    summaries=[]
    for db in ['mimic','eicu']:
        sample=prepare(db);raw,path=query(db,sample);row=compare(db,sample,raw);summaries.append(row)
        print(json.dumps(row),flush=True)
    write(A/'phase6_source_audit_summary.csv',summaries)
    receipt={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'parameter_sha256':hashlib.sha256((A/'phase6_params.yaml').read_bytes()).hexdigest(),'sample_design':'6 known-dose common anchors per stratum; SHA256 ordering, target 3 with and 3 without a reference-window change; no outcome-directed selection','summaries':summaries}
    (A/'phase6_source_audit_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
