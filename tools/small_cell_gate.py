#!/usr/bin/env python3
"""Fail closed before delivery, upload, commit or push.

Counts of 1--9 paired with numeric measurements/derived estimates are blocked.
Patient identifiers and clinical timestamp fields are blocked independently.
Diagnostics contain paths, field names and row locations, never cell values.
JSON audit receipts carry operational metadata; field semantics, not the mere
presence of a number, distinguish them from clinical measurement tables.
ZIP payloads are inspected recursively. Restricted folders are never a release
destination and are excluded only during a whole-workspace inventory scan.
"""
from __future__ import annotations
import argparse,csv,datetime,hashlib,io,json,math,re,subprocess,sys,zipfile
from pathlib import Path

VERSION='1.0'
_CACHE={}
DATA={'.csv','.tsv','.json','.xlsx','.xlsm','.parquet','.zip'}
SNIFF={'.txt','.dat','.jsonl','.ndjson'}
RESTRICTED={'restricted','_restricted'}
VENDOR={'.git','.venv','venv','runtime','_vendor','site-packages','node_modules','__pycache__','repo_tools'}
ID_RE=re.compile(r'(^id$|^pid$|^sid$|^cid$|^patient$|^subject$|(^|_)(subject|patient|person|case|stay|hadm|admission|encounter|record|cluster|pair|row|anon)(id|_id|_index|_key|_number)$|patientunitstayid|patienthealthsystemstayid|admissionid|caseid|recordid)',re.I)
TIME_RE=re.compile(r'(^|_)(charttime|storetime|starttime|endtime|admittime|dischtime|intime|outtime|dob|dateofbirth|timestamp|time|date|measurement_time|cuff_time|art_time|arterial_time|pair_time|event_time|time_offset|offset)$',re.I)
MEASURE_RE=re.compile(r'(^|_)(map|bp|sbp|dbp|mean_map|diff|difference|bias|loa|estimate|risk|rate|probability|proportion|percent|percentage|median|mean|minimum|maximum|q1|q3|iqr|value|age|weight|height|bmi|dose|neeq|ne_eq|temperature|heart_rate|lactate|sensitivity|specificity|event_share|denominator_share|original_estimate|estimate_after_exclusion|difference_pp|lr_minus|lr_plus|lr_negative|lr_positive|ci_lower|ci_upper|lcl|ucl|lower_ci|upper_ci)(_|$)',re.I)
COUNT_EXCLUDE=re.compile(r'(bootstrap|resampl|knots|seed|byte|word|reference|citation|file|page|line|script|column|hash|parameter|timeout|version|attempt|replicate|dpi|font|release)',re.I)
def count_field(k):
 k=str(k).lower()
 if COUNT_EXCLUDE.search(k):return False
 return k in {'n','count','counts','frequency','frequency_n','sample_size','pairs','patients','n_pairs','n_patients','all_pairs','all_patients','numerator','denominator','top_patients','excluded_top_event_contributors','n_missing','n_observed'} or k.startswith(('denominator','cell_')) or bool(re.search(r'(^n_|_pairs$|_patients$|_count$|^num_|^denom_|_frequency$)',k))
def safe_field(k):
 k=str(k)
 if re.fullmatch(r'\d+',k) or re.search(r'\d{5,}',k) or len(k)>100:
  return 'key_sha256_'+hashlib.sha256(k.encode()).hexdigest()[:16]
 return k
def number(x):
 if isinstance(x,bool) or x is None:return None
 try:
  n=float(x)
  return n if math.isfinite(n) else None
 except (TypeError,ValueError):return None
def decode(b):
 for enc in ('utf-8-sig','utf-16','gb18030'):
  try:return b.decode(enc)
  except UnicodeError:pass
 raise ValueError('unsupported_text_encoding')
def strict_json(t):
 def unique(pairs):
  obj={}
  for k,v in pairs:
   if k in obj:raise ValueError('duplicate_json_key')
   obj[k]=v
  return obj
 return json.loads(t,object_pairs_hook=unique)
def json_tables(obj,path='$'):
 if isinstance(obj,list):
  if obj and all(isinstance(r,dict) for r in obj):
   yield path,obj
   for i,r in enumerate(obj):
    for j,(k,v) in enumerate(r.items()):
     if isinstance(v,(dict,list)):yield from json_tables(v,path+'/'+str(i)+'/node'+str(j))
  else:
   for i,r in enumerate(obj):
    if isinstance(r,(dict,list)):yield from json_tables(r,path+'/'+str(i))
 elif isinstance(obj,dict):
  if isinstance(obj.get('columns'),list) and isinstance(obj.get('data'),list):
   cols=obj['columns'];rows=obj['data']
   if all(isinstance(r,list) for r in rows):yield path+'/data',[dict(zip(cols,r)) for r in rows]
  simple={k:v for k,v in obj.items() if not isinstance(v,(list,dict))}
  if simple:yield path,[simple]
  vectors={k:v for k,v in obj.items() if isinstance(v,list) and v and all(not isinstance(x,(dict,list)) for x in v)}
  if len(vectors)>1 and len({len(v) for v in vectors.values()})==1:
   yield path+'/vectors',[dict(zip(vectors,v)) for v in zip(*vectors.values())]
  for j,(k,v) in enumerate(obj.items()):
   if isinstance(v,(dict,list)):yield from json_tables(v,path+'/node'+str(j))
def tables(b,suffix):
 if suffix in {'.csv','.tsv'}:
  text=decode(b)
  try:delimiter=csv.Sniffer().sniff(text[:8192],delimiters=',\t;|').delimiter
  except csv.Error:delimiter='\t' if suffix=='.tsv' else ','
  reader=csv.DictReader(io.StringIO(text),delimiter=delimiter);cols=reader.fieldnames or []
  if len(cols)!=len(set(cols)):raise RuntimeError('duplicate_csv_header')
  if cols and all(number(k) is not None for k in cols):raise RuntimeError('unclassified_numeric_header')
  if any(re.search(r'[,\t;|]',str(k)) for k in cols):raise RuntimeError('ambiguous_delimiter_header')
  rows=list(reader)
  if any(None in r or any(v is None for v in r.values()) for r in rows):raise RuntimeError('ragged_csv_rows')
  yield '',rows
 elif suffix=='.json':
  obj=strict_json(decode(b))
  if isinstance(obj,(int,float)) or (isinstance(obj,list) and obj and all(isinstance(v,(int,float,list)) for v in obj)):
   raise RuntimeError('unclassified_numeric_json')
  def check_keys(x):
   if isinstance(x,dict):
    if any(re.fullmatch(r'\d+',str(k)) for k in x):raise RuntimeError('dynamic_numeric_json_keys')
    for v in x.values():
     if isinstance(v,(dict,list)):check_keys(v)
   elif isinstance(x,list):
    for v in x:
     if isinstance(v,(dict,list)):check_keys(v)
  check_keys(obj)
  yield from json_tables(obj)
 elif suffix in {'.xlsx','.xlsm'}:
  try:import openpyxl
  except ImportError:raise RuntimeError('openpyxl_required')
  wb=openpyxl.load_workbook(io.BytesIO(b),data_only=False,read_only=True)
  for sheet in wb:
   it=iter(sheet.values);cols=[str(v or '') for v in next(it,())]
   if len(cols)!=len(set(cols)):raise RuntimeError('duplicate_spreadsheet_header')
   rows=[dict(zip(cols,r)) for r in it]
   if any(isinstance(v,str) and v.startswith('=') for r in rows for v in r.values()):
    raise RuntimeError('spreadsheet_formula_requires_materialized_values')
   yield sheet.title,rows
 elif suffix=='.parquet':
  try:import pyarrow.parquet as pq
  except ImportError:raise RuntimeError('pyarrow_required')
  table=pq.read_table(io.BytesIO(b))
  if len(table.column_names)!=len(set(table.column_names)):raise RuntimeError('duplicate_parquet_column')
  yield '',table.to_pylist()
def scan_bytes(b,name,depth=0):
 key=(hashlib.sha256(b).digest(),Path(name).suffix.lower(),depth if Path(name).suffix.lower()=='.zip' else 0)
 if key in _CACHE:
  f,e,old=_CACHE[key]
  def replace_path(r):return {**r,'path':name+r['path'][len(old):]}
  return [replace_path(r) for r in f],[replace_path(r) for r in e]
 f,e=_scan_bytes(b,name,depth)
 _CACHE[key]=(f,e,name)
 return f,e
def _scan_bytes(b,name,depth=0):
 suffix=Path(name).suffix.lower();findings=[];entries=[]
 if suffix in {'.csv','.tsv'} and b.lstrip(b'\xef\xbb\xbf \t\r\n').startswith((b'{',b'[')):
  try:strict_json(decode(b))
  except ValueError:pass
  else:return [{'path':name,'rule':'FORMAT_MISMATCH','row_location':'','fields':[]}],entries
 if suffix in SNIFF:
  try:
   t=decode(b).strip()
   if t.startswith(('{','[')):
    strict_json(t)
    virtual=name+'.__as_json__.json'
    f,e=scan_bytes(b,virtual,depth)
    return [{**r,'path':name+r['path'][len(virtual):]} for r in f],[{**r,'path':name+r['path'][len(virtual):]} for r in e]
   first=next((v for v in t.splitlines() if v.strip()),'')
   keys=[v.strip() for v in re.split(r'[,\t;|]',first)]
   if len(keys)>1 and any(count_field(k) or ID_RE.search(k) or TIME_RE.search(k) or MEASURE_RE.search(k) for k in keys):
    return [{'path':name,'rule':'UNSUPPORTED_TABULAR_EXTENSION','row_location':'header','fields':[safe_field(k) for k in keys]}],entries
  except (ValueError,UnicodeError):
   return [{'path':name,'rule':'UNPARSED_STRUCTURED_TEXT','row_location':'','fields':[]}],entries
  return findings,entries
 if suffix=='.zip':
  if depth>=6:return [{'path':name,'rule':'ARCHIVE_DEPTH_LIMIT','row_location':'','fields':[]}],entries
  try:
   with zipfile.ZipFile(io.BytesIO(b)) as z:
    for i in z.infolist():
     if i.is_dir():continue
     parts=Path(i.filename.replace('\\','/')).parts
     if any(p in RESTRICTED or p=='.git' for p in parts):
      findings.append({'path':name+'!/'+i.filename,'rule':'RESTRICTED_MEMBER','row_location':'','fields':[]});continue
     if Path(i.filename).suffix.lower() in DATA|SNIFF:
      f,e=scan_bytes(z.read(i),name+'!/'+i.filename,depth+1);findings+=f;entries+=e
  except Exception as ex:findings.append({'path':name,'rule':'READ_ERROR','row_location':'','fields':[],'error_type':type(ex).__name__})
  return findings,entries
 if suffix not in DATA:return findings,entries
 digest=hashlib.sha256(b).hexdigest()
 try:
  for sheet,rows in tables(b,suffix):
   cols=sorted({str(k) for r in rows for k in r if k is not None});counts=[k for k in cols if count_field(k)]
   ids=[k for k in cols if ID_RE.search(k)];times=[k for k in cols if TIME_RE.search(k)]
   # These literal tokens name source schema clocks in a parameter record;
   # they are neither a recorded time nor an elapsed-time value.
   clock_tokens={'charttime','storetime','unitadmittime','observationoffset','chartoffset'}
   times=[k for k in times if not (rows and all(r.get(k) in clock_tokens for r in rows if r.get(k) is not None) and any(r.get(k) is not None for r in rows))]
   times=[k for k in times if not (len(rows)==1 and k.endswith('_observation_time') and isinstance(rows[0].get(k),str) and not re.search(r'\d',rows[0][k]) and any(token in rows[0][k].lower() for token in clock_tokens))]
   clinical=bool(any(MEASURE_RE.search(k) for k in cols)) or bool(counts)
   # Operational timestamps in audit receipts are not patient timestamps. Clinical
   # tables still block generic timestamp/time fields and all specific time aliases.
   if not clinical:times=[k for k in times if k.lower() not in {'timestamp','time','date','event_time'}]
   if ids:findings.append({'path':name,'table':safe_field(sheet),'rule':'ID_COLUMN','row_location':'header','fields':[safe_field(k) for k in ids]})
   if times:findings.append({'path':name,'table':safe_field(sheet),'rule':'TIME_COLUMN','row_location':'header','fields':[safe_field(k) for k in times]})
   entries.append({'path':name,'table':safe_field(sheet),'rows':len(rows),'columns':[safe_field(k) for k in cols],'sha256':digest})
   for idx,r in enumerate(rows,2):
    # With a count column present, an unrecognized numeric data column is
    # treated conservatively as a measurement instead of silently allowed.
    metadata={'row_number','source_row','line','line_number','bytes','size','page','dpi','version','year','seed','index'}
    measures=[k for k in cols if k not in counts and k not in ids and k not in times and (MEASURE_RE.search(k) or (counts and k.lower() not in metadata)) and number(r.get(k)) is not None]
    small=[k for k in counts if number(r.get(k)) is not None and 1<=number(r.get(k))<=9]
    if small and measures:
     findings.append({'path':name,'table':safe_field(sheet),'rule':'SMALL_CELL_WITH_NUMERIC_VALUE','row_location':str(idx),'fields':[safe_field(k) for k in small],'numeric_measurement_fields':[safe_field(k) for k in measures]})
    masked=[k for k in counts if isinstance(r.get(k),str) and re.search(r'(<\s*10|[≤<]\s*9|1\s*[-–]\s*9|suppressed)',r[k],re.I)]
    if masked and measures:
     findings.append({'path':name,'table':safe_field(sheet),'rule':'MASKED_SMALL_COUNT_WITH_NUMERIC_VALUE','row_location':str(idx),'fields':[safe_field(k) for k in masked],'numeric_measurement_fields':[safe_field(k) for k in measures]})
 except Exception as ex:
  findings.append({'path':name,'rule':'READ_ERROR','row_location':'','fields':[],'error_type':type(ex).__name__})
 return findings,entries
def scan_paths(paths,include_restricted=False):
 findings=[];entries=[];excluded=[];seen=set()
 for raw in paths:
  p=Path(raw).resolve()
  fs=sorted(p.rglob('*')) if p.is_dir() else [p]
  for q in fs:
   if not q.is_file() or q.suffix.lower() not in DATA|SNIFF or q in seen:continue
   seen.add(q)
   if any(v in VENDOR for v in q.parts):continue
   if not include_restricted and any(v in RESTRICTED for v in q.parts):excluded.append(str(q));continue
   f,e=scan_bytes(q.read_bytes(),str(q));findings+=f;entries+=e
 return findings,entries,excluded
def staged(repo):
 cmd=['git','-C',str(repo)]
 p=subprocess.run(cmd+['diff','--cached','--name-only','--diff-filter=ACMR','-z'],capture_output=True,check=True)
 findings=[];entries=[]
 for v in p.stdout.decode('utf-8').split('\0'):
  if not v or Path(v).suffix.lower() not in DATA|SNIFF:continue
  b=subprocess.run(cmd+['show',':'+v],capture_output=True,check=True).stdout
  f,e=scan_bytes(b,'staged:'+v);findings+=f;entries+=e
 return findings,entries,[]
def pre_push(repo,lines):
 """Check blobs newly reachable from the actual refs Git is about to push."""
 cmd=['git','-C',str(repo)];findings=[];entries=[];seen=set()
 for raw in lines.splitlines():
  if not raw.strip():continue
  parts=raw.split()
  if len(parts)!=4:raise RuntimeError('invalid_pre_push_ref_record')
  local_ref,local_oid,remote_ref,remote_oid=parts
  if not re.fullmatch(r'[0-9a-f]{40,64}',local_oid) or not re.fullmatch(r'[0-9a-f]{40,64}',remote_oid):raise RuntimeError('invalid_pre_push_oid')
  if set(local_oid)=={'0'}:continue
  args=['rev-list','--objects',local_oid]
  if set(remote_oid)!={'0'}:
   exists=subprocess.run(cmd+['cat-file','-e',remote_oid],capture_output=True).returncode==0
   if exists:args.append('^'+remote_oid)
  objects=subprocess.run(cmd+args,capture_output=True,check=True).stdout.decode('utf-8').splitlines()
  for item in objects:
   oid,sep,path=item.partition(' ')
   if not sep or Path(path).suffix.lower() not in DATA|SNIFF or (oid,path) in seen:continue
   seen.add((oid,path));b=subprocess.run(cmd+['cat-file','blob',oid],capture_output=True,check=True).stdout
   f,e=scan_bytes(b,'push:'+local_oid+':'+path)
   findings+=f;entries+=e
 # Also check the current workspace so an uncommitted unsafe payload is not
 # accidentally treated as ready by a successful push hook invocation.
 f,e,x=scan_paths([repo]);return findings+f,entries+e,x
def main():
 a=argparse.ArgumentParser(description=__doc__);a.add_argument('paths',nargs='*');a.add_argument('--report',required=True);mode=a.add_mutually_exclusive_group();mode.add_argument('--staged',action='store_true');mode.add_argument('--pre-push',action='store_true');a.add_argument('--repo',default='.');a.add_argument('--include-restricted',action='store_true')
 args=a.parse_args()
 try:
  if args.staged:f,e,x=staged(Path(args.repo).resolve())
  elif args.pre_push:f,e,x=pre_push(Path(args.repo).resolve(),sys.stdin.read())
  else:
   if not args.paths:a.error('one or more paths are required')
   f,e,x=scan_paths(args.paths,args.include_restricted)
  report={'gate_version':VERSION,'status':'FAIL' if f else 'PASS','files_scanned':len({r['path'] for r in e}),'findings_count':len(f),'findings':f,'file_inventory':e,'restricted_files_not_read':x,'cell_values_in_report':False}
 except Exception as ex:report={'gate_version':VERSION,'status':'FAIL','findings_count':1,'findings':[{'rule':'EXECUTION_ERROR','error_type':type(ex).__name__}],'cell_values_in_report':False}
 out=Path(args.report).resolve();out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
 print(json.dumps({'status':report['status'],'files_scanned':report.get('files_scanned',0),'findings_count':report['findings_count'],'report':str(out)},ensure_ascii=False))
 return 0 if report['status']=='PASS' else 1
if __name__=='__main__':raise SystemExit(main())
