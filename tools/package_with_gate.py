#!/usr/bin/env python3
"""Package a curated payload only after its current bytes pass the data gate."""
from pathlib import Path
import argparse,hashlib,json,os,sys,zipfile
from small_cell_gate import scan_paths,scan_bytes,RESTRICTED
def main():
 p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--output',required=True);p.add_argument('--report',required=True);a=p.parse_args()
 src=Path(a.source).resolve();out=Path(a.output).resolve();report=Path(a.report).resolve()
 if not src.is_dir() or out.is_relative_to(src) or report.is_relative_to(src):raise RuntimeError('payload_output_and_report_must_be_separate')
 fs=sorted(f for f in src.rglob('*') if f.is_file())
 if any(any(v in RESTRICTED or v=='.git' for v in f.relative_to(src).parts) for f in fs):raise RuntimeError('restricted_or_git_payload_refused')
 hits,entries,_=scan_paths([src],include_restricted=True)
 receipt={'status':'FAIL' if hits else 'PASS','source':str(src),'findings':hits,'file_inventory':entries,'cell_values_in_report':False}
 report.parent.mkdir(parents=True,exist_ok=True)
 if hits:
  report.write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps({'status':'FAIL','findings_count':len(hits),'report':str(report)}));return 1
 # The files are hash-pinned and reread just once for ZIP creation. Mutation
 # between preflight and packaging fails without replacing an existing archive.
 pin={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in fs}
 out.parent.mkdir(parents=True,exist_ok=True);temp=out.with_name(out.name+'.gate_tmp')
 try:
  with zipfile.ZipFile(temp,'w',compression=zipfile.ZIP_DEFLATED) as z:
   for f in fs:
    b=f.read_bytes()
    if hashlib.sha256(b).hexdigest()!=pin[str(f)]:raise RuntimeError('payload_changed_after_gate')
    z.writestr(f.relative_to(src).as_posix(),b)
  zbytes=temp.read_bytes();hits,checked=scan_bytes(zbytes,str(out))
  if hits:raise RuntimeError('archive_postflight_failed')
  receipt.update({'status':'PASS','output':str(out),'sha256':hashlib.sha256(zbytes).hexdigest(),'member_files':len(fs),'source_sha256':pin,'postflight_passed':True})
  os.replace(temp,out);report.write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
 except Exception:
  if temp.exists():temp.unlink()
  raise
 print(json.dumps({'status':'PASS','member_files':len(fs),'output':str(out),'sha256':receipt['sha256'],'report':str(report)},ensure_ascii=False));return 0
if __name__=='__main__':raise SystemExit(main())
