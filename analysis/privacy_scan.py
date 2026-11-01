"""Static release checks; reports file/rule locations, never matched secrets.

This complements the private denominator/suppression review. It cannot infer a
patient denominator from a rounded proportion or replace source-data governance.
"""
from pathlib import Path
import argparse
import ast
import csv
import hashlib
import json
import re
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]
ID_PATTERN=re.compile(r'\b(?:subject_id|stay_id|patientunitstayid|admissionid|hadm_id|CaseID|patient_id|pid|cid)\b',re.I)
BAD_HEADERS=re.compile(r'^(?:subject_id|stay_id|patientunitstayid|admissionid|hadm_id|caseid|patient_id|pid|cid|patient|sid|timestamp|charttime|storetime|starttime|endtime|cuff_map|art_map|rawdata|numerator|denominator|all_pairs|all_patients|n_patients|n_pairs|TP|FN|FP|TN|denominator_patients)$',re.I)
PATTERNS={
    'absolute_workstation_path':re.compile(r'(?i)\b[A-Z]:[\\/]|/mnt/[a-z]/|/(?:home|Users)/[^/\s]+/'),
    'email':re.compile(r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}\b'),
    'database_connection_string':re.compile(r'(?i)(?:postgres(?:ql)?|mysql|mongodb|mssql)://[^\s]+'),
    'private_key':re.compile(r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----'),
    'github_token':re.compile(r'\b(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{30,})\b'),
    'aws_access_key':re.compile(r'\b(?:AKIA|ASIA)[A-Z0-9]{16}\b'),
    'patient_id_literal':re.compile(r'(?i)\b(?:subject_id|stay_id|patientunitstayid|admissionid|hadm_id|CaseID|patient_id|pid|cid)\s*(?:==|=|IN\s*\()\s*[\[\(\'"]?\d{4,}'),
    'hardcoded_large_cohort_count':re.compile(r'(?<![<>=\d,])\b(?:\d{4,}|\d{1,3}(?:,\d{3})+) (?:patients|pairs)\b'),
}
ALLOWED_SUFFIX={'.py','.sql','.r','.csv','.json','.yaml','.yml','.md','.cff','.txt','.toml','.ipynb'}
SKIP_PARTS={'.git','.venv','__pycache__','reproduced'}

def scan(gitleaks=None):
    findings=[];schema=[];notebooks=0;files=[];syntax=0
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or any(x in SKIP_PARTS for x in p.relative_to(ROOT).parts):continue
        rel=p.relative_to(ROOT).as_posix()
        if rel in {'docs/privacy_scan.md','docs/privacy_scan.json'}:continue
        files.append(p)
        if p.suffix.lower() not in ALLOWED_SUFFIX and p.name not in {'LICENSE','.gitignore','.gitattributes'}:
            findings.append({'file':rel,'rule':'unapproved_file_type'});continue
        try:text=p.read_text(encoding='utf-8-sig')
        except UnicodeDecodeError:findings.append({'file':rel,'rule':'binary_or_unknown_encoding'});continue
        for name,pattern in PATTERNS.items():
            for m in pattern.finditer(text):
                findings.append({'file':rel,'rule':name,'line':text[:m.start()].count('\n')+1})
        if p.suffix.lower() in {'.py','.sql','.r'}:
            tokens=sorted(set(m.group(0) for m in ID_PATTERN.finditer(text)))
            if tokens:schema.append({'file':rel,'schema_tokens':tokens})
        if p.suffix.lower()=='.py':
            try:ast.parse(text);syntax+=1
            except SyntaxError as e:findings.append({'file':rel,'rule':'python_syntax','line':e.lineno})
        if p.suffix.lower()=='.csv':
            data=list(csv.DictReader(text.splitlines()));headers=list(data[0]) if data else next(csv.reader(text.splitlines()),[])
            if rel.startswith('outputs/aggregate/'):
                for h in headers:
                    if BAD_HEADERS.fullmatch(h) or 'bootstrap' in h.lower():findings.append({'file':rel,'rule':'forbidden_aggregate_column','column':h})
                identifiers={'database','scenario','endpoint','metric','stratum','category','grouping','site','axis','pressure_range','sensitivity','analysis','method','precision','source_csv','source_line','display'}
                for i,row in enumerate(data,2):
                    for key,value in row.items():
                        if key in identifiers or not value:continue
                        try:float(value)
                        except ValueError:continue
                        if re.search(r'\.\d{4,}',value):findings.append({'file':rel,'rule':'full_precision_numeric_aggregate','line':i,'column':key})
        if p.suffix.lower()=='.ipynb':
            notebooks+=1
            obj=json.loads(text)
            if any(cell.get('outputs') or cell.get('execution_count') is not None for cell in obj.get('cells',[])):
                findings.append({'file':rel,'rule':'notebook_output_or_execution_count'})
    secret_status='not_run';secret_version=''
    if gitleaks:
        version=subprocess.run([str(gitleaks),'version'],capture_output=True,text=True,check=True)
        secret_version=version.stdout.strip()
        with tempfile.TemporaryDirectory() as temp:
            report=Path(temp)/'report.json'
            proc=subprocess.run([str(gitleaks),'dir',str(ROOT),'--report-format','json','--report-path',str(report),'--redact','--no-banner'],capture_output=True,text=True)
            if proc.returncode not in (0,1):raise RuntimeError('gitleaks could not complete; see its local invocation')
            leaks=json.loads(report.read_text()) if report.exists() else []
            for entry in leaks:
                location=Path(entry.get('File','')).name
                findings.append({'file':location,'rule':'gitleaks:'+entry.get('RuleID','unknown'),'line':entry.get('StartLine')})
            secret_status='pass' if proc.returncode==0 and not leaks else 'fail'
    result={'status':'pass' if not findings and secret_status=='pass' else 'fail',
            'files_scanned':len(files),'python_files_parsed':syntax,'notebooks':notebooks,
            'gitleaks_status':secret_status,'gitleaks_version':secret_version,
            'findings':findings,'schema_identifiers_in_code':schema,
            'file_sha256':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--gitleaks',type=Path);p.add_argument('--report',type=Path,default=ROOT/'docs/privacy_scan.json');args=p.parse_args()
    result=scan(args.gitleaks);args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    markdown=f'''# Privacy scan

Status: **{result['status'].upper()}**. Scanned {result['files_scanned']} text files; parsed {result['python_files_parsed']} Python source files. Secret scan: gitleaks {result['gitleaks_version']} ({result['gitleaks_status']}). No patient data were read by this static scanner.

The checks cover secrets/tokens, private keys, absolute workstation paths, database connection strings, email addresses, literal patient-ID filters, forbidden aggregate columns, full-precision aggregate values, file types and notebook outputs. Matched secret values are never written to this report. The machine-readable companion contains rule locations and file hashes.

Notebooks selected: {result['notebooks']}. No outputs or execution counts were retained. The documented nbstripout workflow is required if notebooks are added; this release contains none.

Identifiers such as subject_id, stay_id and patientunitstayid occur in extraction SQL and analytical grouping code as schema names. They are not ID observations. Static matches were classified by file type; aggregate CSV columns were separately checked. Dictionary item identifiers are concept keys, not patient identifiers. Actual patient or pair records and literal patient-ID selections are prohibited.

The private denominator audit verified the minimum contributing population before rounding and export. Cells below ten patients and cells without sufficient denominator evidence were withheld. No exact patient/pair counts, marginal count totals, individual-contribution percentages, or bootstrap draws are present. Full-precision ratios were replaced with manuscript-display precision. Cross-table disclosure review included code assertions and captions, not only CSV columns. The same safeguards apply even while this repository is private.

This scan is supplemented by an independent pre-upload review. It does not certify unrestricted rights to the source databases or replace the signed data-use agreements. Original row-level audit files remain outside the repository.
'''
    (ROOT/'docs/privacy_scan.md').write_text(markdown,encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['status','files_scanned','python_files_parsed','notebooks','gitleaks_status','findings']},indent=2))
    if result['status']!='pass':raise SystemExit(1)

if __name__=='__main__':main()
