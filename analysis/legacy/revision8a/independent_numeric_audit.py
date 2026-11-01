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
"""Recount existing frozen paired files; never refit or re-bootstrap."""
from pathlib import Path
import sys, csv, math, json, hashlib, datetime

import duckdb

ROOT=Path(r'@CUFF_WORKSPACE@')
OUT=ROOT/'phase3/revision8a'
SRC=ROOT/'phase3/revision4/data'
A=ROOT/'phase3/icm_revision4/analysis'
DB={'MIMIC-IV':'mimic','eICU':'eicu','AmsterdamUMCdb':'amsterdam','SICdb':'sicdb'}
STRATA={'all':'true','known':'nee IS NOT NULL','unknown':'nee IS NULL','<=0.1':'nee<=0.1','>0.1-0.25':'nee>0.1 AND nee<=0.25','>0.25-0.5':'nee>0.25 AND nee<=0.5','>0.5-1.0':'nee>0.5 AND nee<=1.0','>1.0':'nee>1.0','<=0.25':'nee<=0.25','>0.25-1.0':'nee>0.25 AND nee<=1.0'}
CAT={'60-64':'art_map>=60 AND art_map<65','55-59':'art_map>=55 AND art_map<60','<55':'art_map<55'}
checks=[]
optional_unreported=[]
def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(path,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def numeric(x):
    try:return float(x)
    except (TypeError,ValueError):return float('nan')
def equal(a,b):
    x=numeric(a);y=numeric(b)
    if math.isnan(x) or math.isnan(y):return math.isnan(x) and math.isnan(y)
    if math.isinf(x) or math.isinf(y):return x==y
    return abs(x-y)<=1e-10*(1+abs(y))
def check(db,st,metric,source,row,field,reported,recounted,note=''):
    checks.append(dict(check_id=f'N{len(checks)+1:04d}',database=db,stratum=st,metric=metric,source_csv=source,source_row=row,field=field,reported=reported,recounted=recounted,passed=equal(reported,recounted),note=note))
def ratio(a,b):return a/b if b else float('nan')
con=duckdb.connect()
con.execute("SET memory_limit='1GB'")
con.execute('SET threads=2')
counts={}
for db,short in DB.items():
    file=ROOT/f'phase1/_restricted/pairs_{short}.parquet'
    con.execute(f"CREATE VIEW {short} AS SELECT * FROM read_parquet('{file.as_posix()}')")
    for st,where in STRATA.items():
        q=f"""SELECT count(*) all_pairs,count(DISTINCT patient) all_patients,
        count(*) FILTER(WHERE art_map<65 AND cuff_map<65) TP,
        count(*) FILTER(WHERE art_map<65 AND cuff_map>=65) FN,
        count(*) FILTER(WHERE art_map>=65 AND cuff_map<65) FP,
        count(*) FILTER(WHERE art_map>=65 AND cuff_map>=65) TN,
        count(DISTINCT patient) FILTER(WHERE art_map<65) low_patients,
        count(DISTINCT patient) FILTER(WHERE art_map>=65) normal_patients,
        count(*) FILTER(WHERE art_map<55 AND cuff_map>=65) deep_num,
        count(*) FILTER(WHERE art_map<55) deep_den,
        count(DISTINCT patient) FILTER(WHERE art_map<55) deep_patients,
        count(*) FILTER(WHERE art_map<60 AND cuff_map>=65) lt60_num,
        count(*) FILTER(WHERE art_map<60) lt60_den,
        count(DISTINCT patient) FILTER(WHERE art_map<60) lt60_patients,
        count(*) FILTER(WHERE art_map<65 AND cuff_map>=70) A1_num,
        count(*) FILTER(WHERE cuff_map>=70) A1_den,
        count(DISTINCT patient) FILTER(WHERE cuff_map>=70) A1_patients,
        count(*) FILTER(WHERE cuff_map>=65 AND cuff_map<70) cuff_65_69,
        count(*) FILTER(WHERE art_map<65 AND cuff_map>=65 AND cuff_map<70) cuff_65_69_low
        FROM {short} WHERE {where}"""
        r=con.execute(q)
        d=dict(zip([col[0] for col in r.description],r.fetchone()))
        for cat,condition in CAT.items():
            r=con.execute(f"SELECT count(*),count(*) FILTER(WHERE cuff_map>=65),count(DISTINCT patient) FROM {short} WHERE ({where}) AND ({condition})").fetchone()
            d[cat+'_den'],d[cat+'_num'],d[cat+'_patients']=r
        counts[(db,st)]=d
        check(db,st,'four-cell conservation','frozen paired parquet','independent recount','TP+FN+FP+TN',d['TP']+d['FN']+d['FP']+d['TN'],d['all_pairs'],'No new model or bootstrap')
        check(db,st,'pressure-band conservation','frozen paired parquet','independent recount','three arterial bands',sum(d[cat+'_den'] for cat in CAT),d['TP']+d['FN'])
    d=counts[(db,'all')]
    check(db,'all','dose-stratum conservation','frozen paired parquet','independent recount','all = five known strata + unknown',d['all_pairs'],sum(counts[(db,st)]['all_pairs'] for st in list(STRATA)[3:8])+counts[(db,'unknown')]['all_pairs'])

three_rates=[]
for filename in ['S2_missed_and_false_alarm.csv','S2_deep_missed_hypotension.csv','S4_dose_by_pressure.csv']:
    rows=read(SRC/filename)
    for row in rows:
        db=row['database'];st=row['stratum'];d=counts[(db,st)]
        metric=row.get('endpoint',row.get('metric',''))
        if metric=='Missed hypotension rate':
            num,den,patients=d['FN'],d['TP']+d['FN'],d['low_patients']
        elif metric=='False-alarm rate':
            num,den,patients=d['FP'],d['FP']+d['TN'],d['normal_patients']
        elif metric.startswith('Deep missed'):
            num,den,patients=d['deep_num'],d['deep_den'],d['deep_patients']
        elif metric.startswith('Missed hypotension conditional'):
            num,den,patients=d['lt60_num'],d['lt60_den'],d['lt60_patients']
        elif metric.startswith('Missed hypotension rate within'):
            cat=row['category'];num,den,patients=d[cat+'_num'],d[cat+'_den'],d[cat+'_patients']
        else:raise ValueError('Unrecognized aggregate metric')
        expected={'numerator':num,'denominator':den,'denominator_patients':patients,'estimate':ratio(num,den)}
        for key in ['all_pairs','all_patients']:
            if key in row:expected[key]=d[key]
        if 'layer_low_arterial_pairs' in row:
            optional={'layer_low_arterial_pairs':d['TP']+d['FN'],'layer_pairs':d['all_pairs'],'layer_low_arterial_prevalence':ratio(d['TP']+d['FN'],d['all_pairs'])}
            for key,value in optional.items():
                if row[key]!='':expected[key]=value
                else:optional_unreported.append({'source_csv':filename,'source_row':row['row_id'],'field':key,'status':'Optional metadata not reported in the source aggregate; not a numerical discrepancy'})
        for key,value in expected.items():
            check(db,st,metric,filename,row['row_id'],key,row[key],value)
        lo,hi=numeric(row['ci_lower']),numeric(row['ci_upper'])
        valid=math.isnan(lo) and math.isnan(hi) or 0<=lo<=hi<=1
        check(db,st,metric,filename,row['row_id'],'CI range/order only',int(valid),1,'Existing CI retained; no bootstrap recomputation')
        if st=='all' and (metric in ('Missed hypotension rate','False-alarm rate') or metric.startswith('Deep missed')):
            three_rates.append(dict(database=db,metric=metric,numerator=num,denominator=den,denominator_patients=patients,all_pairs=d['all_pairs'],estimate=ratio(num,den),ci_lower=row['ci_lower'],ci_upper=row['ci_upper'],source_csv=filename,source_row=row['row_id'],CI_status='Reused unchanged; bootstrap not rerun'))

weight_den=sum(counts[(db,'all')]['TP']+counts[(db,'all')]['FN'] for db in ['MIMIC-IV','eICU'])
weights={cat:sum(counts[(db,'all')][cat+'_den'] for db in ['MIMIC-IV','eICU'])/weight_den for cat in CAT}
for row in read(ROOT/'phase2/phase2_C_standard_weights.csv'):
    check('main databases','all','direct-standardization weights','phase2_C_standard_weights.csv',row['arterial_bin'],'weight',row['weight'],weights[row['arterial_bin']],'All low arterial MAP pairs, including unknown dose')
skipped=[]
for row in read(SRC/'S5_standardized_missed_hypotension.csv'):
    if row['scenario']!='Primary five-minute pairing':
        skipped.append(row['row_id']);continue
    d=counts[(row['database'],row['stratum'])]
    estimate=sum(weights[cat]*ratio(d[cat+'_num'],d[cat+'_den']) for cat in CAT)
    for key,value in {'estimate':estimate,'raw_missed_rate':ratio(d['FN'],d['TP']+d['FN']),'low_arterial_pairs':d['TP']+d['FN'],'denominator_patients':d['low_patients'],'minimum_cell_pairs':min(d[cat+'_den'] for cat in CAT)}.items():
        check(row['database'],row['stratum'],'direct-standardized missed rate','S5_standardized_missed_hypotension.csv',row['row_id'],key,row[key],value,'Point estimate only; original CI unchanged')

for filename in ['A1_cuff_ge70.csv','A2_likelihood_ratios.csv','A_contingency_counts.csv']:
    for lineno,row in enumerate(read(A/filename),2):
        db={'mimic':'MIMIC-IV','eicu':'eICU'}[row['database']];st=row['stratum'];d=counts[(db,st)]
        if filename.startswith('A1'):
            expected={'all_pairs':d['all_pairs'],'all_patients':d['all_patients'],'numerator':d['A1_num'],'denominator':d['A1_den'],'denominator_patients':d['A1_patients'],'estimate':ratio(d['A1_num'],d['A1_den']),'cuff_ge65_risk':ratio(d['FN'],d['FN']+d['TN'])}
            metric='P(arterial MAP <65 | cuff MAP >=70)'
        elif filename.startswith('A2'):
            specificity=ratio(d['TN'],d['FP']+d['TN']);miss=ratio(d['FN'],d['TP']+d['FN'])
            lr=ratio(miss,specificity) if row['metric']=='LR_negative' else ratio(1-miss,1-specificity)
            expected={key:d[key] for key in ['TP','FN','FP','TN','all_pairs','all_patients']}
            expected.update(estimate=lr,missed_rate=miss,specificity=specificity)
            metric=row['metric']
        else:
            expected={key:d[key] for key in ['TP','FN','FP','TN','A1_num','A1_den','all_pairs','all_patients','cuff_65_69','cuff_65_69_low']}
            expected.update(all_normal_cuff=d['FN']+d['TN'],A1_den_patients=d['A1_patients'])
            metric='four-cell and cuff-threshold counts'
        for key,value in expected.items():check(db,st,metric,filename,lineno,key,row[key],value)

for db,field,value in [('AmsterdamUMCdb','FN',1323),('AmsterdamUMCdb','low_den',2598),('AmsterdamUMCdb','FP',618),('AmsterdamUMCdb','normal_den',4351),('SICdb','FN',481),('SICdb','low_den',1342),('SICdb','FP',449),('SICdb','normal_den',1126)]:
    d=counts[(db,'all')];actual=d['TP']+d['FN'] if field=='low_den' else d['FP']+d['TN'] if field=='normal_den' else d[field]
    check(db,'all','user supplied overall counts','user instruction','rev8a',field,value,actual)

write(OUT/'reviews/independent_numeric_checks.csv',checks)
write(OUT/'data/independent_three_rates.csv',three_rates)
quad=[dict(database=db,**{key:counts[(db,'all')][key] for key in ['all_pairs','all_patients','TP','FN','FP','TN']}) for db in DB]
write(OUT/'data/independent_four_cell_counts.csv',quad)
sources=[SRC/f for f in ['S2_missed_and_false_alarm.csv','S2_deep_missed_hypotension.csv','S4_dose_by_pressure.csv','S5_standardized_missed_hypotension.csv']]+[A/f for f in ['A1_cuff_ge70.csv','A2_likelihood_ratios.csv','A_contingency_counts.csv']]
receipt=dict(completed_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),checks=len(checks),passed=sum(x['passed'] for x in checks),failed=sum(not x['passed'] for x in checks),recounted_existing_pair_files=True,source_verification='Separate SICdb raw-file audit only',models_rerun=False,bootstrap_rerun=False,existing_CIs='Retained unchanged; numerical bootstrap validity not independently rerun',optional_blank_metadata_not_treated_as_discrepancy=len(optional_unreported),skipped_recovered_dose_standardized_rows=skipped,standardization_weights=weights,sources_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
(OUT/'qa/independent_numeric_receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in receipt.items() if k!='sources_sha256'},ensure_ascii=False))
print(json.dumps(quad,ensure_ascii=False))
print(json.dumps([row for row in checks if not row['passed']],ensure_ascii=False))
