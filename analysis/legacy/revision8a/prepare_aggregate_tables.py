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
import csv,json,hashlib,shutil,sys
sys.stdout.reconfigure(encoding='utf-8')
O=Path('@CUFF_WORKSPACE@')/'phase3/revision8a'
BASE=O.parent
D=O/'data'
D.mkdir(parents=True,exist_ok=True)
(O/'qa').mkdir(parents=True,exist_ok=True)
def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(p,rows):
 with p.open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
for p in (BASE/'revision4/data').glob('*.csv'):shutil.copy2(p,D/p.name)
for name in ['A1_cuff_ge70.csv','A2_likelihood_ratios.csv']:
 shutil.copy2(BASE/'icm_revision4/analysis'/name,D/name)
dbs=['MIMIC-IV','eICU','AmsterdamUMCdb','SICdb'];strata=['all','<=0.1','>0.1-0.25','>0.25-0.5','>0.5-1.0','>1.0']
miss=read(D/'S2_missed_and_false_alarm.csv');deep=read(D/'S2_deep_missed_hypotension.csv')
pooled=read(D/'S2_pooled_missed_and_false_alarm.csv');pdeep=read(D/'S2_pooled_deep_missed_hypotension.csv')
def one(rows,**filt):
 r=[x for x in rows if all(x.get(k)==str(v) for k,v in filt.items())]
 assert len(r)==1,(filt,len(r))
 return r[0]
metrics=[('missed','Missed hypotension rate'),('deep','Deep missed hypotension rate (arterial reference MAP <55 mmHg)'),('false_alarm','False alarm rate')]
print('ENDPOINTS',sorted(set(x['endpoint'] for x in miss)))
# Source terminology is retained in the selection rather than inferred by value.
fa=[x for x in sorted(set(x['endpoint'] for x in miss)) if 'false' in x.lower()][0]
metrics[2]=('false_alarm',fa)
table2=[]
for db in dbs+['Random-effects pooled']:
 for st in strata:
  for short,metric in metrics:
   source='S2_deep_missed_hypotension.csv' if short=='deep' else 'S2_missed_and_false_alarm.csv'
   src_rows=deep if short=='deep' else miss
   if db=='Random-effects pooled':
    source='S2_pooled_deep_missed_hypotension.csv' if short=='deep' else 'S2_pooled_missed_and_false_alarm.csv'
    src_rows=pdeep if short=='deep' else pooled
    r=one(src_rows,**({'metric':metric} if short=='deep' else {'endpoint':metric}),stratum='known' if st=='all' else st)
    dp=r.get('denominator_patients','')
    components=[]
    if not dp:
     included=r.get('included','').split(';')
     names=dict(zip(['mimic','eicu','amsterdam','sicdb'],dbs))
     comp=[one(deep,database=names[k],scenario='Primary five-minute pairing',metric=metric,stratum='known' if st=='all' else st) for k in included]
     dp=str(sum(int(x['denominator_patients']) for x in comp));components=[x['row_id'] for x in comp]
   else:
    r=one(src_rows,database=db,scenario='Primary five-minute pairing',stratum=st,**({'metric':metric} if short=='deep' else {'endpoint':metric}))
    dp=r.get('denominator_patients','0');components=[]
   denom=int(float(r['denominator'] or 0));num=int(float(r['numerator'] or 0))
   assert 0<=num<=denom,(db,st,short)
   est=r.get('estimate','');lo=r.get('ci_lower','');hi=r.get('ci_upper','')
   status='not estimable' if denom==0 else 'adequate'
   if denom and (denom<30 or (lo and hi and float(hi)-float(lo)>0.30)):
    status='insufficient precision'
   if db=='AmsterdamUMCdb' and st=='>1.0': status='insufficient precision' if denom else 'not estimable'
   boundary=num in (0,denom) if denom else False
   table2.append(dict(row_id=f'T2-{len(table2)+1:03d}',database=db,stratum='known' if db=='Random-effects pooled' and st=='all' else st,metric=short,numerator=num,denominator=denom,denominator_patients=dp,estimate=est,ci_lower=lo,ci_upper=hi,I2=r.get('I2',''),precision=status,boundary=boundary,source_csv=source,source_row_id=r['row_id'],patients_component_row_ids=';'.join(components),included_databases=r.get('included',''),note='Patient count sums disjoint database-specific denominator patients' if components else ''))
write(D/'Table2_rev8a.csv',table2)
# Complete false-alarm display for Supplementary Table 2, including unknown dose.
s2=[dict(x,precision_display='not estimable' if int(x['denominator'])==0 else ('insufficient precision' if (int(x['denominator'])<30 or (x['ci_lower'] and x['ci_upper'] and float(x['ci_upper'])-float(x['ci_lower'])>0.30) or (x['database']=='AmsterdamUMCdb' and x['stratum']=='>1.0')) else 'adequate')) for x in miss if x['endpoint']==fa and x['scenario']=='Primary five-minute pairing']
for x in s2:
 if int(x['denominator'])>0 and int(x['numerator']) in (0,int(x['denominator'])):
  x['precision_display']+='; boundary'
write(D/'S2_false_alarm_display.csv',s2)
cross=read(D/'S4_dose_by_pressure.csv');depth=[]
for db in dbs[:2]:
 for st in strata[1:]:
  rows=[one(cross,database=db,scenario='Primary five-minute pairing',stratum=st,category=cat) for cat in ['<55','55-59','60-64']]
  total=sum(int(x['denominator']) for x in rows)
  assert total==int(rows[0]['layer_low_arterial_pairs'])
  for r in rows:
   depth.append(dict(row_id=f'S5-depth-{len(depth)+1:03d}',database=db,stratum=st,arterial_MAP_band=r['category'],pairs=int(r['denominator']),low_arterial_pairs=total,fraction=int(r['denominator'])/total,source_csv='S4_dose_by_pressure.csv',source_row_id=r['row_id'],numerator_column='denominator',denominator_column='layer_low_arterial_pairs'))
write(D/'S5_hypotension_depth_distribution.csv',depth)
for file,out,fil in [('A1_cuff_ge70.csv','S17_cuff_ge70.csv',None),('A2_likelihood_ratios.csv','S18_likelihood_ratios.csv','five_strata')]:
 rows=[]
 for line,r in enumerate(read(D/file),2):
  if fil and r['grouping']!=fil:continue
  rows.append(dict(row_id=f'{out[:3]}-{len(rows)+1:03d}',**r,source_csv=file,source_line=line))
 write(D/out,rows)
receipt={'table2_rows':len(table2),'s2_false_alarm_rows':len(s2),'depth_rows':len(depth),
         'new_S17_rows':len(read(D/'S17_cuff_ge70.csv')),'new_S18_rows':len(read(D/'S18_likelihood_ratios.csv')),
         'bootstrap_recalculated':False}
(O/'qa/aggregate_build_receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps(receipt))
