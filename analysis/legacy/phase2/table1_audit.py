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
"""Read source tables/files independently; never imports the Table 1 conversion."""
from pathlib import Path
import sys,re,json
ROOT=Path('@CUFF_WORKSPACE@');R=ROOT/'phase2/_restricted';OUT=ROOT/'phase2';S0=ROOT/'phase0/_restricted';SIC=Path('@SICDB_ROOT@')

import duckdb,pandas as pd,numpy as np
def n(x):
 try:return float(x)
 except (TypeError,ValueError):return np.nan
def midpoint(s):
 z=re.findall('[0-9]+(?:\\.[0-9]+)?',str(s));return (float(z[0])+float(z[1]))/2 if len(z)==2 else np.nan
def conversion(rate,unit,wt,drug):
 r=n(rate);u=str(unit).lower().replace(' ','').replace('µg','mcg').replace('μg','mcg').replace('nanograms','ng').replace('micrograms','mcg').replace('units','u').replace('unit','u').replace('hour','hr').replace('uur','hr').replace('minutes','min').replace('minute','min').replace('ie','u')
 if not np.isfinite(r) or r<=0:return np.nan
 if not (u.endswith('/min') or u.endswith('/hr')):return np.nan
 perkg='/kg/' in u;tm=60 if u.endswith('/hr') else 1
 if drug=='vasopressin':
  if not u.startswith('u/') or (perkg and not 30<=wt<=300):return np.nan
  return r/tm*(wt if perkg else 1)*2.5
 mass={'mcg':1,'mg':1000,'g':1000000,'ng':.001}.get(u.split('/')[0]);fac={'norepinephrine':1,'epinephrine':1,'phenylephrine':.1,'dopamine':.01,'angiotensin_ii':10}.get(drug)
 if mass is None or fac is None or (not perkg and not 30<=wt<=300):return np.nan
 return r*mass*fac/tm/(1 if perkg else wt)
def eq(a,b):return (pd.isna(a) and pd.isna(b)) or (not pd.isna(a) and not pd.isna(b) and (str(a)==str(b) or (np.isfinite(n(a)) and np.isfinite(n(b)) and np.isclose(n(a),n(b),rtol=1e-7,atol=1e-8))))
def sqlcsv(p):return f"read_csv_auto('{p.as_posix()}',all_varchar=true)"
summ=[];details=[]
for db in ['mimic','eicu','amsterdam','sicdb']:
 e=pd.read_csv(R/f'table1_audit_expected_{db}.csv',dtype={'sid':str,'patient':str});stored=pd.read_csv(R/f'table1_patients_{db}.csv',dtype={'sid':str}).set_index('sid');c=duckdb.connect();c.execute("SET memory_limit='700MB';SET threads=1")
 c.execute(f"CREATE TEMP VIEW aa AS SELECT * FROM {sqlcsv(R/f'table1_audit_expected_{db}.csv')}")
 if db in ['mimic','eicu']:
  b=pd.read_csv(R/f'table1_audit_source_{db}_baseline.csv',dtype=str,keep_default_na=False);raw=pd.read_csv(R/f'table1_audit_source_{db}_drugs.csv',dtype=str,keep_default_na=False)
 else:
  path=S0/'amsterdam_stays.csv' if db=='amsterdam' else SIC/'cases.csv.gz';key='admissionid' if db=='amsterdam' else 'CaseID'
  b=c.execute(f'SELECT x.* FROM {sqlcsv(path)} x JOIN aa ON x.{key}=aa.sid').df();b['sid']=b[key]
  rpath=S0/'amsterdam_drugs.csv.gz' if db=='amsterdam' else SIC/'medication.csv.gz';rk='admissionid' if db=='amsterdam' else 'CaseID'
  raw=c.execute(f'SELECT x.* FROM {sqlcsv(rpath)} x JOIN aa ON x.{rk}=aa.sid').df()
  b.to_csv(R/f'table1_audit_source_{db}_baseline.csv',index=False);raw.to_csv(R/f'table1_audit_source_{db}_drugs.csv',index=False)
  mp=pd.read_csv(ROOT/f'phase0/map_vaso_{db}.csv',dtype=str);mp=mp[mp.selected.eq('1')].set_index('itemid');ref=pd.read_csv(SIC/'d_references.csv.gz',dtype=str,keep_default_na=False).set_index('ReferenceGlobalID').ReferenceValue.to_dict()
 for a in e.itertuples(index=False):
  bb=b[b.sid.astype(str).eq(a.sid)].iloc[0];vals=[]
  if db=='mimic':
   sex={'M':'M','F':'F'}.get(bb.gender);death=n(bb.hospital_expire_flag);score=n(bb.sofa);adm=bb.admission_type
   rr=raw[raw.sid.astype(str).eq(a.sid)].copy();rr['start']=pd.to_datetime(rr.starttime).astype('datetime64[s]').astype('int64');rr['stop']=pd.to_datetime(rr.endtime).astype('datetime64[s]').astype('int64');rr=rr[(rr.start<=a.max_nee_time)&(rr.stop>a.max_nee_time)]
   for x in rr.itertuples(index=False):vals.append(conversion(x.rate,x.rateuom,n(x.patientweight),x.drug))
  elif db=='eicu':
   sex={'male':'M','female':'F'}.get(bb.gender.lower());death={'alive':0.,'expired':1.}.get(bb.hospitaldischargestatus.lower(),np.nan);sc=pd.to_numeric(b[b.sid.eq(a.sid)].apachescore,errors='coerce');sc=sc[sc>=0];score=float(sc.iloc[0]) if sc.nunique()==1 else np.nan;adm=bb.unitstaytype
   rr=raw[raw.patientunitstayid.astype(str).eq(a.sid)].copy();rr['tm']=pd.to_numeric(rr.infusionoffset)*60;gr=[]
   for (drug,tm),g in rr.groupby(['drug','tm']):
    doses=[];ongoing=False
    for x in g.itertuples(index=False):
     us=[z for z in re.findall(r'\(([^()]*)\)',x.drugname) if re.search(r'/min|/hr|/hour',z,re.I)];u=us[-1] if us else '';r=n(x.drugrate);ongoing=ongoing or r>0 or ('ml/hr' in u.lower() and n(x.infusionrate)>0);v=conversion(r,u,n(bb.admissionweight),drug)
     if np.isfinite(v):doses.append(v)
    v=float(np.median(doses)) if doses and max(doses)-min(doses)<=max(1e-6,max(doses)*.01) else np.nan;gr.append(dict(drug=drug,tm=tm,dose=v,ongoing=ongoing))
   if gr:
    z=pd.DataFrame(gr).sort_values(['drug','tm']);z['stop']=np.minimum(z.groupby('drug').tm.shift(-1).fillna(np.inf),z.tm+3600);z=z[z.ongoing&(z.tm<=a.max_nee_time)&(z.stop>a.max_nee_time)];vals=z.dose.tolist()
  elif db=='amsterdam':
   sex={'man':'M','vrouw':'F'}.get(str(bb.gender).lower());death=np.nan;score=np.nan;adm={'0':'planned','1':'unplanned'}.get(str(bb.urgency));wt=midpoint(bb.weightgroup)
   rr=raw[raw.admissionid.astype(str).eq(a.sid)].copy();rr=rr[rr.itemid.isin(mp.index)&pd.to_numeric(rr.iscontinuous,errors='coerce').eq(1)&rr.ordercategory.str.contains('Spuitpompen',case=False)]
   for x in rr.itertuples(index=False):
    if n(x.start)/1000<=a.max_nee_time<n(x.stop)/1000 and n(x.dose)>0:
     u=x.doseunit+('/kg/' if n(x.doserateperkg)==1 else '/')+x.doserateunit;vals.append(conversion(x.dose,u,wt,mp.loc[x.itemid,'drug']))
  else:
   sex={'male':'M','female':'F'}.get(ref.get(str(bb.Sex),''));disp=ref.get(str(bb.HospitalDischargeType),'');death={'Deceased':1.,'Sterbefall':1.,'Survived':0.,'Entlassung aus dem Krankenhaus':0.,'Transfer in ein anderes Krankenhaus':0.,'Entlassung gegen Revers':0.}.get(disp,np.nan);score=np.nan;adm=ref.get(str(bb.SurgicalAdmissionType));wt=n(bb.WeightOnAdmission)/1000
   rr=raw[raw.CaseID.astype(str).eq(a.sid)];rr=rr[rr.DrugID.isin(mp.index)&pd.to_numeric(rr.IsSingleDose,errors='coerce').eq(0)]
   for x in rr.itertuples(index=False):
    if n(x.Offset)<=a.max_nee_time<n(x.OffsetDrugEnd) and n(x.AmountPerMinute)>0:vals.append(conversion(x.AmountPerMinute,mp.loc[x.DrugID,'unitname']+'/min',wt,mp.loc[x.DrugID,'drug']))
  total=float(sum(vals)) if vals and all(np.isfinite(v) for v in vals) else np.nan
  checks={'sex':eq(a.sex,sex),'hospital_death':eq(a.hospital_death,death),'severity':eq(a.severity,score),'admission_type':eq(a.admission_type,adm)}
  peak_checked=np.isfinite(a.max_nee_observed);peak_match=eq(a.max_nee_observed,total) if peak_checked else None
  age_proxy_match=eq(stored.loc[a.sid,'age'],midpoint(bb.agegroup)) if db=='amsterdam' else None
  details.append(dict(database=db,sid=a.sid,**{k+'_match':v for k,v in checks.items()},age_proxy_match=age_proxy_match,peak_checked=peak_checked,peak_matches=peak_match,peak_expected=a.max_nee_observed,peak_source=total))
 sub=pd.DataFrame(details);sub=sub[sub.database.eq(db)];summ.append(dict(database=db,patients=len(sub),baseline_field_checks=len(sub)*4+(len(sub) if db=='amsterdam' else 0),baseline_field_discrepancies=int((~sub[['sex_match','hospital_death_match','severity_match','admission_type_match']]).sum().sum())+int(sub.age_proxy_match.eq(False).sum()),peak_points_checked=int(sub.peak_checked.sum()),peak_point_discrepancies=int(sub.loc[sub.peak_checked,'peak_matches'].eq(False).sum()),source_method='independent_PostgreSQL' if db in ['mimic','eicu'] else 'independent_DuckDB_original_CSV',scope='baseline fields and source dose at reported peak; not independent whole-trajectory maximum reconstruction'))
 c.close()
pd.DataFrame(details).to_csv(R/'table1_audit_comparison.csv',index=False);pd.DataFrame(summ).to_csv(OUT/'table1_source_audit.csv',index=False)
print(json.dumps(summ,ensure_ascii=False))
