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
import sys,json,re,importlib.util
ROOT=Path('@CUFF_WORKSPACE@');SRC=ROOT/'phase0/_restricted';OUT=ROOT/'phase2';RES=OUT/'_restricted';RES.mkdir(parents=True,exist_ok=True)

import duckdb,pandas as pd,numpy as np
spec=importlib.util.spec_from_file_location('old_core',ROOT/'sql/08_core.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
DBS=['mimic','eicu','amsterdam','sicdb'];DRUGS=['norepinephrine','epinephrine','phenylephrine','dopamine','vasopressin','angiotensin_ii']
SIC=Path('@SICDB_ROOT@')
ROWS=[];QC=[];FLOW=[]
def num(x):return pd.to_numeric(x,errors='coerce')
def reg(c,name,d):
 d=d.copy()
 for col in d:
  if isinstance(d[col].dtype,pd.StringDtype):d[col]=d[col].astype(object).where(d[col].notna(),None)
 c.register(name,d)
def csvq(p):return f"read_csv_auto('{p.as_posix()}',all_varchar=true)"
def clean_category(s):
 s=s.astype('string').str.strip();return s.mask(s.str.lower().isin(['','unknown','nan','none','other/unknown','null']))
def interval_mid(s):
 z=re.findall(r'\d+(?:\.\d+)?',str(s));return (float(z[0])+float(z[1]))/2 if len(z)==2 else np.nan
def add(db,pt,var,label,vals,source,definition,status='observed',binary=False,categorical=False):
 v=pd.Series(vals,index=pt.index);N=len(pt);n=int(v.notna().sum());base=dict(database=db,variable=var,label=label,patients=N,n_observed=n,n_missing=N-n,missing_fraction=(N-n)/N,source=source,definition=definition,status=status)
 if categorical or binary:
  cats=[1.,0.] if binary else sorted(v.dropna().astype(str).unique())
  for cat in cats:
   count=int(v.eq(cat).sum()) if binary else int(v.astype('string').eq(cat).fillna(False).sum())
   ROWS.append(dict(**base,category=('yes' if cat==1 else 'no') if binary else cat,numerator=count,denominator=n,percent=100*count/n if n else np.nan,median=np.nan,q1=np.nan,q3=np.nan))
  if not cats: ROWS.append(dict(**base,category='not_available',numerator=np.nan,denominator=0,percent=np.nan,median=np.nan,q1=np.nan,q3=np.nan))
 else:
  y=num(v);z=y.dropna();ROWS.append(dict(**base,category='median_IQR',numerator=np.nan,denominator=n,percent=np.nan,median=z.median() if len(z) else np.nan,q1=z.quantile(.25) if len(z) else np.nan,q3=z.quantile(.75) if len(z) else np.nan))

def invalid_weight_intervals(db,c,pt,st):
 bad=[];w=st.set_index('sid').weight
 if db=='mimic':
  reg(c,'dose_stays',st[['sid','dose_fallback_weight']])
  r=c.execute(f'''SELECT x.* FROM {csvq(SRC/'mimic_drugs.csv')} x JOIN dose_stays s ON x.stay_id=s.sid WHERE (CASE WHEN try_cast(x.patientweight AS DOUBLE)>0 THEN try_cast(x.patientweight AS DOUBLE) ELSE s.dose_fallback_weight END <30 OR CASE WHEN try_cast(x.patientweight AS DOUBLE)>0 THEN try_cast(x.patientweight AS DOUBLE) ELSE s.dose_fallback_weight END >300) AND try_cast(x.rate AS DOUBLE)>0 AND x.ordercategorydescription='Continuous Med' AND x.statusdescription NOT IN ('Rewritten','Bolus')''').df()
  mp=pd.read_csv(ROOT/'phase0/map_vaso_mimic.csv',dtype=str).set_index('itemid').drug.to_dict();fw=st.set_index('sid').dose_fallback_weight.to_dict()
  for r in r.itertuples(index=False):
   dr=mp.get(str(r.itemid));ww=old.number(r.patientweight);ww=ww if ww>0 else fw.get(str(r.stay_id),np.nan);u=str(r.rateuom).lower().replace(' ','');kg='/kg' in u;need=(not kg and dr!='vasopressin') or (kg and dr=='vasopressin')
   if dr and need and np.isfinite(old.convert_rate(r.rate,r.rateuom,ww,dr)[0]):bad.append(dict(sid=str(r.stay_id),start=pd.Timestamp(r.starttime).timestamp(),stop=pd.Timestamp(r.endtime).timestamp()))
 elif db=='eicu':
  ids=st[(st.weight<30)|(st.weight>300)][['sid','weight']];reg(c,'bad_stays',ids)
  r=c.execute(f'''SELECT x.*,s.weight FROM {csvq(SRC/'eicu_drugs.csv')} x JOIN bad_stays s ON x.patientunitstayid=s.sid''').df()
  mp=pd.read_csv(ROOT/'phase0/map_vaso_eicu.csv',dtype=str);mp=mp[mp.selected.eq('1')].set_index('itemid').drug.to_dict();allrows=[]
  for r in r.itertuples(index=False):
   dr=mp.get(r.drugname)
   if not dr:continue
   us=[v for v in re.findall(r'\(([^()]*)\)',r.drugname) if re.search(r'/min|/hr|/hour',v,re.I)];u=us[-1] if us else '';kg='/kg' in u;need=(not kg and dr!='vasopressin') or (kg and dr=='vasopressin')
   entered=need and old.number(r.drugrate)>0 and np.isfinite(old.convert_rate(r.drugrate,u,r.weight,dr)[0])
   allrows.append(dict(sid=str(r.patientunitstayid),drug=dr,start=old.number(r.infusionoffset)*60,bad=entered))
  if allrows:
   z=pd.DataFrame(allrows).groupby(['sid','drug','start'],as_index=False).bad.any().sort_values(['sid','drug','start']);z['stop']=np.minimum(z.groupby(['sid','drug']).start.shift(-1).fillna(np.inf),z.start+3600);bad=z.loc[z.bad,['sid','start','stop']].to_dict('records')
 elif db in ['amsterdam','sicdb']:
  ids=st[(st.weight<30)|(st.weight>300)].sid.tolist()
  if ids:
   # Only source-converted mass rates can require this invalid admission weight.
   z=c.execute("SELECT sid,drug,start,\"end\" stop,dose,rawid FROM drugs").df();z=z[z.sid.astype(str).isin(ids)&z.dose.notna()]
   if db=='sicdb':bad=z[z.drug.ne('vasopressin')][['sid','start','stop']].to_dict('records')
   else:
    r=pd.read_csv(SRC/'amsterdam_drugs.csv.gz',dtype=str,keep_default_na=False);r=r[r.admissionid.isin(ids)];mp=pd.read_csv(ROOT/'phase0/map_vaso_amsterdam.csv',dtype=str).set_index('itemid').drug.to_dict()
    for a in r.itertuples(index=False):
     dr=mp.get(a.itemid);kg=old.number(a.doserateperkg)==1;need=(not kg and dr!='vasopressin')or(kg and dr=='vasopressin');u=a.doseunit+('/kg/' if kg else '/')+a.doserateunit
     if dr and need and old.number(a.iscontinuous)==1 and np.isfinite(old.convert_rate(a.dose,u,w.get(a.admissionid,np.nan),dr)[0]):bad.append(dict(sid=a.admissionid,start=old.number(a.start)/1000,stop=old.number(a.stop)/1000))
 return pd.DataFrame(bad,columns=['sid','start','stop'])

def maxdose(db,c,pt,st):
 reg(c,'chosen',pt[['sid']]);dr=c.execute('SELECT d.* FROM drugs d JOIN chosen k USING(sid)').df();ex=c.execute('SELECT d.* FROM drug_episodes60 d JOIN chosen k USING(sid)').df();ep=c.execute('SELECT d.* FROM episodes60 d JOIN chosen k USING(sid)').df();bw=invalid_weight_intervals(db,c,pt,st)
 dr.to_csv(RES/f'table1_drugs_{db}.csv',index=False);ex.to_csv(RES/f'table1_expected_{db}.csv',index=False);ep.to_csv(RES/f'table1_episodes_{db}.csv',index=False);bw.to_csv(RES/f'table1_invalid_weight_intervals_{db}.csv',index=False)
 cols=['dose','unknown','expected','therapy','badweight']+DRUGS;events=[]
 def interval_events(z,a,b,extra):
  f=pd.DataFrame({'sid':z.sid.astype(str),'t':z[a].to_numpy()});g=pd.DataFrame({'sid':z.sid.astype(str),'t':z[b].to_numpy()})
  for col in cols:
   x=extra.get(col,0);f[col]=x;g[col]=-np.asarray(x)
  events.extend([f,g])
 interval_events(dr,'start','end',{'dose':dr.dose.fillna(0).to_numpy(),'unknown':dr.dose.isna().astype(int).to_numpy(),**{x:dr.drug.eq(x).astype(int).to_numpy() for x in DRUGS}})
 interval_events(ex,'start','stop',{'expected':1});interval_events(ep,'start','stop',{'therapy':1})
 if len(bw):interval_events(bw,'start','stop',{'badweight':1})
 ee=pd.concat(events,ignore_index=True).groupby(['sid','t'],sort=True)[cols].sum().reset_index();out=[]
 for sid,g in ee.groupby('sid',sort=False):
  tt=g.t.to_numpy(float);dur=np.diff(tt);v=np.cumsum(g[cols].to_numpy(float),axis=0)[:-1];within=(v[:,3]>0)&(dur>0);active=(v[:,5:]>0).sum(axis=1);total=np.round(v[:,0],12)
  known=within&(v[:,2]>0)&(v[:,2]==active)&(v[:,1]==0)&(v[:,4]==0)&(total>=0)&(total<=5)
  therapy=float(dur[within].sum());kt=float(dur[known].sum());maxval=float(total[known].max()) if known.any() else np.nan
  idx=np.where(known&(total==maxval))[0] if known.any() else []
  maxt=float(tt[idx[0]]) if len(idx) else np.nan
  out.append(dict(sid=str(sid),max_nee_observed=maxval,max_nee_time=maxt,treatment_seconds=therapy,known_seconds=kt,unknown_seconds=therapy-kt,known_duration_fraction=kt/therapy if therapy else np.nan,max_nee_status='complete_within_defined_episodes' if kt>0 and np.isclose(kt,therapy,rtol=0,atol=1e-6) else 'partial_observation' if kt>0 else 'no_known_segment',over5_seconds=float(dur[within&(total>5)].sum())))
 d=pd.DataFrame(out);expected=ep.assign(seconds=ep.stop-ep.start).groupby('sid').seconds.sum();d['expected_seconds']=d.sid.map(expected)
 assert np.allclose(d.treatment_seconds,d.expected_seconds,rtol=0,atol=1e-5)
 assert d.sid.is_unique and len(d)==len(pt)
 d.to_csv(RES/f'table1_maxdose_{db}.csv',index=False)
 use=dr.groupby('drug').sid.apply(lambda x:set(x.astype(str))).to_dict()
 return d,use

def main():
 refs=pd.read_csv(SIC/'d_references.csv.gz',dtype=str,keep_default_na=False);ref=refs.set_index('ReferenceGlobalID').ReferenceValue.to_dict()
 sic=pd.read_csv(SIC/'cases.csv.gz',dtype=str,keep_default_na=False).set_index('CaseID')
 for db in DBS:
  c=duckdb.connect(str(SRC/f'{db}.duckdb'),read_only=True);c.execute("SET memory_limit='700MB';SET threads=1")
  d=c.execute(f"SELECT * FROM read_parquet('{(ROOT/f'phase1/_restricted/pairs_{db}.parquet').as_posix()}')").df();d.sid=d.sid.astype(str)
  assert d.cid.is_unique and d.groupby('patient').sid.nunique().max()==1
  pt=d.sort_values(['patient','t','cid']).drop_duplicates('patient').copy().reset_index(drop=True);pc=d.groupby('patient').size();pt['pair_count']=pt.patient.map(pc)
  reg(c,'chosen',pt[['sid']]);st=c.execute('SELECT s.* FROM stays s JOIN chosen k USING(sid)').df();st.sid=st.sid.astype(str)
  defs={'age':'ICU admission age; de-identification/top coding retained','weight':'phase1 range-checked admission/earliest weight, kg','bmi':'phase1 range-checked BMI from height and weight'};statuses={'age':'observed','weight':'observed','bmi':'observed'}
  if db=='mimic':
   x=pd.read_csv(RES/'table1_mimic_source.csv',dtype=str,keep_default_na=False).set_index('sid');assert x.index.is_unique
   pt['hospital_death']=num(pt.sid.map(x.hospital_expire_flag));pt['severity']=num(pt.sid.map(x.sofa));pt.loc[~pt.severity.between(0,24),'severity']=np.nan
   pt['unit_type']=pt.sid.map(x.first_careunit);pt['admission_type']=pt.sid.map(x.admission_type);pt['admission_source']=pt.sid.map(x.admission_location)
   severity_label='SOFA：入 ICU 首日';severity_source='mimiciv_derived.first_day_sofa.sofa';severity_status='existing_derived_first_day'
   death_source='mimiciv_hosp.admissions.hospital_expire_flag';admitdef='Hospital admission_type, not harmonized ICU medical/surgical admission'
  elif db=='eicu':
   x=pd.read_csv(RES/'table1_eicu_source.csv',dtype=str,keep_default_na=False);base=x.drop_duplicates('sid').set_index('sid');ap=x[x.apacheversion.eq('IVa')].copy();ap['score']=num(ap.apachescore).where(num(ap.apachescore)>=0);conf=ap.groupby('sid').score.nunique().gt(1);sco=ap.groupby('sid').score.first();sco.loc[conf[conf].index]=np.nan
   pt['hospital_death']=pt.sid.map(base.hospitaldischargestatus).str.lower().map({'expired':1.,'alive':0.});pt['severity']=pt.sid.map(sco);pt['unit_type']=pt.sid.map(base.unittype);pt['admission_type']=pt.sid.map(base.unitstaytype);pt['admission_source']=pt.sid.map(base.unitadmitsource)
   severity_label='APACHE IVa';severity_source='public.apachepatientresult[apacheversion=IVa].apachescore';severity_status='existing_derived_IVa';death_source='public.patient.hospitaldischargestatus';admitdef='Source unitstaytype (admit/readmit/etc.), not medical/surgical classification';statuses['age']='top_coded_age_retained'
   QC.append(dict(database=db,check='conflicting_IVa_scores',value=int(conf.sum()),unit='patients'))
  elif db=='amsterdam':
   x=pd.read_csv(SRC/'amsterdam_stays.csv',dtype=str,keep_default_na=False).set_index('admissionid');pt['age_group']=pt.sid.map(x.agegroup);pt['age']=pt.age_group.map(interval_mid)
   pt['hospital_death']=np.nan;pt['severity']=np.nan;pt['unit_type']=pt.sid.map(x.location);pt['admission_type']=pt.sid.map(x.urgency).map({'0':'planned','1':'unplanned'});pt['admission_source']=pt.sid.map(x.origin)
   for v in ['age','weight','bmi']:statuses[v]='closed_interval_midpoint_proxy';defs[v]='Closed released intervals only; midpoint proxy; open intervals unknown. BMI derives from weight/height midpoints.'
   severity_label='严重程度评分（本轮指定量表）';severity_source='not_available_in_selected_existing_sources';severity_status='not_extracted_no_new_score';death_source='not_available_complete_hospital_disposition';admitdef='Admission urgency in admissions: 0 planned,1 unplanned';pt['bmi40']=np.nan
  else:
   x=sic;pt['sex']=pt.sid.map(x.Sex).map(ref).map({'male':'M','female':'F'});disp=pt.sid.map(x.HospitalDischargeType).map(ref)
   pt['hospital_death']=disp.map({'Deceased':1.,'Sterbefall':1.,'Survived':0.,'Entlassung aus dem Krankenhaus':0.,'Transfer in ein anderes Krankenhaus':0.,'Entlassung gegen Revers':0.})
   pt['hospital_disposition_source']=disp;pt['severity']=np.nan;pt['unit_type']=pt.sid.map(x.HospitalUnit).map(ref);pt['admission_type']=pt.sid.map(x.SurgicalAdmissionType).map(ref);pt['admission_source']=pt.sid.map(x.ReferringUnit).map(ref)
   severity_label='严重程度评分（本轮指定量表）';severity_source='no_SOFA_or_APACHE_IVa_used;existing_cases_saps3_not_substituted';severity_status='not_requested_incomparable_SAPS3_not_substituted';death_source='cases.HospitalDischargeType+d_references';admitdef='SurgicalAdmissionType source categories';statuses.update(age='rounded_5_years_top_coded90',weight='released_grams_to_kg_rounded',bmi='derived_from_rounded_weight_height')
  for v,label in [('age','年龄（岁）'),('weight','体重（kg）'),('bmi','BMI（kg/m²）')]:add(db,pt,v,label,pt[v],'phase1 baseline; source interval exception for Amsterdam',defs[v],statuses[v])
  add(db,pt,'sex','性别',pt.sex,'phase1 baseline; SICdb cases.Sex dictionary','Recorded sex; unknown explicitly missing',categorical=True)
  add(db,pt,'bmi_gt40','BMI >40',pt.bmi40,'phase1 range-checked BMI','Amsterdam not classified from interval proxy' if db=='amsterdam' else 'BMI>40 among available BMI','not_applicable_interval_proxy' if db=='amsterdam' else 'observed_or_rounded_proxy',binary=True)
  add(db,pt,'atrial_fibrillation','房颤记录',pt.af_model,'phase1 af_model','MIMIC rhythm or diagnosis; eICU diagnosis-covered; Amsterdam rhythm; SIC not ascertained','not_ascertained' if db=='sicdb' else 'source_specific_ever_recorded',binary=True)
  for v,label in [('unit_type','ICU 单元类别（源口径）'),('admission_type','入院类别（源口径）'),('admission_source','入院来源（源口径）')]:add(db,pt,v,label,clean_category(pt[v]),'index stay source fields',admitdef if v=='admission_type' else ('SIC hospital unit is last case unit; others source ICU unit' if v=='unit_type' else 'Source admission origin; not harmonized'),categorical=True)
  if db=='amsterdam':add(db,pt,'age_group','原始年龄区间',pt.age_group,'admissions.agegroup','Original released interval labels, including open intervals',categorical=True)
  add(db,pt,'severity_score',severity_label,pt.severity,severity_source,'Existing database-specific score; not comparable across score systems',severity_status)
  add(db,pt,'hospital_death','院内死亡',pt.hospital_death,death_source,'Death during index hospital stay; no ICU-death or unspecified death-date substitution','not_observable' if db=='amsterdam' else 'source_hospital_disposition',binary=True)
  add(db,pt,'pairs_per_patient','每患者配对数',pt.pair_count,'phase1 pairs parquet','All phase1 primary +/-5 min pairs for the patient')
  ma,use=maxdose(db,c,pt,st);pt=pt.merge(ma,on='sid',validate='one_to_one')
  for dr in DRUGS:
   pt['use_'+dr]=pt.sid.isin(use.get(dr,set())).astype(float)
   add(db,pt,'vasopressor_use',dr,pt['use_'+dr],'phase0 positive treatment segments in index ICU stay','Any source-recorded positive treatment segment; categories overlap; no record is not proof of no treatment','documented_treatment',binary=True)
   ROWS[-2]['category']=dr;ROWS.pop() # one yes row per drug, denominator all paired patients
  add(db,pt,'max_nee_observed','可判定治疗片段中观测到的最大 NE-eq',pt.max_nee_observed,'all phase0 treatment segments with phase1 range rules','Maximum of summed active components over all positive-length known segments; not limited to BP pairs; unknown segments not imputed','partial_observation_reported')
  add(db,pt,'known_duration_fraction','治疗期总 NE-eq 可判定时长比例',pt.known_duration_fraction,'all episode and drug interval endpoints','Per patient known dose duration / all defined treatment episode duration')
  add(db,pt,'dose_observation_status','治疗期剂量观察完整性',pt.max_nee_status,'all treatment intervals','Complete only within defined recorded/inferred treatment episodes',categorical=True)
  add(db,pt,'treatment_hours','已定义治疗期时长（小时）',pt.treatment_seconds/3600,'phase0 episodes60','Merged defined episodes; eICU inferred; includes short bridged gaps')
  pt.to_csv(RES/f'table1_patients_{db}.csv',index=False)
  # New-field audit sample: baseline fields plus the maximum-dose source time.
  pt.sample(n=min(12,len(pt)),random_state=20261027)[['sid','patient','sex','hospital_death','severity','admission_type','max_nee_observed','max_nee_time','max_nee_status']].to_csv(RES/f'table1_audit_expected_{db}.csv',index=False)
  QC.extend([dict(database=db,check='patients_unique',value=int(pt.patient.is_unique),unit='boolean'),dict(database=db,check='patients',value=len(pt),unit='patients'),dict(database=db,check='pairs_sum',value=int(pt.pair_count.sum()),unit='pairs'),dict(database=db,check='treatment_duration_conserved',value=1,unit='boolean'),dict(database=db,check='maxdose_no_known_patients',value=int(pt.max_nee_observed.isna().sum()),unit='patients')])
  c.close();print(json.dumps({'database':db,'patients':len(pt),'table1':'computed','maxdose_known_patients':int(pt.max_nee_observed.notna().sum())}),flush=True)
  # Deliver progressively; no row-level data enter these aggregate files.
  pd.DataFrame(ROWS).to_csv(OUT/'phase2_table1.csv',index=False,encoding='utf-8-sig');pd.DataFrame(QC).to_csv(OUT/'table1_quality_checks.csv',index=False)
 # Existing flow products; no overlapping exclusions are added.
 unitmap={'all_source_stays':'stays','adult_first_icu':'patients','adult_first_order_unknown':'stays','vasopressor_patients':'patients','arterial_during_therapy':'patients','both_bp_during_therapy':'patients','paired5_patients':'patients','paired5_at_least10_patients':'patients'}
 for db in DBS:
  ff=pd.read_csv(ROOT/f'phase0/phase0_feasibility_{db}.csv')
  for r in ff[ff.section.eq('flow')].itertuples(index=False):FLOW.append(dict(database=db,stage=r.metric,value=int(r.value),unit=unitmap.get(r.metric,'count'),source=f'phase0_feasibility_{db}.csv',definition='Existing phase0 cohort; branches not additive',status='source_scope_retained'))
  for met,u in [('eligible_cuff_readings','readings'),('paired5_readings','pairs'),('paired2_readings','pairs')]:
   val=ff[(ff.section=='primary')&(ff.metric==met)].iloc[0].value;FLOW.append(dict(database=db,stage=met,value=int(val),unit=u,source=f'phase0_feasibility_{db}.csv',definition='Primary cohort anchors/pairs; units distinct',status='source_scope_retained'))
  f1=pd.read_csv(ROOT/'phase1/phase1_flow.csv');zz=f1[f1.database.eq(db)]
  for r in zz.itertuples(index=False):FLOW.append(dict(database=db,stage='phase1_dose_'+str(r.stratum),value=int(r.pairs),unit='pairs',source='phase1_flow.csv',definition='Mutually exclusive clean dose strata; patient counts across dose strata overlap',status='phase1_clean'))
 pd.DataFrame(FLOW).to_csv(OUT/'phase2_flow_numbers.csv',index=False,encoding='utf-8-sig')
 print('TABLE1_AND_FLOW_COMPLETE',flush=True)
if __name__=='__main__':main()
