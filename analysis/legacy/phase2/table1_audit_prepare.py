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
from pathlib import Path
import pandas as pd
ROOT=Path('@CUFF_WORKSPACE@');R=ROOT/'phase2/_restricted';CODE=ROOT/'sql/phase2'
for db in ['mimic','eicu','amsterdam','sicdb']:
 p=pd.read_csv(R/f'table1_patients_{db}.csv',dtype={'sid':str,'patient':str})
 un=p[p.max_nee_observed.isna()];kn=p[p.max_nee_observed.notna()];nu=min(4,len(un));a=pd.concat([un.sample(n=nu,random_state=20261027),kn.sample(n=12-nu,random_state=20261027)])
 a[['sid','patient','sex','hospital_death','severity','admission_type','max_nee_observed','max_nee_time','max_nee_status']].to_csv(R/f'table1_audit_expected_{db}.csv',index=False)
 if db in ['mimic','eicu']:
  m=pd.read_csv(ROOT/f'phase0/map_vaso_{db}.csv',dtype=str);m=m[m.selected.eq('1')];m[['itemid','drug']].to_csv(R/f'table1_audit_drugmap_{db}.csv',index=False)
  pre=f"""\set ON_ERROR_STOP on
CREATE TEMP TABLE a(sid text,patient text,sex text,hospital_death text,severity text,admission_type text,max_nee_observed text,max_nee_time text,max_nee_status text);
\copy a FROM '@CUFF_WORKSPACE@/phase2/_restricted/table1_audit_expected_{db}.csv' CSV HEADER
CREATE TEMP TABLE dm(itemid text,drug text);
\copy dm FROM '@CUFF_WORKSPACE@/phase2/_restricted/table1_audit_drugmap_{db}.csv' CSV HEADER
BEGIN TRANSACTION READ ONLY;
"""
  if db=='mimic':
   q1="SELECT i.stay_id AS sid,p.gender,h.hospital_expire_flag,h.admission_type,s.sofa FROM a JOIN mimiciv_icu.icustays i ON i.stay_id=a.sid::bigint JOIN mimiciv_hosp.patients p USING(subject_id) JOIN mimiciv_hosp.admissions h ON h.hadm_id=i.hadm_id LEFT JOIN mimiciv_derived.first_day_sofa s ON s.stay_id=i.stay_id"
   q2="SELECT x.stay_id AS sid,x.orderid,x.itemid,x.starttime,x.endtime,x.rate,x.rateuom,x.patientweight,x.ordercategorydescription,x.statusdescription,dm.drug FROM a JOIN mimiciv_icu.inputevents x ON x.stay_id=a.sid::bigint JOIN dm ON x.itemid::text=dm.itemid WHERE x.rate>0 AND x.ordercategorydescription='Continuous Med' AND x.statusdescription NOT IN ('Rewritten','Bolus')"
  else:
   q1="SELECT p.patientunitstayid AS sid,p.gender,p.hospitaldischargestatus,p.unitstaytype,p.admissionweight,r.apachescore,r.apacheversion FROM a JOIN public.patient p ON p.patientunitstayid=a.sid::bigint LEFT JOIN public.apachepatientresult r ON r.patientunitstayid=p.patientunitstayid AND r.apacheversion='IVa'"
   q2="SELECT x.*,dm.drug FROM a JOIN public.infusiondrug x ON x.patientunitstayid=a.sid::bigint JOIN dm ON x.drugname=dm.itemid"
  suffix=''
  for name,q in [('baseline',q1),('drugs',q2)]:suffix+=f"\\copy ({q}) TO '@CUFF_WORKSPACE@/phase2/_restricted/table1_audit_source_{db}_{name}.csv' CSV HEADER\n"
  (CODE/f'table1_audit_{db}.sql').write_text(pre+suffix+'ROLLBACK;\n',encoding='utf-8')
print('48 source-audit patients sampled; 8 or more have known-dose maxima per database.')
