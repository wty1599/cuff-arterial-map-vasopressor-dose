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
import csv
from datetime import datetime
from zoneinfo import ZoneInfo
from importlib import import_module
maps=import_module('02_prepare_maps'); ROOT=maps.ROOT;OUT=maps.OUT
def literal(s):return "'"+str(s).replace("'","''")+"'"
def values(rows,key):return ','.join(literal(x[key]) for x in rows)
def export(db,queries):
    text='\\set ON_ERROR_STOP on\nSET default_transaction_read_only=on;\nSET max_parallel_workers_per_gather=2;\n'
    for name,q in queries:
        text+=f"\\copy ({q}) TO '@CUFF_WORKSPACE@/phase0/_restricted/{db}_{name}.csv' CSV HEADER\n"
    (ROOT/'sql'/f'03_export_{db}.sql').write_text(text,encoding='utf-8')
b=maps.read(OUT/'map_bp_mimic.csv');v=maps.read(OUT/'map_vaso_mimic.csv');c=maps.read(OUT/'map_covariates_mimic.csv')
bpids=values([x for x in b if x['selected']=='1']+c,'itemid');drugids=values([x for x in v if x['selected']=='1'],'itemid')
export('mimic',[
('stays',"SELECT i.*,p.gender,p.anchor_age+(extract(year from i.intime)::int-p.anchor_year) AS age,p.dod,a.hospital_expire_flag,a.deathtime,a.dischtime,row_number() OVER(PARTITION BY i.subject_id ORDER BY i.intime,i.stay_id) AS patient_icu_rank FROM mimiciv_icu.icustays i JOIN mimiciv_hosp.patients p USING(subject_id) JOIN mimiciv_hosp.admissions a USING(hadm_id)"),
('drugs',f"SELECT stay_id,itemid,starttime,endtime,rate,rateuom,amount,amountuom,patientweight,orderid,linkorderid,ordercategoryname,ordercategorydescription,statusdescription FROM mimiciv_icu.inputevents WHERE itemid IN ({drugids})"),
('bp_cov',f"SELECT stay_id,charttime,itemid,value,valuenum,valueuom FROM mimiciv_icu.chartevents WHERE itemid IN ({bpids})"),
('omr',"SELECT o.* FROM mimiciv_hosp.omr o WHERE result_name ~* 'weight|height'"),
('af_dx',"SELECT DISTINCT x.hadm_id FROM mimiciv_hosp.diagnoses_icd x JOIN mimiciv_hosp.d_icd_diagnoses d USING(icd_code,icd_version) WHERE d.long_title ~* 'atrial fibrillation'"),
('arterial_line',"SELECT p.stay_id,p.itemid,p.starttime,p.endtime FROM mimiciv_icu.procedureevents p JOIN mimiciv_icu.d_items d USING(itemid) WHERE d.label ~* 'arterial.*line|arterial.*catheter|art.*line'"),
])
e=maps.read(OUT/'map_vaso_eicu.csv');ed=values([x for x in e if x['selected']=='1'],'itemid')
export('eicu',[
('stays','SELECT * FROM public.patient'),
('drugs',f'SELECT * FROM public.infusiondrug WHERE drugname IN ({ed})'),
('art','SELECT vitalperiodicid,patientunitstayid,observationoffset,systemicsystolic,systemicdiastolic,systemicmean FROM public.vitalperiodic WHERE systemicmean IS NOT NULL OR systemicsystolic IS NOT NULL OR systemicdiastolic IS NOT NULL'),
('cuff','SELECT vitalaperiodicid,patientunitstayid,observationoffset,noninvasivesystolic,noninvasivediastolic,noninvasivemean FROM public.vitalaperiodic WHERE noninvasivemean IS NOT NULL OR noninvasivesystolic IS NOT NULL OR noninvasivediastolic IS NOT NULL'),
('rhythm',"SELECT patientunitstayid,nursingchartoffset,nursingchartcelltypevallabel,nursingchartcelltypevalname,nursingchartvalue FROM public.nursecharting WHERE concat_ws(' ',nursingchartcelltypevallabel,nursingchartcelltypevalname) ~* 'rhythm'"),
('af_dx',"SELECT DISTINCT patientunitstayid,diagnosisstring FROM public.diagnosis WHERE diagnosisstring ~* 'atrial fibrillation|atrial fibr'"),
('arterial_line',"SELECT patientunitstayid,treatmentoffset,treatmentstring FROM public.treatment WHERE treatmentstring ~* 'arterial.*catheter|arterial.*line'"),
('hospital_infusion_coverage',"SELECT p.hospitalid,count(DISTINCT p.patientunitstayid) AS all_stays,count(DISTINCT i.patientunitstayid) AS infusion_any_stays FROM public.patient p LEFT JOIN (SELECT DISTINCT patientunitstayid FROM public.infusiondrug) i USING(patientunitstayid) GROUP BY p.hospitalid"),
])
print('Readonly source-export SQL generated from dictionary mappings.')
