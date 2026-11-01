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
import sys,json,datetime,hashlib
ROOT=Path('@CUFF_WORKSPACE@')

import pandas as pd
OUT=ROOT/'phase3'; RES=OUT/'_restricted'; CODE=ROOT/'sql/phase3'
for p in [OUT,RES,OUT/'_logs',CODE]:p.mkdir(parents=True,exist_ok=True)
for db in ['mimic','eicu']:
    d=pd.read_parquet(ROOT/f'phase1/_restricted/pairs_{db}.parquet',columns=['sid','patient'])
    z=d.drop_duplicates('patient').sort_values('sid');assert z.sid.is_unique
    z[['sid']].to_csv(RES/f'population_ids_{db}.csv',index=False)

mimic=r'''\set ON_ERROR_STOP on
CREATE TEMP TABLE population_ids(sid bigint);
\copy population_ids FROM '@CUFF_WORKSPACE@/phase3/_restricted/population_ids_mimic.csv' CSV HEADER
BEGIN TRANSACTION READ ONLY;
\copy (SELECT i.stay_id AS sid,i.first_careunit,count(s.stay_id) sepsis_rows,count(*) FILTER(WHERE s.sepsis3 IS TRUE) sepsis_positive_rows,count(*) FILTER(WHERE s.sepsis3 IS FALSE) sepsis_negative_rows,count(*) FILTER(WHERE s.stay_id IS NOT NULL AND s.sepsis3 IS NULL) sepsis_null_rows FROM population_ids k JOIN mimiciv_icu.icustays i ON i.stay_id=k.sid LEFT JOIN mimiciv_derived.sepsis3 s ON s.stay_id=i.stay_id GROUP BY i.stay_id,i.first_careunit ORDER BY i.stay_id) TO '@CUFF_WORKSPACE@/phase3/_restricted/population_source_mimic.csv' CSV HEADER
\copy (SELECT first_careunit AS careunit,count(*) AS stays FROM mimiciv_icu.icustays GROUP BY first_careunit ORDER BY first_careunit) TO '@CUFF_WORKSPACE@/phase3/population_dictionary_mimic_units.csv' CSV HEADER
\copy (SELECT sepsis3,count(*) AS records,count(DISTINCT stay_id) AS stays FROM mimiciv_derived.sepsis3 GROUP BY sepsis3 ORDER BY sepsis3) TO '@CUFF_WORKSPACE@/phase3/population_dictionary_mimic_sepsis.csv' CSV HEADER
ROLLBACK;
'''
eicu=r'''\set ON_ERROR_STOP on
CREATE TEMP TABLE population_ids(sid bigint);
\copy population_ids FROM '@CUFF_WORKSPACE@/phase3/_restricted/population_ids_eicu.csv' CSV HEADER
BEGIN TRANSACTION READ ONLY;
\copy (SELECT p.patientunitstayid AS sid,p.unittype,p.apacheadmissiondx FROM population_ids k JOIN public.patient p ON p.patientunitstayid=k.sid ORDER BY p.patientunitstayid) TO '@CUFF_WORKSPACE@/phase3/_restricted/population_source_eicu.csv' CSV HEADER
\copy (SELECT unittype AS careunit,count(*) AS stays FROM public.patient GROUP BY unittype ORDER BY unittype) TO '@CUFF_WORKSPACE@/phase3/population_dictionary_eicu_units.csv' CSV HEADER
\copy (SELECT apacheadmissiondx AS diagnosis,count(*) AS stays FROM public.patient WHERE apacheadmissiondx ILIKE '%sepsis%' GROUP BY apacheadmissiondx ORDER BY apacheadmissiondx) TO '@CUFF_WORKSPACE@/phase3/population_dictionary_eicu_sepsis.csv' CSV HEADER
ROLLBACK;
'''
(CODE/'population_source_mimic.sql').write_text(mimic,encoding='utf-8')
(CODE/'population_source_eicu.sql').write_text(eicu,encoding='utf-8')
print('Prepared ID files and source SQL; no sensitivity endpoints computed.')
