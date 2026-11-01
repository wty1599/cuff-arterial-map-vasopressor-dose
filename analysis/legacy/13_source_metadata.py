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
import csv,json,datetime
ROOT=Path('@CUFF_WORKSPACE@');OUT=ROOT/'phase0'
with open(OUT/'map_bp_mimic.csv',encoding='utf-8-sig',newline='') as f:b=list(csv.DictReader(f))
ids=','.join(x['itemid'] for x in b if x.get('linksto')=='chartevents')
sql=f"""\\set ON_ERROR_STOP on
SET default_transaction_read_only=on;
SET max_parallel_workers_per_gather=2;
\\copy (SELECT itemid,count(*) AS n_rows FROM mimiciv_icu.chartevents WHERE itemid IN ({ids}) GROUP BY itemid) TO '@CUFF_WORKSPACE@/phase0/mimic_bp_dictionary_source_counts.csv' CSV HEADER
\\copy (SELECT * FROM mimiciv_hosp.d_labitems WHERE label ~* 'creatin|kreatin') TO '@CUFF_WORKSPACE@/phase0/mimic_creatinine_dictionary.csv' CSV HEADER
\\copy (SELECT EXISTS(SELECT 1 FROM mimiciv_hosp.labevents l JOIN mimiciv_hosp.d_labitems d USING(itemid) WHERE d.label ~* '^creatinine$' AND d.fluid='Blood' AND l.valuenum IS NOT NULL AND l.charttime IS NOT NULL) AS timed_serum_creatinine_rows_exist) TO '@CUFF_WORKSPACE@/phase0/mimic_creatinine_presence.csv' CSV HEADER
"""
(ROOT/'sql/13_mimic_source_metadata.sql').write_text(sql,encoding='utf-8')
sql="""\\set ON_ERROR_STOP on
SET default_transaction_read_only=on;
\\copy (SELECT labname,labmeasurenameinterface,labmeasurenamesystem,count(*) AS n_rows,count(labresultoffset) AS timed_rows FROM public.lab WHERE labname ~* 'creatin|kreatin' GROUP BY 1,2,3) TO '@CUFF_WORKSPACE@/phase0/eicu_creatinine_dictionary_presence.csv' CSV HEADER
"""
(ROOT/'sql/13_eicu_source_metadata.sql').write_text(sql,encoding='utf-8')
p=OUT/'deviation_log.csv'
with p.open(encoding='utf-8-sig',newline='') as f:r=csv.DictReader(f);fields=r.fieldnames;rows=list(r)
if not any(r['id']=='D006' for r in rows):
    rows.append(dict(id='D006',change='SICdb治疗观察起点采用cases.ICUOffset，排除前置手术期；首次ICU按OffsetAfterFirstAdmission+ICUOffset排序',reason='当前官方文档明确case可能含ICU之前的手术；0秒不是ICU起点',timing='源文档复核和血压解码期间，SICdb配对、剂量层和一致性结果产生之前；'+datetime.datetime.now().astimezone().isoformat(),authorization='源时间语义纠错，预设ICU人群和配对阈值不变',affected_outputs='SICdb治疗区间、首次ICU与源回查'))
    with p.open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fields);w.writeheader();w.writerows(rows)
print('Source metadata queries prepared.')
