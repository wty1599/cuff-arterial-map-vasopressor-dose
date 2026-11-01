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

"""Correct a source-time restriction on admission height; preserve all pairs/doses."""
from importlib import import_module
g=import_module('08_core');globals().update({k:v for k,v in vars(g).items() if not k.startswith('_')})
import shutil
c=connect('mimic');cq=csvq(R/'mimic_bp_cov.csv');f=OUT/'phase0_feasibility_mimic.csv';backup=OUT/'_logs/phase0_feasibility_mimic_before_height_correction.csv'
if not backup.exists():shutil.copy2(f,backup)
before=pd.read_csv(f);pairguard=c.execute('SELECT count(*),sum(cuff_map),sum(art_map),sum(nee),sum(t) FROM pairs').fetchone()
c.execute(f"""CREATE OR REPLACE TABLE height_admission_corrected AS WITH x AS (
 SELECT x.stay_id sid,epoch(try_cast(x.charttime AS TIMESTAMP)) t,s.t0,m.label,try_cast(x.valuenum AS DOUBLE)*CASE WHEN m.label='Height' THEN 2.54 ELSE 1 END height
 FROM {cq} x JOIN cov_map m ON x.itemid=m.itemid JOIN stays s ON x.stay_id=s.sid
 WHERE m.label ILIKE '%height%' AND try_cast(x.valuenum AS DOUBLE)>0)
 SELECT sid,arg_min(height,struct_pack(distance := abs(t-t0),source_t := t,label_name := label)) height FROM x GROUP BY sid""")
c.execute('UPDATE stays SET height=h.height,bmi=stays.weight/(h.height/100)^2 FROM height_admission_corrected h WHERE stays.sid=h.sid')
c.execute('UPDATE pairs SET bmi=s.bmi FROM stays s WHERE pairs.sid=s.sid')
c.execute('CREATE OR REPLACE TABLE treated_primary AS SELECT * FROM stays s WHERE s.is_primary AND EXISTS(SELECT 1 FROM drugs d WHERE d.sid=s.sid)')
nt=c.execute('SELECT count(*) FROM treated_primary').fetchone()[0];missing=c.execute('SELECT count(*) FROM treated_primary WHERE height IS NULL').fetchone()[0];pat,n=c.execute('SELECT count(DISTINCT pid),count(*) FROM pairs WHERE is_primary AND bmi>40').fetchone()
after=before.copy()
for section,metric,v in [('missingness','height',missing),('missingness','height_fraction',missing/nt),('stratum_patients','bmi_gt40',pat),('stratum_pairs','bmi_gt40',n)]:after.loc[after.section.eq(section)&after.metric.eq(metric),'value']=v
changed=before.value.ne(after.value)&~(before.value.isna()&after.value.isna());assert set(before.loc[changed,'metric'])<={'height','height_fraction','bmi_gt40'}
assert pairguard==c.execute('SELECT count(*),sum(cuff_map),sum(art_map),sum(nee),sum(t) FROM pairs').fetchone()
after.to_csv(f,index=False);pd.DataFrame([{'corrected_field':'MIMIC admission height','height_missing_before':float(before.loc[before.metric.eq('height'),'value'].iloc[0]),'height_missing_after':missing,'denominator':nt,'pairs_values_times_doses_unchanged':True,'rule':'closest valid height linked to stay, including before ICU intake; inch converted to cm'}]).to_csv(OUT/'mimic_height_correction_qc.csv',index=False)
p=OUT/'deviation_log.csv'
with p.open(encoding='utf-8-sig',newline='') as z:r=csv.DictReader(z);fields=r.fieldnames;rows=list(r)
if not any(r['id']=='D007' for r in rows):
    rows.append(dict(id='D007',change='MIMIC入院身高改取同一stay关联的最近有效身高，允许记录早于ICU intime',reason='97.1%缺失触发源核查；绝大多数该stay身高记录发生在intime之前，先前后入ICU限制错误丢失入院身高',timing='首次全套聚合计数已生成后；身高红旗源核查后；R4正式判定之前；'+datetime.datetime.now().astimezone().isoformat(),authorization='提取纠错；保留纠正前CSV，患者/时间/血压/NEeq与信号样本不变',affected_outputs='MIMIC身高缺失率、BMI与BMI>40计数'))
    with p.open('w',encoding='utf-8-sig',newline='') as z:w=csv.DictWriter(z,fields);w.writeheader();w.writerows(rows)
j=ROOT/'sql/phase0_parameters.json';x=json.loads(j.read_text(encoding='utf-8'));x['implementation_details']['mimic_height_correction']='Closest valid height linked to the same ICU stay, including charttime before ICU intime; inches converted to cm. Correction after missingness red flag, logged D007.';j.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'height_missing_before':int(before.loc[before.metric.eq('height'),'value'].iloc[0]),'height_missing_after':missing,'denominator':nt,'pairs_unchanged':True}));c.close()
