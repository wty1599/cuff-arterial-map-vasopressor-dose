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

"""Bounded post-result care-level check; not a redefinition of the primary analysis."""
from importlib import import_module
g=import_module('08_core');globals().update({k:v for k,v in vars(g).items() if not k.startswith('_')})
u=pd.read_csv(OUT/'mimic_transfer_careunit_dictionary.csv',keep_default_na=False)
u['classification']=np.where(u.careunit.str.contains(r'Intensive Care|\bCCU\b|\bTSICU\b',case=False,regex=True),'ICU',np.where(u.careunit.str.contains('intermediate|stepdown',case=False),'MC',np.where(u.careunit.str.contains(r'^unknown$|^$',case=False),'unknown','other')))
u['human_review']='pending';u.to_csv(OUT/'map_careunit_mimic_transfers.csv',index=False)
c=connect('mimic');save(c,'care_probe_map',u)
q=csvq(R/'mimic_transfer_careunit_probe.csv')
c.execute(f"""CREATE OR REPLACE TEMP TABLE care_probe_intervals AS
 SELECT s.sid,s.pid,m.classification,greatest(s.t0,epoch(try_cast(x.intime AS TIMESTAMP))) start_time,
 least(s.t1,epoch(try_cast(x.outtime AS TIMESTAMP))) finish
 FROM stays s JOIN {q} x ON s.hadm=x.hadm_id JOIN care_probe_map m ON x.careunit=m.careunit
 WHERE epoch(try_cast(x.outtime AS TIMESTAMP))>s.t0 AND epoch(try_cast(x.intime AS TIMESTAMP))<s.t1
 AND try_cast(x.outtime AS TIMESTAMP)>try_cast(x.intime AS TIMESTAMP)""")
c.execute("""CREATE OR REPLACE TEMP TABLE care_probe_first AS WITH x AS (
 SELECT sid,pid,min(start_time) first_icu_time FROM care_probe_intervals WHERE classification='ICU' GROUP BY sid,pid)
 SELECT *,row_number() OVER(PARTITION BY pid ORDER BY first_icu_time,sid) rn FROM x""")
c.execute("""CREATE OR REPLACE TEMP TABLE care_probe_anchors AS SELECT cf.cid,cf.sid,s.pid,s.is_primary,
 EXISTS(SELECT 1 FROM care_probe_first f WHERE f.sid=cf.sid AND f.rn=1) first_confirmed,
 EXISTS(SELECT 1 FROM care_probe_intervals t WHERE t.sid=cf.sid AND t.classification='ICU' AND cf.t>=t.start_time AND cf.t<t.finish) has_icu,
 EXISTS(SELECT 1 FROM care_probe_intervals t WHERE t.sid=cf.sid AND t.classification<>'ICU' AND cf.t>=t.start_time AND cf.t<t.finish) overlap_other_or_unknown,
 EXISTS(SELECT 1 FROM pairs_windows p WHERE p.cid=cf.cid AND p.win=5) any_existing_pair
 FROM cuff60 cf JOIN stays s USING(sid) WHERE s.adult""")
rows=[]
for scope,cond in [('current_source_first','is_primary'),('first_with_confirmed_icu_transfer','first_confirmed')]:
    z=c.execute(f"SELECT count(*),count(*) FILTER(WHERE has_icu),count(*) FILTER(WHERE has_icu AND overlap_other_or_unknown),count(*) FILTER(WHERE has_icu AND NOT overlap_other_or_unknown),count(*) FILTER(WHERE has_icu AND NOT overlap_other_or_unknown AND any_existing_pair),count(DISTINCT pid) FILTER(WHERE has_icu AND NOT overlap_other_or_unknown) FROM care_probe_anchors WHERE {cond}").fetchone()
    row=dict(zip(['all_cuff_anchors','anchors_with_icu','icu_overlap_other_or_unknown','unambiguous_icu_anchors','existing_paired_icu_anchors_upper_bound','patients_with_unambiguous_icu_anchors'],z));row['scope']=scope;row['pair_fraction_upper_bound']=z[4]/z[3] if z[3] else np.nan;row['interpretation']='numerator still permits arterial values outside strict ICU; actual restricted-window pair fraction cannot exceed this upper bound for this anchor set';rows.append(row)
pd.DataFrame(rows).to_csv(OUT/'phase0_mimic_careunit_bound.csv',index=False)
c.close()
p=OUT/'deviation_log.csv'
with p.open(encoding='utf-8-sig',newline='') as f:r=csv.DictReader(f);fields=r.fieldnames;old=list(r)
if not any(x['id']=='D011' for x in old):
    old.append(dict(id='D011',change='纠正报告将四库均称严格ICU的措辞；MIMIC主分析为来源icustays口径，另独立查询transfers作限定的R3上界核查',reason='末次单元名称复核发现MIMIC first_careunit含中间监护/stepdown/PACU；起始名称不能代表整段，不能直接整例删除',timing='四库配对、初步R1-R4、初稿和第一次正文独立复核之后；交付前；'+datetime.datetime.now().astimezone().isoformat(),authorization='纠正来源语义与评估是否可能反转R3；不改原主分析分母、阈值、样本或一致性',affected_outputs='报告照护级别表述；额外源转科映射与配对覆盖上界；原有配对/信号/分层CSV不变'))
    with p.open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fields);w.writeheader();w.writerows(old)
j=ROOT/'sql/phase0_parameters.json';p=json.loads(j.read_text(encoding='utf-8'));p['implementation_details']['mimic_care_scope']='Primary uses source icustays as provided, not proven pure ICU at every time. Post-result bounded transfers probe recorded D011; original main counts unchanged.';j.write_text(json.dumps(p,ensure_ascii=False,indent=2),encoding='utf-8')
print(pd.DataFrame(rows).to_json(orient='records'))
