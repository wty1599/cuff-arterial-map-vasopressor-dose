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

from importlib import import_module
g=import_module('08_core');globals().update({k:v for k,v in vars(g).items() if not k.startswith('_')})
import shutil
rows=[]
for db in ['mimic','eicu','amsterdam','sicdb']:
    p=OUT/f'phase0_bp_quality_{db}.csv';b=OUT/'_logs'/f'phase0_bp_quality_{db}_before_disjoint_summary.csv'
    if not b.exists():shutil.copy2(p,b)
    c=connect(db)
    d=exported(c,p.name,"SELECT modality,count(*) candidate_time_family_rows,count(map_value) map_nonmissing,sum((map_value IS NULL)::INT) missing_map,sum((map_value IS NOT NULL AND NOT valid)::INT) observed_artifact_exclusions,sum((map_value IS NOT NULL AND incomplete_artifact_components)::INT) incomplete_artifact_components,sum(valid::INT) valid_rows FROM bp_quality GROUP BY modality")
    assert (d.candidate_time_family_rows==d.missing_map+d.observed_artifact_exclusions+d.valid_rows).all()
    for x in d.itertuples():rows.append({'db':db,'modality':x.modality,'candidate':x.candidate_time_family_rows,'missing_plus_artifact_plus_valid':x.missing_map+x.observed_artifact_exclusions+x.valid_rows,'conserved':True})
    c.close()
pd.DataFrame(rows).to_csv(OUT/'phase0_bp_summary_conservation.csv',index=False)
p=OUT/'deviation_log.csv'
with p.open(encoding='utf-8-sig',newline='') as f:r=csv.DictReader(f);fields=r.fieldnames;old=list(r)
if not any(x['id']=='D010' for x in old):
    old.append(dict(id='D010',change='血压质量汇总中伪迹剔除数只计MAP非缺失记录，使缺失/伪迹/有效三类互斥守恒',reason='原摘要NOT valid会把MAP缺失且SBP/DBP异常的记录也计入伪迹，造成摘要重复计数；实际清洗规则和配对未受影响',timing='全部配对和信号已生成后；报告汇总核查时；'+datetime.datetime.now().astimezone().isoformat(),authorization='聚合摘要纠错；不改数据、队列、阈值或任何配对结果',affected_outputs='phase0_bp_quality四库CSV；保留原摘要'))
    with p.open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fields);w.writeheader();w.writerows(old)
print(json.dumps({'disjoint_bp_summary_conserved':len(rows),'pairs_unchanged':True}))
