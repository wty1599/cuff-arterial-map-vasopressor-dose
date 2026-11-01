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

"""Bounded source-unit correction after first aggregates; BP pairs and sample preserved."""
from importlib import import_module
g=import_module('08_core');globals().update({k:v for k,v in vars(g).items() if not k.startswith('_')})
import shutil
from pathlib import Path
db='sicdb';c=connect(db)
guard_sql='SELECT count(*),bit_xor(hash(cid,sid,t,cuff_map,art_map,nearest_diff,nearest_art_t,art_min_t,art_max_t,episode_start,episode_stop)) FROM (SELECT p.*,cf.episode_start,cf.episode_stop FROM pairs p JOIN cuff60 cf USING(cid))'
before_guard=c.execute(guard_sql).fetchone()
logdir=OUT/'_logs/sicdb_before_unit_af_correction';logdir.mkdir(exist_ok=True)
old=pd.read_csv(logdir/f'phase0_feasibility_{db}.csv' if (logdir/f'phase0_feasibility_{db}.csv').exists() else OUT/f'phase0_feasibility_{db}.csv')
for f in OUT.glob('*sicdb.csv'):
    dest=logdir/f.name
    if not dest.exists():shutil.copy2(f,dest)
for name in ['audit_input_sicdb.csv','audit_file_comparison_sicdb.csv']:
    src=R/name;dst=R/(src.stem+'_before_unit_af_correction.csv')
    if src.exists() and not dst.exists():shutil.copy2(src,dst)
c.close()
import_module('08_build_stays_drugs').build(db)
c=connect(db)
c.execute("""CREATE OR REPLACE TABLE pair_actual AS SELECT p.cid,count(DISTINCT d.drug) n_active,
 sum(d.dose) sum_dose,count(*) FILTER(WHERE d.drug IS NOT NULL AND d.dose IS NULL) n_unknown
 FROM pair_base p LEFT JOIN drugs d ON p.sid=d.sid AND p.t>=d.start AND p.t<d.end GROUP BY p.cid""")
c.execute("""UPDATE pairs SET nee=CASE WHEN e.n_expected>0 AND e.n_expected=a.n_active AND a.n_unknown=0 THEN a.sum_dose END,
 nee_status=CASE WHEN e.n_expected=0 OR e.n_expected<>a.n_active THEN 'undocumented_gap_component' WHEN a.n_unknown>0 THEN 'unit_weight_or_rate_unknown' ELSE 'known' END
 FROM pair_actual a JOIN pair_expected e USING(cid) WHERE pairs.cid=a.cid""")
c.execute('UPDATE pairs SET bmi=s.bmi,af=s.af,af_source=s.af_source FROM stays s WHERE pairs.sid=s.sid')
c.execute('CREATE OR REPLACE TABLE treated_primary AS SELECT * FROM stays s WHERE s.is_primary AND EXISTS(SELECT 1 FROM drugs d WHERE d.sid=s.sid)')
assert before_guard==c.execute(guard_sql).fetchone()
updated=old.copy()
strata=[('le0.1','nee<=0.1'),('gt0.1_le0.25','nee>0.1 AND nee<=0.25'),('gt0.25_le0.5','nee>0.25 AND nee<=0.5'),('gt0.5_le1','nee>0.5 AND nee<=1'),('gt1','nee>1'),('gt0.5_total','nee>0.5'),('unknown','nee IS NULL'),('bmi_gt40','bmi>40'),('af_any','af')]
for label,cond in strata:
    pat,n=c.execute(f'SELECT count(DISTINCT pid),count(*) FROM pairs WHERE is_primary AND ({cond})').fetchone()
    for sec,v in [('stratum_patients',pat),('stratum_pairs',n)]:updated.loc[updated.section.eq(sec)&updated.metric.eq(label),'value']=v
updated.to_csv(OUT/f'phase0_feasibility_{db}.csv',index=False)
exported(c,'phase0_af_sources_sicdb.csv','SELECT af_source,count(DISTINCT pid) patients,count(*) pairs FROM pairs WHERE is_primary AND af GROUP BY af_source')
exported(c,'phase0_nee_unknown_reasons_sicdb.csv','SELECT nee_status,count(*) pairs,count(DISTINCT pid) patients FROM pairs WHERE is_primary GROUP BY nee_status')
# Refresh expected values for the exact saved random audit sample; do not resample.
sample=pd.read_csv(R/'audit_input_sicdb.csv',dtype={'sid':str});ids=sample[['audit_id','sid','t']].copy()
actual=c.execute('SELECT sid,t,nee FROM pairs').df();assert not actual.duplicated(['sid','t']).any()
look=actual.set_index(['sid','t']).nee
sample['nee']=look.reindex(pd.MultiIndex.from_frame(sample[['sid','t']])).to_numpy();assert ids.equals(sample[['audit_id','sid','t']]);sample.to_csv(R/'audit_input_sicdb.csv',index=False)
changed=old.value.ne(updated.value)&~(old.value.isna()&updated.value.isna())
diff=old.loc[changed,['section','metric','value']].rename(columns={'value':'before'}).merge(updated[['section','metric','value']],on=['section','metric']).rename(columns={'value':'after'})
diff['bp_pairs_times_patients_unchanged']=True;diff.to_csv(OUT/'sicdb_units_af_correction_qc.csv',index=False)
c.close();import_module('16_complete_aggregates').finish(db)
path=OUT/'deviation_log.csv'
with path.open(encoding='utf-8-sig',newline='') as f:r=csv.DictReader(f);fields=r.fieldnames;rows=list(r)
for rec in [
 dict(id='D008',change='SICdb WeightOnAdmission 先从克除以1000转换为kg；从原始用药重新换算NE-eq并更新BMI',reason='最终跨库合理性核查发现BMI几乎均超40、体重中位数75000；本地既往源代码和报告确认发布值为克，官方5kg表述为分箱宽度；初次源审计也复制了单位误读',timing='初步SICdb配对、剂量/BMI计数、总体信号及初次源审计之后；正式R2/R4结论之前；'+datetime.datetime.now().astimezone().isoformat(),authorization='数据单位纠错；保留初版聚合与源审计；不改研究阈值',affected_outputs='SICdb剂量分层、BMI>40；同一100例源审计重新核对；患者/血压/时间/信号样本不变'),
 dict(id='D009',change='SICdb房颤仅纳入明确I48.0/.1/.2或不含房扑的明确房颤主诊断文本；I48/I48.9不再当确证AF',reason='I48类包含房扑及未分型房颤/房扑，原宽泛匹配会高估AF',timing='初步房颤计数和总体信号之后；最终内部报告之前；'+datetime.datetime.now().astimezone().isoformat(),authorization='概念匹配纠错，未计算房颤分层一致性',affected_outputs='SICdb房颤患者/配对计数；不改血压、队列和信号')]:
    if not any(x['id']==rec['id'] for x in rows):rows.append(rec)
with path.open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fields);w.writeheader();w.writerows(rows)
q=ROOT/'sql/phase0_parameters.json';p=json.loads(q.read_text(encoding='utf-8'));p['implementation_details']['sicdb_weight']='Released WeightOnAdmission is grams; divide by 1000 to kg before BMI and mass/kg dosing. Corrected after preliminary results, D008.';p['implementation_details']['sicdb_af']='Main diagnosis only; explicit AF I48.0/.1/.2 or unambiguous AF text. I48 and I48.9 represent AF/flutter ambiguity and are not confirmed AF.';p['implementation_details']['first_stay_sicdb']=p['implementation_details']['first_stay_sicdb'].replace('by OffsetAfterFirstAdmission is','by OffsetAfterFirstAdmission + ICUOffset is');q.write_text(json.dumps(p,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'db':db,'correction_completed':True,'blood_pressure_pairs_unchanged':True,'audit_sample_unchanged':True,'changed_aggregate_cells':len(diff)}))
