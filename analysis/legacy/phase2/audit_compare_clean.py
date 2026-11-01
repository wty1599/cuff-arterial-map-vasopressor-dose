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
"""Compare independent raw-source evidence with read-only phase-1 clean pairs."""
from pathlib import Path
import sys,json,datetime,hashlib,traceback
ROOT=Path('@CUFF_WORKSPACE@');
import pandas as pd,duckdb
OUT=ROOT/'phase2';R=OUT/'_restricted';DBS=['mimic','eicu','amsterdam','sicdb']
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def matches(a,b):
    a=pd.to_numeric(a,errors='coerce');b=pd.to_numeric(b,errors='coerce')
    return (a.isna()&b.isna())|((a-b).abs()<=1e-6*(1+a.abs()))
def main():
    samples=pd.read_csv(OUT/'audit_sample_manifest.csv')
    artifacts=pd.read_csv(OUT/'audit_artifact_summary.csv').set_index('database').to_dict('index')
    stamps=[getattr(p.stat(),'st_birthtime',p.stat().st_ctime) for p in OUT.glob('phase2_A_*.csv')]
    first_result=min(stamps) if stamps else None
    def relative(ts):
        if first_result is None:return '阶段2结果文件尚未出现'
        return '首份阶段2 A结果文件创建之前' if datetime.datetime.fromisoformat(ts).timestamp()<first_result else '首份阶段2 A结果文件创建之后'
    summaries=[];timeline=[]
    for db in DBS:
        sample=samples[samples.database.eq(db)].iloc[0]
        m=json.loads((R/f'audit_summary_{db}.json').read_text(encoding='utf-8'))
        assert m['source_queries_executed'] is True
        a=pd.read_csv(R/(f'audit_comparison_{db}.csv' if db in ['mimic','eicu'] else f'audit_file_comparison_{db}.csv'),dtype={'sid':str})
        if len(a)!=int(sample.sampled_pairs):raise RuntimeError('Evidence/sample count mismatch')
        p=ROOT/'phase1'/'_restricted'/f'pairs_{db}.parquet';before=digest(p)
        if before!=sample.phase1_input_sha256:raise RuntimeError('Protected phase1 input changed')
        c=duckdb.connect();c.execute("SET threads=1;SET memory_limit='700MB'")
        c.register('audit_keys',a[['cid']]);z=c.execute(f"SELECT p.cid,p.sid,p.t,p.cuff_map,p.art_map,p.nee,p.weight,p.height,p.bmi FROM read_parquet('{p.as_posix()}') p JOIN audit_keys k USING(cid)").df();c.close()
        z=z.rename(columns={x:'production_'+x for x in z if x!='cid'})
        cmp=a.merge(z,on='cid',how='left',validate='one_to_one',indicator=True)
        cuff,art=('source_cuff_map','source_art_map') if db in ['mimic','eicu'] else ('source_cuff','source_art')
        fields={'cuff_map':cuff,'art_map':art,'nee':'source_clean_nee','weight':'source_clean_weight','height':'source_clean_height','bmi':'source_clean_bmi'}
        flags=[]
        for name,source in fields.items():
            col=name+'_clean_match';cmp[col]=matches(cmp[source],cmp['production_'+name]);flags.append(col)
        cmp['pair_identity_match']=cmp._merge.eq('both') & cmp.sid.astype(str).eq(cmp.production_sid.astype(str)) & ((cmp.t-cmp.production_t).abs()<.0001)
        cmp['deep_missed_qualification_match']=cmp[art].lt(55)&cmp[cuff].ge(65)
        flags+=['pair_identity_match','deep_missed_qualification_match'];cmp['clean_all_match']=cmp[flags].all(axis=1)
        cmp.to_csv(R/f'audit_clean_comparison_{db}.csv',index=False)
        assert before==digest(p)
        end=datetime.datetime.now(datetime.timezone.utc).isoformat()
        row=dict(database=db,sampled_pairs=len(a),eligible_deep_missed_pairs=int(sample.eligible_deep_missed_pairs),
            shortfall=int(sample.shortfall),overlap_phase0_phase1_audit_pairs=int(sample.overlap_old_sample_pairs),
            source_query_completed=True,source_discrepancy_pairs=m['discrepancy_pairs'],
            clean_comparison_completed=True,clean_discrepancy_pairs=int((~cmp.clean_all_match).sum()),
            matched_production_pairs=int(cmp._merge.eq('both').sum()),source_clean_known_nee=int(a.source_clean_nee.notna().sum()),
            source_clean_unknown_nee=int(a.source_clean_nee.isna().sum()),production_file_sha256=before,comparison_completed_at=end)
        for f in flags:row[f+'_discrepancies']=int((~cmp[f]).sum())
        for k,v in artifacts[db].items():
            if k.startswith('artifact_') or k=='quality_records':row[k]=v
        summaries.append(row)
        m.update(phase2_clean_output_comparison='completed',clean_discrepancy_pairs=row['clean_discrepancy_pairs'],clean_comparison_completed_at=end,production_file_sha256=before)
        (R/f'audit_summary_{db}.json').write_text(json.dumps(m,ensure_ascii=False,indent=2),encoding='utf-8')
        events={'sample_fixed':sample.created_at,'source_query_started_at':m['source_query_started_at'],'source_query_completed_at':m['source_query_completed_at'],'clean_comparison_completed':end}
        for event,at in events.items():timeline.append(dict(database=db,event=event,at_utc=at,relative_to_phase1='已见阶段1结果之后',relative_to_phase2=relative(at),basis='执行时间；首份阶段2 A聚合CSV创建时间，仅说明时序'))
    pd.DataFrame(summaries).to_csv(OUT/'phase2_source_audit_summary.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(timeline).to_csv(OUT/'audit_execution_timeline.csv',index=False,encoding='utf-8-sig')
    lines=['# 阶段2独立源审计（内部）','',
      '按 seed=20261027，从只读阶段1主配对中每库抽取30个动脉MAP<55且袖带MAP≥65配对，共120对；所有库候选均足够，全部避开阶段0/1旧审计配对。配对是抽样单位，不保证每库30名不同患者。抽样及源查询均在已见阶段1结果后；相对阶段2结果时序见audit_execution_timeline.csv。','',
      '| 库 | 样本 | 源重建不一致 | 清理后不一致 | NE-eq未知 | suspected | none_supported | unassessable |',
      '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for r in summaries:
        n=r['sampled_pairs'];pct=lambda k:f"{int(r[k])}/{n} ({r[k]/n:.1%})"
        lines.append(f"| {r['database']} | {n} | {r['source_discrepancy_pairs']} | {r['clean_discrepancy_pairs']} | {r['source_clean_unknown_nee']} | {pct('artifact_suspected')} | {pct('artifact_none_supported')} | {pct('artifact_unassessable')} |")
    lines+=['',
      '源证据：MIMIC/eICU在临时表导入样本/映射后，以只读PostgreSQL事务查询真实chartevents/inputevents/OMR、vitalperiodic/vitalaperiodic/infusiondrug/patient。Amsterdam重新扫描原始numericitems暂存副本及原始admissions/drugitems；SICdb重新扫描data_float_h、medication、cases，独立按分钟解码。仅从派生表取得期望值和抽样身份；源值不由生产pairs/drugs重算，也不导入生产转换函数。',
      '',
      '回查同刻同来源族SBP/DBP/MAP伪迹、原治疗episode内±5分钟动脉中位数、最近动脉时间及候选数量、源药物组成、单位、每个成分gap未知与实际换算体重。继续使用已确认的AngII ng/kg/min×0.01、eICU最长60分钟有界延续和SICdb入院体重克÷1000。无明确浓度的ml/hr不猜剂量。清理沿用phase1体重30–300 kg、身高120–230 cm、BMI12–80及NE-eq>5置未知；配对保留。',
      '',
      '阻尼/冲洗界限：现有字典和固定样本±5分钟记录中的明确动脉波形异常才标suspected；明确正常波形且无异常记录标none_supported，仅表示当前记录未支持异常，不能排除真实伪迹。缺少相关质量记录、只有普通冲洗液或校零操作、仅有低art/高cuff，以及其他器官/IABP波形均不能据此裁定；记unassessable。Amsterdam/SICdb本次可用字典无可等同的外周动脉质量通道，邻近数值记录不足以评估波形。未新增波形重建或临床裁定。',
      'MIMIC的Positional仅作具名“位置依赖动脉波形疑点”，不称已确证阻尼、冲洗或当前MAP错误。该分类在已看到阶段2初步结果、审计读到该源标签之后明确。Sharp不自动认作正常或异常；普通校零也不认作伪迹。',
      '',
      'suspected的受限抽样cid清单位于phase2/_restricted/audit_suspected_artifact_cids.csv。抽样中源支持的疑点不能直接外推为全库伪迹率，不能只删除抽样阳性后称全库伪迹已排除。本审计不改主分析。',
      '',
      '五项QA分工：',
      '1. Table1新增变量缺失/分布引用已完成的phase2_table1.csv、phase2_table1_notes.md及table1_quality_checks.csv；独立基线源证据见table1_source_audit.csv。本审计未重复该分支，只以源记录核对固定样本体格及剂量。',
      '2. 跨库事件率由phase2_A_deep_miss.csv、phase2_A_miss_distribution.csv、phase2_B_cuff_perspective.csv提供，分母分别按明确目标事件定义；本审计不重新估计总体事件率。',
      '3. 总体可判定比例及有效bootstrap重复由主任务表列出，原未知剂量继续单列；本报告另列样本剂量未知及伪迹证据三状态，未知不当作无伪迹。',
      '4. 人群/剂量层计数守恒参照主任务phase2_validation.csv；本审计单独核身份、样本资格及清理保留。患者可跨剂量层，不能直接加总患者数。',
      '5. 新独立源回查见本报告和phase2_source_audit_summary.csv。阶段0负例审计仅作为未改提取范围的既有证据，不冒称本轮重做。',
      '',
      '实施与局限：复制既有独立审计至sql/phase2/audit_*后最小适配，新样本固定前写定深度漏判范围；额外保存两文件库最近时间候选和动脉记录数量，增强时间核对。未修改阶段0/1文件。所有逐条证据仅在phase2/_restricted。哈希用于确认只读输入版本，不代替源验证；样本内一致不代表全库零误差、临床映射已签核或袖带策略安全。']
    lines+=['','源查询期间主任务仅更新阶段2 C分析子流/CI说明，未改变抽样或清理语义。审计的全参数文件守卫因此被改为只核真正相关定义，并保留完整旧新指纹；eICU有界源查询重执行一次。此为检查程序修正，不是源数据不一致，时序见audit_implementation_note.json。']
    coverage=OUT/'audit_mimic_positional_coverage.csv'
    if coverage.exists():
        v=pd.read_csv(coverage).iloc[0]
        lines+=['',f"追加授权的MIMIC全主配对查询：以精确标签Arterial line Waveform Appear、值Positional和±5分钟窗口，只读回查全部{int(v.primary_pairs):,}对。{int(v.pairs_with_any_waveform_record):,}对（{v.waveform_record_coverage:.2%}）有该波形字段；{int(v.pairs_with_positional_record):,}对、{int(v.patients_with_positional_flag):,}名患者带Positional标记，占全部配对{v.positional_flag_fraction:.2%}。{int(v.pairs_without_quality_record):,}对无该窗口内质量记录，不当作正常。完整受限清单：phase2/_restricted/audit_suspected_artifact_cids_all_pairs.csv；聚合依据：audit_mimic_positional_coverage.csv。主任务据此执行具名探索敏感性，不称已排除所有阻尼/冲洗。该追加查询在已见阶段2结果、发现抽样Positional之后获授权并执行，开始/完成时间保存在覆盖CSV。"]
    (OUT/'phase2_source_audit.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(summaries,ensure_ascii=False))
if __name__=='__main__':
    try:main()
    except Exception:
        (R/'audit_clean_comparison_error.log').write_text(traceback.format_exc(),encoding='utf-8')
        print(json.dumps({'status':'failed','details':'restricted audit_clean_comparison_error.log'}));raise SystemExit(1)
