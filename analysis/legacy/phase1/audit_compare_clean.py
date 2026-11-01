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
"""Compare source-reconstructed audit values with final phase-1 paired files.

This consumes actual source-query evidence. It does not call derived-pair
agreement an independent source query. Every per-pair comparison is restricted.
"""
from pathlib import Path
import sys,json,math,datetime,hashlib,traceback
ROOT=Path('@CUFF_WORKSPACE@')

import pandas as pd,duckdb
OUT=ROOT/'phase1';RES=OUT/'_restricted'
DBS=['mimic','eicu','amsterdam','sicdb']
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def matches(a,b):
    a=pd.to_numeric(a,errors='coerce');b=pd.to_numeric(b,errors='coerce')
    return (a.isna()&b.isna())|((a-b).abs()<=1e-6*(1+a.abs()))
def main():
    summaries=[];timeline=[]
    stamps=[]
    for p in OUT.glob('phase1_primary_*.csv'):
        st=p.stat();stamps.append(getattr(st,'st_birthtime',st.st_ctime))
    first_result=min(stamps) if stamps else None
    def relative(ts):
        if first_result is None:return '本轮主分层CSV尚未出现'
        v=datetime.datetime.fromisoformat(ts).timestamp()
        return '首份本轮主分层CSV创建之前' if v<first_result else '首份本轮主分层CSV创建之后'
    for db in DBS:
        manifest=RES/f'audit_summary_{db}.json'
        if not manifest.exists():
            summaries.append(dict(database=db,sampled_pairs=50,source_query_completed=False,source_discrepancy_pairs=None,clean_comparison_completed=False,clean_discrepancy_pairs=None));continue
        m=json.loads(manifest.read_text(encoding='utf-8'))
        evidence=RES/(f'audit_comparison_{db}.csv' if db in ['mimic','eicu'] else f'audit_file_comparison_{db}.csv')
        a=pd.read_csv(evidence,dtype={'sid':str})
        p=RES/f'pairs_{db}.parquet';has=p.exists();source_col={'mimic':('source_cuff_map','source_art_map'),'eicu':('source_cuff_map','source_art_map'),'amsterdam':('source_cuff','source_art'),'sicdb':('source_cuff','source_art')}[db]
        row=dict(database=db,sampled_pairs=len(a),source_query_completed=True,source_discrepancy_pairs=m['discrepancy_pairs'],source_clean_known_nee=int(a.source_clean_nee.notna().sum()),source_clean_unknown_nee=int(a.source_clean_nee.isna().sum()),clean_comparison_completed=has,clean_discrepancy_pairs=None)
        if has:
            before=digest(p)
            c=duckdb.connect();c.execute("SET threads=1;SET memory_limit='700MB'")
            c.register('audit_keys',a[['cid']]);z=c.execute(f"SELECT p.cid,p.sid,p.t,p.cuff_map,p.art_map,p.nee,p.weight,p.height,p.bmi FROM read_parquet('{p.as_posix()}') p JOIN audit_keys k USING(cid)").df();c.close()
            z=z.rename(columns={x:'production_'+x for x in z if x!='cid'})
            cmp=a.merge(z,on='cid',how='left',validate='one_to_one',indicator=True)
            fields={'cuff_map':source_col[0],'art_map':source_col[1],'nee':'source_clean_nee','weight':'source_clean_weight','height':'source_clean_height','bmi':'source_clean_bmi'}
            checkcols=[]
            for name,source in fields.items():
                col=name+'_clean_match';cmp[col]=matches(cmp[source],cmp['production_'+name]);checkcols.append(col)
            cmp['pair_identity_match']=cmp._merge.eq('both') & cmp.sid.astype(str).eq(cmp.production_sid.astype(str)) & ((cmp.t-cmp.production_t).abs()<.0001)
            checkcols.append('pair_identity_match');cmp['clean_all_match']=cmp[checkcols].all(axis=1)
            cmp.to_csv(RES/f'audit_clean_comparison_{db}.csv',index=False)
            if before!=digest(p):raise RuntimeError('Production paired file changed during bounded comparison')
            row.update(clean_discrepancy_pairs=int((~cmp.clean_all_match).sum()),matched_production_pairs=int(cmp._merge.eq('both').sum()),production_file_sha256=before,comparison_completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
            for col in checkcols:row[col+'_discrepancies']=int((~cmp[col]).sum())
            m.update(phase1_clean_output_comparison='completed',clean_discrepancy_pairs=row['clean_discrepancy_pairs'],clean_comparison_completed_at=row['comparison_completed_at'],production_file_sha256=before)
            manifest.write_text(json.dumps(m,ensure_ascii=False,indent=2),encoding='utf-8')
        summaries.append(row)
        for key in ['source_query_started_at','source_query_completed_at']:
            if key in m:timeline.append(dict(database=db,event=key,at_utc=m[key],relative_to_phase0='阶段0结果之后',relative_to_phase1_main=relative(m[key]),basis='源执行记录；分层CSV的文件创建时间，仅用于时序核对'))
        if row.get('comparison_completed_at'):timeline.append(dict(database=db,event='clean_comparison_completed',at_utc=row['comparison_completed_at'],relative_to_phase0='阶段0结果之后',relative_to_phase1_main=relative(row['comparison_completed_at']),basis='清理后生产文件与独立源重建比较'))
    pd.DataFrame(summaries).to_csv(OUT/'phase1_source_audit_summary.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(timeline).to_csv(OUT/'audit_execution_timeline.csv',index=False,encoding='utf-8-sig')
    lines=['# 阶段1独立源审计（内部）','',
      '固定 seed=20261027，从各库阶段0主配对中抽取50个动脉MAP<65的配对，共200个。采用带种子的SHA256排序抽样；全部避开阶段0旧源审计配对。抽样单位为配对，不保证每库50名不同患者。样本于阶段0结果之后、阶段1主分层结果之前固定；执行和完成时序见 audit_execution_timeline.csv。','',
      '| 数据库 | 样本配对 | 原始源回查不一致 | 清理后不一致 | 清理后NE-eq未知 |','| --- | --- | --- | --- | --- |']
    for r in summaries:
        lines.append('| '+r['database']+' | '+str(r['sampled_pairs'])+' | '+str(r.get('source_discrepancy_pairs') if r.get('source_query_completed') else '待完成')+' | '+str(r.get('clean_discrepancy_pairs') if r.get('clean_comparison_completed') else '待完成')+' | '+str(r.get('source_clean_unknown_nee','待完成'))+' |')
    lines+=['',
      'MIMIC/eICU 使用独立 PostgreSQL 源查询，临时表仅导入固定抽样与字典映射，随后在只读事务内读取真实 chartevents/inputevents/OMR 或 vitalperiodic/vitalaperiodic/infusiondrug/patient。Amsterdam 重新读取原始 numericitems 暂存副本及原始 admissions/drugitems；SICdb 重新读取原始 data_float_h、medication、cases 并独立按分钟解码。血压和药物源重建不导入生产转换函数，不以本地生产 pairs/drugs 代替源记录。阶段0 DuckDB仅用于固定抽样及期望值，不作为源验证。',
      '',
      '核查同刻同来源族 SBP/DBP/MAP 伪迹、原治疗episode内±5分钟动脉中位数、最近时间差、各组成药及其未记录gap、原单位及实际换算体重。SICdb克除以1000、AngII的ng/kg/min×0.01沿用已确认规则。eICU有界延续沿用已确认60分钟规则；无明确浓度的ml/hr保留未知。源单位或体重不足不得填零。',
      '',
      '本轮独立清理：体重30–300 kg、身高120–230 cm，重算BMI后保留12–80；界外置未知。药物换算实际需要体重且该源体重界外时剂量未知；直接每公斤质量速率及vasopressin U/min不因入院体重异常失效。eICU若原同刻已知中位数包含异常体重换算成分，该成分不能靠删除后重取中位数恢复为已知；本来不可换算的别名不参与该中位数。总NE-eq>5置未知。配对身份与血压保留。',
      '',
      '五项合理性检查的分工与依据：',
      '- 新增清理后的缺失/分布：phase1_quality_ranges.csv、phase1_covariate_missingness.csv；本审计仅对固定样本以源体重、身高和剂量复核。',
      '- 跨库事件/记录率：phase1_all_rates.csv、phase1_secondary_false_positive.csv、phase1_covariate_prevalence.csv；既有治疗/监测使用率参照阶段0报告。',
      '- 可判定比例：phase1_flow.csv、phase1_quality_ranges.csv、phase1_pair_coverage.csv；未知剂量单列。',
      '- 人数与配对守恒：phase1_flow.csv与原主配对数核对；本样本额外核对每库50个身份和血压保留。剂量层患者可跨层重叠。',
      '- 独立源数据：本报告及phase1_source_audit_summary.csv；沿用的阶段0缺失负例证据仅解释既有源覆盖，不冒充本轮新增查询。',
      '',
      '界限：本审计支持随机样本内的源转换与清理一致性，不证明全库无遗漏、全部概念映射已临床签核，也不验证选择偏倚、因果效应或袖带监测策略安全性。随机样本可能未覆盖每一种极端值；清理分支另以合成输入检查，不能把合成检查称为源数据验证。所有患者身份、逐条源证据及比较结果留在 phase1/_restricted。哈希用于证据完整性，不代替源验证。',
      '',
      '实施记录：复制阶段0独立审计到sql/phase1/audit_*；阶段0文件未覆盖。脚本准备时修正Python运行环境和字符串语法，均在新源查询前。阶段1分层结果产生前发现并纠正了vasopressin体重依赖语义；随后明确同刻别名中仅原本可换算的参与成分适用体重质量规则，固定样本及源导出不变。具体阶段1生产修正时序以主任务日志为准。']
    (OUT/'phase1_source_audit.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(summaries,ensure_ascii=False))
if __name__=='__main__':
    try:main()
    except Exception:
        (RES/'audit_clean_comparison_failure.log').write_text(traceback.format_exc(),encoding='utf-8')
        print(json.dumps({'status':'failed','details':'restricted diagnostic log; no row-level output'}));raise SystemExit(1)
