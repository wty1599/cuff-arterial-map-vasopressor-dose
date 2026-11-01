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
import pandas as pd,numpy as np,json,hashlib
ROOT=Path('@CUFF_WORKSPACE@');O=ROOT/'phase2';R=O/'_restricted';s=pd.read_csv(O/'phase2_table1.csv');flow=pd.read_csv(O/'phase2_flow_numbers.csv');au=pd.read_csv(O/'table1_source_audit.csv');q=pd.read_csv(O/'table1_quality_checks.csv')
names={'mimic':'MIMIC-IV','eicu':'eICU','amsterdam':'AmsterdamUMCdb','sicdb':'SICdb'}
errors=[]
for db in names:
 t=s[s.database.eq(db)];N=int(t.iloc[0].patients)
 if not ((t.n_observed+t.n_missing)==N).all():errors.append(db+': observed+missing')
 if not np.allclose(t.missing_fraction,t.n_missing/N):errors.append(db+': missing fraction')
 zz=t[t.percent.notna()]
 if not np.allclose(zz.percent,zz.numerator/zz.denominator*100):errors.append(db+': percentage')
 pf=flow[(flow.database==db)&flow.stage.str.startswith('phase1_dose_')]
 total=int(flow[(flow.database==db)&(flow.stage=='paired5_readings')].value.iloc[0])
 if pf.value.sum()!=total:errors.append(db+': flow dose sum')
 for v in ['sex','hospital_death','atrial_fibrillation','bmi_gt40','dose_observation_status','unit_type','admission_type','admission_source','age_group']:
  z=t[t.variable.eq(v)]
  if len(z) and z.numerator.notna().any() and not np.isclose(z.numerator.sum(),z.n_observed.iloc[0]):errors.append(db+': '+v+' categorical sum')
print(json.dumps({'aggregate_checks_errors':errors,'table1_rows':len(s),'flow_rows':len(flow),'source_audit_patients':int(au.patients.sum()),'source_peak_points':int(au.peak_points_checked.sum()),'source_discrepancies':int(au.baseline_field_discrepancies.sum()+au.peak_point_discrepancies.sum())}))
assert not errors
def line(db,var,category=None):
 z=s[(s.database==db)&(s.variable==var)];return z[z.category.eq(category)].iloc[0] if category is not None else z.iloc[0]
def fi(x):return f'{int(x):,}'
def med(r,scale=1):return '不可得' if pd.isna(r['median']) else f"{r['median']*scale:.3g} [{r.q1*scale:.3g}, {r.q3*scale:.3g}]"
def table(head,rows):return '| '+' | '.join(head)+' |\n|'+'|'.join(['---']*len(head))+'|\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in rows)
parts=['# 阶段 2 Table 1 与四库流程说明','内部中文，映射待签核版。仅描述基线，不做死亡或评分关联。阶段 0/1 文件均未写入。','## 人群与字段口径','主体为阶段 1 ±5 min 配对患者，每患者一个索引 ICU stay、一个基线行。分类百分比以该字段可判定患者为分母；patients 始终保留全体基线患者数，n_missing 与 missing_fraction 另列。药物使用分类可重叠，不能相加为总人数。CSV 中 n_observed=0 的二分类行其 0/0 不代表事件率为 0，percent 留空并由 status 解释。','长表：[phase2_table1.csv](phase2_table1.csv)；流程：[phase2_flow_numbers.csv](phase2_flow_numbers.csv)。']
rows=[]
for db in names:
 r=line(db,'age');rows.append([names[db],fi(r.patients),med(r),med(line(db,'weight')),med(line(db,'bmi')),med(line(db,'pairs_per_patient'))])
parts+= [table(['数据库','配对患者','年龄：中位数 [IQR]','体重 kg','BMI','每患者配对数'],rows),'Amsterdam 的年龄、体重和 BMI 数值是闭区间中点代理；开放区间未知。原年龄区间另有 age_group 行，不能把中点解释为精确年龄。SICdb 年龄及体格保留原数据粗化语义，体重从源文件克换算为 kg。MIMIC/eICU 年龄保留源脱敏或顶端编码。','## 最大剂量与可判定性','本表的变量完整名称为**可判定治疗片段中观测到的最大 NE-eq**。计算覆盖索引 stay 的所有既有升压药治疗片段，使用药物、预期药物 episode 和总治疗 episode 的全部起止点分割时间，合并同时活动成分；没有限定在袖带或动脉读数时点。任一预期药处于无剂量间隙、活动成分未知、实际换算体重异常或合计 >5，则总剂量未知。没有插补或将未知当零，也没有并入恢复版 eICU 剂量。','known_duration_fraction 的分母为当前定义治疗 episodes 总时长，分子为总剂量可判定时长。complete 仅表示这些已记录/推定 episodes 内可判定，不证明医院所有未记录治疗都被完整观察。eICU 继续采用正剂量至下一同药记录、最长 60 分钟的推定治疗期。']
rows=[]
for db in names:
 r=line(db,'max_nee_observed');obs=line(db,'known_duration_fraction');zz=s[(s.database==db)&(s.variable=='dose_observation_status')].set_index('category').numerator.to_dict();rows.append([names[db],f'{fi(r.n_observed)}/{fi(r.patients)}',med(r),med(obs,100),fi(zz.get('complete_within_defined_episodes',0)),fi(zz.get('partial_observation',0)),fi(zz.get('no_known_segment',0))])
parts +=[table(['数据库','可得最大值患者','最大 NE-eq：中位数 [IQR]','患者可判定时长 %：中位数 [IQR]','定义内完整','部分观察','无已知片段'],rows),'eICU 许多患者没有任何可判定总剂量片段，所给最大值分布仅描述可得患者；不能把该分布外推为全部治疗患者的真实最高剂量。所有患者仍保留在 Table 1。','## 评分、死亡和不可得原因']
rows=[]
for db in names:
 se=line(db,'severity_score');de=line(db,'hospital_death','yes');rows.append([names[db],se.label,med(se),f'{fi(se.n_missing)}/{fi(se.patients)}',('不可得' if not de.n_observed else f'{fi(de.numerator)}/{fi(de.denominator)} ({de.percent:.1f}%)'),f'{fi(de.n_missing)}/{fi(de.patients)}'])
parts +=[table(['数据库','评分','评分：中位数 [IQR]','评分缺失','院内死亡/可判定','死亡状态缺失'],rows),
'- MIMIC 使用现成 first_day_sofa，没有重建评分。eICU 仅使用 APACHE IVa，排除缺失/负值，不与 IV 混用；同版本冲突值置未知。缺评分表示源表无可用值，不根据其他变量补算。',
'- Amsterdam 源 admissions 没有完整医院出院状态，dateofdeath 及 ICU 去向不能替代院内死亡，因此本项全部不可得；本轮也没有对应的现成 SOFA/APACHE IVa 字段。',
'- SICdb 用 HospitalDischargeType 与本版字典确认医院死亡/活着离院；不使用 ICU DischargeState 或含院外死亡的 OffsetOfDeath。院内转科结案和不明结案不归为活着离院。本样本医院结案均可判定。源 cases 虽有 SAPS3，但未替换本轮指定的 SOFA/APACHE IVa，避免把不同量表混为同列。',
'- 房颤保留阶段 1 定义：MIMIC 心律/诊断覆盖，eICU 诊断，Amsterdam 心律；SICdb 不做房颤判定。Amsterdam 不从区间 BMI 推断 BMI>40。',
'- ICU 单元与入院类别按各源含义列示。MIMIC admission_type 是住院入院类别，eICU unitstaytype 是 ICU stay 类型，Amsterdam urgency 是计划/非计划，SICdb 是 SurgicalAdmissionType；不可直接当作一套统一内外科分组。SIC HospitalUnit 是病例最后单元，不标为首入单元。',
'- 体重/BMI 等缺失沿用阶段 1 合理范围处理；详细逐字段数量见长表 n_missing。未发现需要超过两小时的新评分或结局提取工作，因此没有新建评分流水线。',
'依据：[eICU APACHE 结果表](https://eicu.mit.edu/eicutables/apachepatientresult/)；[SICdb cases 字段说明](https://www.sicdb.com/Documentation/Table:_Cases)，2026-10-27 核对。',
'## 五项 QA 与独立源查',
'1. 每项缺失率及中位数/分布均保存在长表。年龄区间代理、SIC 性别新提取和源体重单位均有明确状态；未把阶段 1 未提取等同源不存在。',
'2. 药物使用比例和院内死亡仅描述各库，不做模型或结局关联；死亡来源的不可观察状态单列。',
'3. 剂量观察完整性分为定义内完整、部分观察和无已知片段，已知时长与未知时长相加守恒。',
'4. 每库患者/索引 stay 唯一，每患者配对数之和等于阶段 1 配对数；五剂量层加未知层等于所有配对。分类观测数加缺失数等于患者数，百分比分母重新核对。流程同时保留患者、stay、读数和配对单位；首次顺序未知等旁支与其他排除不作互斥相加。',
f'5. 每库抽 12 名患者，共 {fi(au.patients.sum())} 名。MIMIC/eICU 另写 PostgreSQL 查询，Amsterdam/SICdb 另用 DuckDB SQL 读取源 CSV；独立解析性别、医院死亡、评分和入院类别，Amsterdam 还核对年龄组中点。共 {fi(au.baseline_field_checks.sum())} 个基线字段比较，差异 {fi(au.baseline_field_discrepancies.sum())}；在 {fi(au.peak_points_checked.sum())} 个已报告最大值时点，用独立单位公式回算源活动药物总量，差异 {fi(au.peak_point_discrepancies.sum())}。该检查确认所报峰值时点的源剂量，并非独立重建每例整个轨迹的最大值。行级证据留在受限目录。',
'源查摘要：[table1_source_audit.csv](table1_source_audit.csv)；计数校验：[table1_quality_checks.csv](table1_quality_checks.csv)。',
'## 展示适配与时序',
'本次规则写定发生在已看到阶段 1 主结果、核对本次源字段之后，Table 1 分布与治疗片段最大值计算之前。随后见到 Table 1 后执行有界独立源查，未据此调整人群或最大值定义。',
'Amsterdam 的阶段 1 年龄字段保留年龄组下界用于既有资格；本轮 Table 1 为避免把下界误当实际年龄，单独展示闭区间中点并保留原始组别，开放区间未知。该适配仅作用于 Table 1，没有改成人资格、配对、阶段 1/2 主分析或模型，已记入 D035。',
 '一次技术修复：本分支最初直接调用 pandas parquet 引擎而未加载本地依赖，改用已有 DuckDB 只读 parquet 后成功；未改变任何数值定义。',
'需 Han 判断但不阻塞：映射人工签核；投稿时是否保留数据库特有入院分类；若需要跨库同一严重程度量表或 Amsterdam 完整院内死亡，需要另行定义与提取，不能从本表现有字段补造。']
idx=parts.index('## 五项 QA 与独立源查')+1
missrows=[];drrows=[]
for db in names:
 missrows.append([names[db]]+[f'{line(db,v).missing_fraction*100:.1f}%' for v in ['age','weight','bmi','sex','atrial_fibrillation','severity_score']])
 drrows.append([names[db]]+[f'{line(db,"vasopressor_use",v).percent:.1f}%' for v in ['norepinephrine','epinephrine','phenylephrine','dopamine','vasopressin','angiotensin_ii']])
parts[idx:idx]=[table(['数据库','年龄缺失','体重缺失','BMI 缺失','性别缺失','房颤状态缺失','指定评分缺失'],missrows),table(['数据库','去甲肾上腺素','肾上腺素','苯肾上腺素','多巴胺','血管加压素','血管紧张素 II'],drrows),'主要跨库差异：MIMIC 的 BMI 缺失高于 eICU；Amsterdam 的区间开放端影响年龄/体格数值可得性。记录房颤比例 MIMIC 33.6%、eICU 10.0%、Amsterdam 39.3%，心律与诊断覆盖不同，不解释为可直接比较的患病率。MIMIC 苯肾上腺素、Amsterdam 多巴胺的记录使用比例较高，保留数据库及单元构成，不把差异推断为用药效果。药物表的 0 表示选定映射和治疗定义下未观察到使用，不能证明库外不存在用药。']
(O/'phase2_table1_notes.md').write_text('\n\n'.join(parts)+'\n',encoding='utf-8')
pd.DataFrame([dict(check='all_aggregate_invariants',passed=True),dict(check='source_field_discrepancies_zero',passed=au.baseline_field_discrepancies.sum()==0),dict(check='source_peak_discrepancies_zero',passed=au.peak_point_discrepancies.sum()==0)]).to_csv(O/'table1_final_checks.csv',index=False)
print('Table 1 notes and read-only aggregate checks complete.')
