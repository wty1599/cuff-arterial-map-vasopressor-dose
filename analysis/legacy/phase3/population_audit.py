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
from pathlib import Path
import sys,json,re
ROOT=Path('@CUFF_WORKSPACE@');OUT=ROOT/'phase3';RES=OUT/'_restricted'
import pandas as pd,numpy as np,duckdb
rules=json.loads((OUT/'population_parameters.json').read_text(encoding='utf-8'))
rows=[];details=[]
for db in ['mimic','eicu']:
    expected=pd.read_csv(RES/f'population_audit_expected_{db}.csv',dtype=str,keep_default_na=False).set_index('sid')
    raw=pd.read_csv(RES/f'population_audit_raw_{db}.csv',dtype=str,keep_default_na=False)
    for sid,z in raw.groupby('sid',sort=True):
        e=expected.loc[sid];unit=z.iloc[0]['first_careunit' if db=='mimic' else 'unittype'].strip()
        if db=='mimic':
            s1=unit not in ['Cardiac Vascular Intensive Care Unit (CVICU)'] if unit else False
            s2=any(v.lower() in ('t','true','1') for v in z.sepsis3)
            stat='confirmed_sepsis3_positive' if s2 else ('not_listed_in_positive_sepsis3_table' if all(not v for v in z.sepsis_stay_id) else 'recorded_false_in_existing_sepsis3_table' if any(v.lower() in ('f','false','0') for v in z.sepsis3) else 'null_only_in_existing_sepsis3_table')
        else:
            s1=unit not in ['CSICU','CTICU','CCU-CTICU'] if unit else False
            dx=z.iloc[0].apacheadmissiondx.strip();unknown=dx.lower() in ('','unknown','none','null','nan','other/unknown')
            s2=(re.search('sepsis',dx,re.IGNORECASE) is not None) and not unknown
            stat='recorded_sepsis_APACHE_diagnosis' if s2 else 'unknown_APACHE_admission_diagnosis' if unknown else 'recorded_other_APACHE_diagnosis'
        comparison={'raw_index_unit_matches':unit==e['first_careunit' if db=='mimic' else 'unittype'].strip(),'S1_retention_matches':s1==(e.S1=='True'),'S2_retention_matches':s2==(e.S2=='True'),'S2_ascertainment_state_matches':stat==e.s2_status}
        if db=='eicu':comparison['raw_APACHE_diagnosis_matches']=dx==e.apacheadmissiondx.strip()
        else:comparison['source_positive_row_count_matches']=sum(v.lower() in ('t','true','1') for v in z.sepsis3)==int(e.sepsis_positive_rows)
        for name,passed in comparison.items():details.append(dict(database=db,sid=sid,check=name,passed=passed))
    dd=pd.DataFrame(details);dd=dd[dd.database.eq(db)]
    rows.append(dict(database=db,sampled_patients=len(expected),source_patients=raw.sid.nunique(),field_checks=len(dd),mismatched_checks=int((~dd.passed).sum()),mismatched_patients=dd[~dd.passed].sid.nunique(),sample_scheme='enriched_joint_S1_S2_states_seed20261027',validation_boundary='Source careunit and existing diagnostic label; no new clinical sepsis adjudication'))
pd.DataFrame(details).to_csv(RES/'population_source_audit_details.csv',index=False)
pd.DataFrame(rows).to_csv(OUT/'population_source_audit.csv',index=False,encoding='utf-8-sig')
assert pd.DataFrame(details).passed.all()

# Independent SQL counts from clean input pairs joined with independently source-backed labels.
c=duckdb.connect();c.execute("SET threads=1;SET memory_limit='700MB'");qc=[]
res=pd.read_csv(OUT/'phase3_S1_S2_results.csv')
for db in ['mimic','eicu']:
    pairpath=(ROOT/f'phase1/_restricted/pairs_{db}.parquet').as_posix();labelpath=(RES/f'population_labels_{db}.csv').as_posix()
    for scenario in ['S1','S2']:
        for r in res[res.database.eq(db)&res.sensitivity.eq(scenario)].itertuples():
            dose='true' if r.stratum=='all' else f"p.dose_stratum='{r.stratum}'"
            condition='p.art_map<65' if r.metric=='missed_hypotension' else 'p.art_map<55' if r.metric=='deep_missed_hypotension' else 'p.cuff_map>=65'
            event='p.art_map<65' if r.metric=='art_lt65_given_cuff_ge65' else 'p.cuff_map>=65'
            q=f"""SELECT count(*) AS layer_pairs,count(*) FILTER(WHERE p.art_map<65) AS layer_low,
            count(*) FILTER(WHERE {condition}) AS denominator,
            count(*) FILTER(WHERE ({condition}) AND ({event})) AS numerator,
            count(DISTINCT p.patient) FILTER(WHERE {condition}) AS denominator_patients
            FROM read_parquet('{pairpath}') p JOIN read_csv_auto('{labelpath}') k ON cast(p.sid AS VARCHAR)=cast(k.sid AS VARCHAR)
            WHERE k.{scenario}=true AND ({dose})"""
            v=c.execute(q).fetchone();wanted=(r.layer_pairs,r.layer_low_arterial_pairs,r.denominator,r.numerator,r.denominator_patients)
            qc.append(dict(database=db,sensitivity=scenario,metric=r.metric,stratum=r.stratum,checks=5,passed=tuple(v)==tuple(wanted)))
pd.DataFrame(qc).to_csv(OUT/'population_independent_count_checks.csv',index=False,encoding='utf-8-sig')
assert pd.DataFrame(qc).passed.all()
print(json.dumps({'source_patients':sum(r['sampled_patients'] for r in rows),'source_field_checks':sum(r['field_checks'] for r in rows),'source_discrepancies':0,'independent_result_rows':len(qc),'independent_count_checks':sum(r['checks'] for r in qc),'status':'PASS'}))
