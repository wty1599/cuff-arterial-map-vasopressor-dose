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
"""Bounded source-supported arterial-quality evidence, never inferred from MAP discordance."""
from pathlib import Path
import json,re,datetime,traceback
import pandas as pd
ROOT=Path('@CUFF_WORKSPACE@');OUT=ROOT/'phase2';R=OUT/'_restricted'
SEARCH=r'damp|flush|waveform|demp|spoel|curve|dämp|spül|artef|kurve|drukgolf'
def dictionary(db):
    if db=='mimic':return pd.read_csv(R/'audit_source_mimic_quality_dictionary.csv',dtype=str).fillna('')
    files={'eicu':['eicu_nursing_dictionary_all.csv','eicu_assessment_dictionary_all.csv','eicu_treatment_dictionary_all.csv'],
           'amsterdam':['amsterdam_dictionary_all.csv']}
    if db=='sicdb':
        d=pd.read_csv('@SICDB_ROOT@/d_references.csv.gz',dtype=str).fillna('')
    else:d=pd.concat([pd.read_csv(ROOT/'phase0'/f,dtype=str).fillna('') for f in files[db]],ignore_index=True).fillna('')
    return d[d.apply(lambda s:s.astype(str).str.contains(SEARCH,case=False,regex=True)).any(axis=1)]
def classify(label,value):
    label=str(label).lower();value=str(value).strip().lower()
    arterial=bool(re.search(r'arterial|\ba[- ]?line\b|\bart line\b|\babp\b',label))
    if not arterial:return 'unassessable','not_an_explicit_arterial_quality_channel'
    # A routine line flush/calibration does not establish a false MAP reading.
    if not re.search(r'waveform|damp|artefact|artifact',label):
        return 'unassessable','line_action_without_waveform_failure_evidence'
    if value=='positional':
        return 'suspected','explicit_position_dependent_arterial_waveform_not_proven_damping'
    if re.fullmatch(r'(?:over|under)?damp(?:ed|ened|ing)?|damp(?:ed|ened)\s+waveform',value):
        return 'suspected','explicit_arterial_damping_record_within_pair_window'
    if re.fullmatch(r'(?:flat|absent|poor|abnormal)(?:\s+waveform)?|artifact|artefact|flush\s+artifact',value):
        return 'suspected','explicit_arterial_waveform_abnormality_within_pair_window'
    if re.fullmatch(r'appropriate|normal|good(?:\s+waveform)?|not damp(?:ed|ened)|undamped',value):
        return 'none_supported','explicit_arterial_waveform_normal_within_pair_window'
    return 'unassessable','quality_value_does_not_establish_damping_or_normality'
def main():
    summaries=[];suspected=[];all_rows=[]
    for db in ['mimic','eicu','amsterdam','sicdb']:
        sample=pd.read_csv(R/f'audit_phase2_{db}.csv',dtype={'sid':str})
        d=dictionary(db);d.to_csv(R/f'audit_artifact_dictionary_{db}.csv',index=False)
        f=R/f'audit_source_{db}_quality.csv'
        if db in ['mimic','eicu'] and not f.exists():raise RuntimeError('Source quality query is not complete')
        q=pd.read_csv(f,dtype=str).fillna('') if f.exists() else pd.DataFrame(columns=['audit_id','source_t','source_label','source_value','source_table'])
        evidence=[]
        for s in sample.itertuples():
            z=q[q.audit_id.eq(s.audit_id)].copy()
            z=z[pd.to_numeric(z.source_t,errors='coerce').sub(s.t).abs().le(300)]
            codes=[classify(x.source_label,x.source_value) for x in z.itertuples()]
            status='suspected' if any(x[0]=='suspected' for x in codes) else 'none_supported' if any(x[0]=='none_supported' for x in codes) else 'unassessable'
            reason=';'.join(sorted({why for what,why in codes})) if codes else 'no_specific_arterial_waveform_quality_record_in_bounded_source'
            row=dict(database=db,audit_id=s.audit_id,cid=s.cid,sid=s.sid,t=s.t,status=status,source_quality_records=len(z),reason=reason)
            evidence.append(row);all_rows.append(row)
            if status=='suspected':suspected.append(row)
        pd.DataFrame(evidence).to_csv(R/f'audit_artifact_evidence_{db}.csv',index=False)
        n=len(evidence);counts=pd.Series([x['status'] for x in evidence]).value_counts()
        row=dict(database=db,audit_pairs=n,quality_records=len(q),dictionary_search_hits=len(d),
            artifact_suspected=int(counts.get('suspected',0)),artifact_none_supported=int(counts.get('none_supported',0)),
            artifact_unassessable=int(counts.get('unassessable',0)),
            artifact_query_scope='fixed-sample +/-5min source quality records; existing full dictionaries; raw BP neighbours',
            completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
        for k in ['artifact_suspected','artifact_none_supported','artifact_unassessable']:row[k+'_fraction']=row[k]/n if n else None
        summaries.append(row)
    pd.DataFrame(all_rows).to_csv(R/'audit_artifact_all_pairs.csv',index=False)
    pd.DataFrame(suspected,columns=['database','audit_id','cid','sid','t','status','source_quality_records','reason']).to_csv(R/'audit_suspected_artifact_cids.csv',index=False)
    pd.DataFrame(summaries).to_csv(OUT/'audit_artifact_summary.csv',index=False,encoding='utf-8-sig')
    print(json.dumps(summaries,ensure_ascii=False))
if __name__=='__main__':
    try:main()
    except Exception:
        (R/'audit_artifact_error.log').write_text(traceback.format_exc(),encoding='utf-8')
        print(json.dumps({'status':'failed','details':'restricted audit_artifact_error.log'}));raise SystemExit(1)
