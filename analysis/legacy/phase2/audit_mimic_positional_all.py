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

"""One bounded source SQL query of an existing quality field around all MIMIC primary pairs."""
from pathlib import Path
import os,sys,json,datetime,traceback,subprocess,hashlib
ROOT=Path('@CUFF_WORKSPACE@')
import duckdb,pandas as pd
from source_query_helpers import source_path,sql_literal,file_digest
OUT=ROOT/'phase2';R=OUT/'_restricted';SQL=ROOT/'sql'/'phase2'
def main():
    R.mkdir(parents=True,exist_ok=True);SQL.mkdir(parents=True,exist_ok=True)
    p=ROOT/'phase1'/'_restricted'/'pairs_mimic.parquet';before=file_digest(p)
    c=duckdb.connect();c.execute("SET threads=1;SET memory_limit='700MB'")
    keys=c.execute(f"SELECT cid,sid,t FROM read_parquet('{p.as_posix()}')").df();c.close()
    frozen=pd.read_csv(ROOT/'phase1'/'phase1_all_rates.csv')
    expected=frozen.loc[(frozen.database=='mimic') & (frozen.scenario=='primary_5min') &
                        (frozen.endpoint=='missed') & (frozen.stratum=='all'),'all_pairs']
    if len(expected)!=1:raise RuntimeError('Expected one frozen primary aggregate count')
    if len(keys)!=int(expected.iloc[0]) or keys.cid.nunique()!=len(keys):
        raise RuntimeError('Unexpected protected primary-pair count')
    keyfile=R/'audit_mimic_all_pair_keys.csv';keys.to_csv(keyfile,index=False)
    evidence=R/'audit_mimic_waveform_all_pairs.csv'
    query="""WITH waveform AS MATERIALIZED (
      SELECT c.stay_id,c.charttime,c.value
      FROM mimiciv_icu.chartevents c JOIN mimiciv_icu.d_items d ON c.itemid=d.itemid
      JOIN (SELECT DISTINCT stay_id FROM audit_targets) s USING(stay_id)
      WHERE d.label='Arterial line Waveform Appear'
    )
    SELECT a.cid,a.stay_id AS sid,extract(epoch FROM a.anchor) AS t,
      count(w.charttime) AS quality_records,
      count(w.charttime) FILTER(WHERE w.value='Positional') AS positional_records,
      min(abs(extract(epoch FROM(w.charttime-a.anchor)))) FILTER(WHERE w.value='Positional') AS nearest_positional_diff_seconds
    FROM audit_targets a LEFT JOIN waveform w ON w.stay_id=a.stay_id
      AND w.charttime BETWEEN a.anchor-interval '5 minutes' AND a.anchor+interval '5 minutes'
    GROUP BY a.cid,a.stay_id,a.anchor"""
    text="\\set ON_ERROR_STOP on\nSET client_min_messages=warning;\nSET TIME ZONE 'UTC';\nCREATE TEMP TABLE audit_input(cid text,sid text,t text);\n"
    text+=f"\\copy audit_input FROM {sql_literal(source_path(keyfile))} CSV HEADER\n"
    text+="CREATE TEMP TABLE audit_targets AS SELECT cid::bigint cid,sid::bigint stay_id,timestamp '1970-01-01'+t::double precision*interval '1 second' anchor FROM audit_input;\nCREATE INDEX ON audit_targets(stay_id,anchor);\nBEGIN READ ONLY;\nSET LOCAL statement_timeout=0;\nSET LOCAL work_mem='16MB';\nSET LOCAL max_parallel_workers_per_gather=0;\n"
    text+=f"\\copy ({query.replace(chr(10),' ')}) TO {sql_literal(source_path(evidence))} CSV HEADER\nCOMMIT;\n"
    script=SQL/'audit_mimic_positional_all.sql';script.write_text(text,encoding='utf-8')
    start=datetime.datetime.now(datetime.timezone.utc).isoformat()
    command=[os.environ.get('CUFF_PSQL','psql'),'-X','-v','ON_ERROR_STOP=1','-d',os.environ.get('CUFF_MIMIC_DATABASE','mimiciv'),'-f',source_path(script)]
    result=subprocess.run(command,capture_output=True,text=True,encoding='utf-8',errors='replace')
    (R/'audit_mimic_positional_all.log').write_text(result.stdout+'\n'+result.stderr,encoding='utf-8')
    if result.returncode:raise RuntimeError('Read-only positional source query failed')
    end=datetime.datetime.now(datetime.timezone.utc).isoformat()
    d=pd.read_csv(evidence)
    if len(d)!=len(keys) or d.cid.nunique()!=len(d):raise RuntimeError('All-pair quality coverage count is impossible')
    a=d[d.positional_records.gt(0)].copy();a.insert(0,'database','mimic')
    a['source_label']='Arterial line Waveform Appear';a['source_value']='Positional';a['status']='suspected_position_dependent_waveform'
    a.to_csv(R/'audit_suspected_artifact_cids_all_pairs.csv',index=False)
    row=dict(database='mimic',primary_pairs=len(d),pairs_with_any_waveform_record=int(d.quality_records.gt(0).sum()),
        pairs_with_positional_record=len(a),pairs_without_quality_record=int(d.quality_records.eq(0).sum()),
        waveform_record_coverage=float(d.quality_records.gt(0).mean()),positional_flag_fraction=float(d.positional_records.gt(0).mean()),
        patients_with_positional_flag=int(a.sid.nunique()),source_query_started_at=start,source_query_completed_at=end,
        source_queries_executed=True,source_label='Arterial line Waveform Appear',source_value='Positional',window_minutes_each_side=5,
        timing='After phase 1 and initial phase 2 results; source positional-waveform sensitivity added after the source-label review',
        interpretation='Source positional-waveform flag only; absence of flag does not establish normality; no adjudication of damping/flush or erroneous MAP',
        phase1_input_sha256=before,sql_sha256=file_digest(script))
    assert before==file_digest(p)
    pd.DataFrame([row]).to_csv(OUT/'audit_mimic_positional_coverage.csv',index=False,encoding='utf-8-sig')
    print(json.dumps(row,ensure_ascii=False))
if __name__=='__main__':
    try:main()
    except Exception:
        (R/'audit_mimic_positional_all_error.log').write_text(traceback.format_exc(),encoding='utf-8')
        print(json.dumps({'status':'failed','details':'restricted audit_mimic_positional_all_error.log'}));raise SystemExit(1)
