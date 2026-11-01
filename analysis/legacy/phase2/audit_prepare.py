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
"""Freeze phase-2 deep-missed-pair samples using read-only phase-1 inputs."""
from pathlib import Path
import sys, json, csv, hashlib, datetime, traceback
ROOT=Path('@CUFF_WORKSPACE@')

import duckdb
OUT=ROOT/'phase2'; RES=OUT/'_restricted'; RES.mkdir(parents=True,exist_ok=True)
P=json.loads((ROOT/'sql'/'phase2'/'phase2_parameters.json').read_text(encoding='utf-8'))
SEED=P['seed']
def main():
    rows=[]
    for db in ['mimic','eicu','amsterdam','sicdb']:
        target=RES/f'audit_phase2_{db}.csv'
        if target.exists():raise RuntimeError('Frozen sample exists; never silently resample')
        source=ROOT/'phase1'/'_restricted'/f'pairs_{db}.parquet'
        c=duckdb.connect(str(ROOT/'phase0'/'_restricted'/f'{db}.duckdb'),read_only=True)
        c.execute("SET threads=1; SET memory_limit='700MB'")
        old0=(ROOT/'phase0'/'_restricted'/f'audit_input_{db}.csv').as_posix()
        old1=(ROOT/'phase1'/'_restricted'/f'audit_phase1_{db}.csv').as_posix()
        c.execute(f"CREATE TEMP TABLE old_sample AS SELECT CAST(sid AS VARCHAR) sid,CAST(t AS DOUBLE) t FROM read_csv_auto('{old0}',all_varchar=true) UNION SELECT CAST(sid AS VARCHAR),CAST(t AS DOUBLE) FROM read_csv_auto('{old1}',all_varchar=true)")
        n0=c.execute(f"SELECT count(*) FROM read_parquet('{source.as_posix()}') WHERE art_map<55 AND cuff_map>=65").fetchone()[0]
        c.execute(f"""CREATE TEMP TABLE candidates AS SELECT p.*,cf.episode_start,cf.episode_stop,
          s.weight expected_phase0_weight,s.height expected_phase0_height,
          EXISTS(SELECT 1 FROM old_sample o WHERE o.sid=p.sid AND abs(o.t-p.t)<.0001) old_sample_pair
          FROM pairs p JOIN cuff60 cf USING(cid) JOIN stays s ON p.sid=s.sid
          JOIN read_parquet('{source.as_posix()}') clean ON p.cid=clean.cid
          WHERE p.is_primary AND clean.art_map<55 AND clean.cuff_map>=65
          AND clean.sid=p.sid AND abs(clean.t-p.t)<.0001""")
        n,new=c.execute('SELECT count(*),count(*) FILTER(WHERE NOT old_sample_pair) FROM candidates').fetchone()
        if n!=n0:raise RuntimeError('Sampling-frame identity count is inconsistent')
        q=f"SELECT * FROM candidates ORDER BY old_sample_pair,sha256(concat('{SEED}|{db}|phase2|',sid,'|',CAST(t AS VARCHAR),'|',CAST(cid AS VARCHAR))) LIMIT 30"
        d=c.execute(q).df();assert len(d)==min(n,30) and d.cid.nunique()==len(d)
        assert ((d.art_map<55)&(d.cuff_map>=65)).all()
        d.insert(0,'audit_id',[f'{db}_phase2_{i:03d}' for i in range(1,len(d)+1)])
        d.to_csv(target,index=False,encoding='utf-8-sig')
        row=dict(database=db,eligible_deep_missed_pairs=n,new_candidates_excluding_old=new,
            sampled_pairs=len(d),shortfall=max(0,30-len(d)),overlap_old_sample_pairs=int(d.old_sample_pair.sum()),
            seed=SEED,sampling='fixed-seed SHA256 order, prior phase0/1 audited pairs last',
            created_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
            timing='已见阶段1结果；阶段2执行时序另列；样本由既有阶段1配对固定',
            sample_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
            phase1_input_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
        rows.append(row);c.close();print(json.dumps(row,ensure_ascii=False),flush=True)
    with (OUT/'audit_sample_manifest.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
if __name__=='__main__':
    try:main()
    except Exception:
        (RES/'audit_prepare_error.log').write_text(traceback.format_exc(),encoding='utf-8')
        print(json.dumps({'status':'failed','details':'restricted audit_prepare_error.log'}));raise SystemExit(1)
