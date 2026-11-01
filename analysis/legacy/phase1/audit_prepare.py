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
"""Freeze independent phase-1 audit samples; never mutate phase-0 databases."""
from pathlib import Path
import sys, json, csv, hashlib, datetime
ROOT=Path('@CUFF_WORKSPACE@')

import duckdb
OUT=ROOT/'phase1'; RES=OUT/'_restricted'; RES.mkdir(parents=True,exist_ok=True)
SEED=20261027
def main():
    rows=[]
    for db in ['mimic','eicu','amsterdam','sicdb']:
        target=RES/f'audit_phase1_{db}.csv'
        if target.exists():
            raise RuntimeError('Frozen sample already exists; do not resample')
        c=duckdb.connect(str(ROOT/'phase0'/'_restricted'/f'{db}.duckdb'),read_only=True)
        c.execute("SET threads=1; SET memory_limit='700MB'")
        old=(ROOT/'phase0'/'_restricted'/f'audit_input_{db}.csv').as_posix()
        c.execute(f"CREATE TEMP TABLE old_sample AS SELECT CAST(sid AS VARCHAR) sid,CAST(t AS DOUBLE) t FROM read_csv_auto('{old}',all_varchar=true)")
        c.execute("""CREATE TEMP TABLE candidates AS SELECT p.*,cf.episode_start,cf.episode_stop,
          s.weight expected_phase0_weight,s.height expected_phase0_height,
          EXISTS(SELECT 1 FROM old_sample o WHERE o.sid=p.sid AND abs(o.t-p.t)<.0001) old_sample_pair
          FROM pairs p JOIN cuff60 cf USING(cid) JOIN stays s ON p.sid=s.sid
          WHERE p.is_primary AND p.art_map<65""")
        n,new=c.execute('SELECT count(*),count(*) FILTER(WHERE NOT old_sample_pair) FROM candidates').fetchone()
        q=f"SELECT * FROM candidates ORDER BY old_sample_pair,sha256(concat('{SEED}|{db}|',sid,'|',CAST(t AS VARCHAR),'|',CAST(cid AS VARCHAR))) LIMIT 50"
        d=c.execute(q).df()
        assert len(d)==50 and d.cid.nunique()==50
        d.insert(0,'audit_id',[f'{db}_phase1_{i:03d}' for i in range(1,51)])
        d.to_csv(target,index=False,encoding='utf-8-sig')
        rows.append(dict(database=db,eligible_low_map_pairs=n,new_candidates_excluding_old=new,sampled_pairs=len(d),overlap_old_sample_pairs=int(d.old_sample_pair.sum()),seed=SEED,sampling='SHA256 ordered sampling of pair identifiers with fixed seed; prefer pairs absent from prior source audit',created_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),timing='阶段0结果之后；阶段1主分层结果之前固定样本',sample_sha256=hashlib.sha256(target.read_bytes()).hexdigest()))
        c.close()
    with (OUT/'phase1_audit_sample_manifest.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps(rows,ensure_ascii=False))
if __name__=='__main__':main()
