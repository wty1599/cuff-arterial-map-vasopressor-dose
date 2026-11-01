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
def build(db):
    c=connect(db)
    pids=c.execute('SELECT DISTINCT pid FROM pairs WHERE is_primary ORDER BY pid').fetchnumpy()['pid'];rng=np.random.default_rng(P['seed']);n=math.ceil(len(pids)*.1)
    sampled=np.sort(rng.choice(pids,n,replace=False)) if n else np.array([],dtype=str)
    save(c,'signal_patient_sample',pd.DataFrame({'pid':sampled}));pd.DataFrame({'pid':sampled}).to_csv(R/f'signal_sample_{db}.csv',index=False)
    d=c.execute('SELECT p.pid,count(*) FILTER(WHERE art_map<65) denom,count(*) FILTER(WHERE art_map<65 AND cuff_map>=65) missed FROM pairs p JOIN signal_patient_sample s USING(pid) WHERE is_primary GROUP BY p.pid ORDER BY p.pid').df()
    total=int(d.denom.sum());missed=int(d.missed.sum());rates=[]
    rng=np.random.default_rng(P['seed'])
    if total and len(d):
        den=d.denom.to_numpy();numr=d.missed.to_numpy()
        for i in range(2000):
            ix=rng.integers(0,len(d),len(d));z=den[ix].sum()
            if z:rates.append(float(numr[ix].sum()/z))
    lo,hi=np.quantile(rates,[.025,.975]) if rates else (np.nan,np.nan)
    row=dict(db=db,eligible_paired_patients=len(pids),sampled_patients=n,patients_with_low_arterial_map=int((d.denom>0).sum()),arterial_map_lt65_pairs=total,cuff_ge65_among_these=missed,unstratified_fraction=missed/total if total else np.nan,ci95_lower=lo,ci95_upper=hi,bootstrap_valid_replicates=len(rates),bootstrap_requested_replicates=2000,seed=P['seed'],ci_method='patient_cluster_bootstrap_percentile',scope='unstratified_only')
    pd.DataFrame([row]).to_csv(OUT/f'phase0_signal_{db}.csv',index=False)
    ids=c.execute('SELECT cid FROM pairs WHERE is_primary ORDER BY sid,t').fetchnumpy()['cid'];rng=np.random.default_rng(P['seed']);chosen=rng.choice(ids,min(100,len(ids)),replace=False) if len(ids) else np.array([],dtype='int64');save(c,'audit_ids',pd.DataFrame({'cid':chosen}))
    audit=c.execute('SELECT p.cid audit_id,p.sid,p.t,p.cuff_map,p.art_map,p.nee,p.nee_status,p.nearest_diff,p.nearest_art_t,p.art_values,cf.episode_start,cf.episode_stop FROM pairs p JOIN audit_ids a USING(cid) JOIN cuff60 cf USING(cid) ORDER BY p.cid').df();audit.to_csv(R/f'audit_input_{db}.csv',index=False)
    print(json.dumps({'db':db,'stage':'unstratified_signal_and_audit_sample_ready','sampled_patients':n,'audit_pairs':len(audit)}));c.close()
if __name__=='__main__':run_guard(sys.argv[1],build)
