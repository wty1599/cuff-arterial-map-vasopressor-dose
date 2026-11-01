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
for db in (sys.argv[1:] or ['mimic','eicu']):
    c=connect(db);parts=[];rng=np.random.default_rng(P['seed'])
    for modality in ['arterial','cuff']:
        missing=c.execute(f"SELECT s.sid FROM treated_primary s WHERE NOT EXISTS(SELECT 1 FROM bp_quality b JOIN episodes60 e ON b.sid=e.sid AND b.t>=e.start AND b.t<e.stop WHERE b.sid=s.sid AND b.modality='{modality}' AND b.map_value IS NOT NULL) ORDER BY s.sid").df()
        chosen=rng.choice(missing.sid.to_numpy(),min(50,len(missing)),replace=False);parts.append(pd.DataFrame({'sid':chosen,'target':modality}))
    sampled=pd.concat(parts,ignore_index=True);save(c,'missing_check_sample',sampled)
    rows=c.execute('SELECT s.sid,s.target,e.start episode_start,e.stop episode_stop FROM missing_check_sample s JOIN episodes60 e USING(sid) ORDER BY s.target,s.sid,e.start').df();rows.to_csv(R/f'missing_check_input_{db}.csv',index=False)
    pd.DataFrame([{'db':db,'sample_patient_target_combinations':len(sampled),'source_intervals':len(rows),'seed':P['seed']}]).to_csv(OUT/f'missing_check_plan_{db}.csv',index=False)
    c.close()
print('Missing-vital verification samples prepared; no rows printed.')
