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
db=sys.argv[1];c=connect(db)
exported(c,f'phase0_rounding_{db}.csv',"""SELECT count(*) unique_arterial_times_used,avg((abs(a.map_value/5-round(a.map_value/5))<0.0000001)::INT) ending_0_or_5_fraction FROM arterial a WHERE EXISTS(SELECT 1 FROM cuff60 cf JOIN stays s USING(sid) WHERE s.is_primary AND cf.sid=a.sid AND cf.episode_start=a.episode_start AND abs(a.t-cf.t)<=300)""")
exported(c,f'phase0_pair_integrity_{db}.csv',"""SELECT count(*) pair_rows,count(DISTINCT cid) unique_cuff_ids,sum((art_map<20 OR art_map>200 OR cuff_map<20 OR cuff_map>200)::INT) out_of_range_pairs,sum((nearest_diff>300)::INT) beyond_window,sum((nee_status='known' AND nee IS NULL)::INT) inconsistent_known_status FROM pairs""")
c.close();print(json.dumps({'db':db,'recording_qc_completed':True}))
