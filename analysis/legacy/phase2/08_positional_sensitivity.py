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

from common import *
from importlib import import_module
path=RES/'audit_suspected_artifact_cids_all_pairs.csv'
if not path.exists():raise RuntimeError('all-pair source Positional evidence not yet available')
flag=pd.read_csv(path)
if 'database' in flag.columns:flag=flag[flag.database.eq('mimic')]
assert 'cid' in flag.columns
assert flag.cid.notna().all()
d=read_pairs('mimic');ids=set(flag.cid.astype(int));eligible=set(d.cid.astype(int))
assert ids.issubset(eligible)
filtered=d[~d.cid.isin(ids)].copy()
ratio_engine=import_module('01_analyses').ratio_engine
r,_=ratio_engine(filtered,'mimic','exclude_source_positional')
save(r[r.domain.isin(['A','B'])],'phase2_sensitivity_exclude_source_positional.csv')
save(pd.DataFrame([dict(database='mimic',pairs_before=len(d),patients_before=d.patient.nunique(),excluded_pairs=len(ids),patients_with_excluded_pairs=d.loc[d.cid.isin(ids),'patient'].nunique(),pairs_after=len(filtered),patients_after=filtered.patient.nunique(),source_label='Arterial line Waveform Appear = Positional',source_time_window_seconds=300,interpretation='source-recorded position-dependent waveform suspicion; not confirmed damping,flush,or erroneous MAP')]),'phase2_positional_sensitivity_counts.csv')
print('SOURCE_POSITIONAL_SENSITIVITY_COMPLETE',flush=True)
