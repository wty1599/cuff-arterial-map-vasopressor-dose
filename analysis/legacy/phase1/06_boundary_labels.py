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
rows=[]
for path in sorted(OUT.glob('phase1_*.csv')):
    d=pd.read_csv(path)
    if not {'boundary_events','ci_lower','ci_upper','precision','denominator'}.issubset(d.columns):continue
    flag=d.boundary_events.eq(True)&d.ci_lower.eq(d.ci_upper)&d.denominator.gt(0)
    d['ci_status']=np.where(flag,'boundary_degenerate',np.where(d.denominator.gt(0),'ok','no_denominator'))
    d.loc[flag,'precision']='boundary_not_informative'
    save(d,path.name)
    if flag.any():
        z=d[flag][['database','scenario','endpoint','stratum','numerator','denominator','estimate','ci_lower','ci_upper']].copy();z['file']=path.name;rows.append(z)
save(pd.concat(rows,ignore_index=True),'phase1_boundary_ci_flags.csv')
print('BOUNDARY_LABELS_COMPLETE')
