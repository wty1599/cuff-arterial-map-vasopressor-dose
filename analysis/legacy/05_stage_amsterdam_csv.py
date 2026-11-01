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

from pathlib import Path
import zipfile,shutil,json,time
R=Path('@CUFF_WORKSPACE@/phase0/_restricted');p=R/'amsterdam_numeric_source.csv'
if p.exists():raise SystemExit('Staged CSV already exists; no overwrite.')
t=time.monotonic()
with zipfile.ZipFile(R/'amsterdam_numeric_inner.zip') as z:
    info=z.getinfo('numericitems.csv')
    assert shutil.disk_usage(R).free>info.file_size+30*1024**3
    with z.open(info) as a,p.with_suffix('.csv.part').open('wb') as b:shutil.copyfileobj(a,b,16*1024*1024)
    assert p.with_suffix('.csv.part').stat().st_size==info.file_size
    p.with_suffix('.csv.part').replace(p)
print(json.dumps({'source_csv_ready':True,'bytes':p.stat().st_size,'elapsed_seconds':round(time.monotonic()-t)}))
