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
import zipfile,shutil,json,time,subprocess
p=Path('@CUFF_WORKSPACE@/phase0/_restricted/amsterdam_numeric_inner.zip')
source=Path('@AMSTERDAM_ROOT@/numericitems.zip')
if not p.exists():
    with zipfile.ZipFile(source) as z:
        members=z.infolist()
        assert len(members)==1 and members[0].filename=='numericitems.zip'
        expected=members[0].file_size
    with p.with_suffix('.zip.part').open('wb') as b:
        run=subprocess.run(['7z','x','-so','-bd',str(source),'numericitems.zip'],stdout=b,stderr=subprocess.PIPE)
    assert run.returncode==0, run.stderr.decode(errors='replace')
    assert p.with_suffix('.zip.part').stat().st_size==expected
    p.with_suffix('.zip.part').replace(p)
with zipfile.ZipFile(p) as z:
    info=[{'name':q.filename,'size':q.file_size,'compressed':q.compress_size} for q in z.infolist()]
Path('@CUFF_WORKSPACE@/phase0/amsterdam_numeric_archive.json').write_text(json.dumps(info,indent=2))
print(json.dumps({'stage':'inner_archive_ready','members':info}))
