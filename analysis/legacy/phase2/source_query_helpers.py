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

"""Minimal query utilities extracted from the independent source audit.

The PostgreSQL client runs locally; database authentication uses PostgreSQL's
standard environment/service configuration and is never embedded here.
"""
from pathlib import Path
import hashlib

def source_path(path):
    return Path(path).as_posix()

def sql_literal(value):
    return "'" + str(value).replace("'", "''") + "'"

def file_digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
