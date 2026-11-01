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
c=connect('mimic');cq=csvq(R/'mimic_bp_cov.csv');oq=csvq(R/'mimic_omr.csv')
exported(c,'mimic_height_source_probe.csv',f"SELECT m.label,count(*) n_records,count(DISTINCT x.stay_id) n_stays,count(*) FILTER(WHERE epoch(try_cast(x.charttime AS TIMESTAMP))<s.t0) before_icu_records,count(*) FILTER(WHERE epoch(try_cast(x.charttime AS TIMESTAMP))>=s.t0) after_icu_records FROM {cq} x JOIN cov_map m ON x.itemid=m.itemid JOIN stays s ON x.stay_id=s.sid WHERE m.label ILIKE '%height%' GROUP BY m.label")
exported(c,'mimic_omr_height_names.csv',f"SELECT result_name,count(*) n_records,count(DISTINCT subject_id) n_patients FROM {oq} WHERE result_name ILIKE '%height%' GROUP BY result_name")
print('Height source availability summaries prepared; no cohort definition changed.')
