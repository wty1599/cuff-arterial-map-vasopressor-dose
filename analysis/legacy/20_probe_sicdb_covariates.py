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
a=read(SIC/'cases.csv.gz',dtype=str)
rows=[]
for v in ['WeightOnAdmission','HeightOnAdmission']:
    x=num(a[v]);q=x.quantile([0,.01,.25,.5,.75,.99,1])
    rows.append({'field':v,'nonnull':int(x.notna().sum()),'positive':int((x>0).sum()),**{str(k):float(v) for k,v in q.items()}})
pd.DataFrame(rows).to_csv(OUT/'sicdb_covariate_unit_probe.csv',index=False)
af=a[a.ICD10Main.str.contains('^I48',regex=True)|a.ICD10MainText.str.contains('atrial fibr|vorhofflim',case=False,regex=True)]
af.groupby(['ICD10Main','ICD10MainText'],dropna=False).size().reset_index(name='source_stays').to_csv(OUT/'sicdb_af_dictionary_probe.csv',index=False)
print(pd.DataFrame(rows).to_json(orient='records'))
print(af.groupby(['ICD10Main','ICD10MainText'],dropna=False).size().reset_index(name='source_stays').to_json(orient='records'))
