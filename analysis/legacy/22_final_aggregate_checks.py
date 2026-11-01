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

"""Aggregate-only sanity summaries. No agreement strata, models or outcome joins."""
from importlib import import_module
g=import_module('08_core');globals().update({k:v for k,v in vars(g).items() if not k.startswith('_')})
rows=[]
for db in ['mimic','eicu','amsterdam','sicdb']:
    c=connect(db)
    for v in ['weight','height','bmi']:
        q=c.execute(f"SELECT count(*),count({v}),min({v}),quantile_cont({v},.01),quantile_cont({v},.25),median({v}),quantile_cont({v},.75),quantile_cont({v},.99),max({v}) FROM stays s WHERE is_primary AND EXISTS(SELECT 1 FROM drugs d WHERE d.sid=s.sid)").fetchone()
        rows.append(dict(db=db,variable=v,unit={'weight':'kg','height':'cm','bmi':'kg/m2'}[v],scope='treated_primary_patients',**dict(zip(['denominator','known','min','p01','q1','median','q3','p99','max'],q))))
    q=c.execute('SELECT count(*),count(nee),min(nee),quantile_cont(nee,.01),quantile_cont(nee,.25),median(nee),quantile_cont(nee,.75),quantile_cont(nee,.99),max(nee) FROM pairs WHERE is_primary').fetchone()
    rows.append(dict(db=db,variable='nee',unit='mcg/kg/min',scope='primary_pairs_not_patients',**dict(zip(['denominator','known','min','p01','q1','median','q3','p99','max'],q))))
    c.close()
pd.DataFrame(rows).to_csv(OUT/'phase0_variable_distribution_qc.csv',index=False)
print(pd.DataFrame(rows)[['db','variable','median','p99','max']].to_json(orient='records'))
