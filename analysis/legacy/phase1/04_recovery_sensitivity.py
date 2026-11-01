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
summarize=import_module('02_rates').summarize
d=pd.read_parquet(RES/'pairs_eicu.parquet')
r=pd.read_csv(RES/'eicu_recovered_pairs.csv',usecols=['cid','recovered_nee'])
d=d.merge(r,on='cid',validate='one_to_one')
restore=d.nee_raw.isna()&d.recovered_nee.notna()&~d.dose_weight_flag&d.recovered_nee.between(0,5)
d.loc[restore,'nee']=d.loc[restore,'recovered_nee']
d['dose_stratum']=strata(d.nee)
result,_=summarize(d,'eicu','eicu_recovered_dose')
save(result,'phase1_sensitivity_eicu_recovered_dose.csv')
save(pd.DataFrame([dict(pairs=len(d),restored_pairs=int(restore.sum()),restored_patients=d.loc[restore,'patient'].nunique(),known_after=int(d.nee.notna().sum()),known_fraction=d.nee.notna().mean())]),'phase1_eicu_recovery_applied.csv')
print('RECOVERY_SENSITIVITY_COMPLETE',flush=True)
