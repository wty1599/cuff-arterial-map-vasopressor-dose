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
checks=[]
rates=pd.read_csv(OUT/'phase1_all_rates.csv')
for db in DBS:
    d=pd.read_parquet(RES/f'pairs_{db}.parquet')
    checks.append(dict(check='unique_cuff_anchor',database=db,passed=bool(d.cid.is_unique),details=str(len(d))))
    for ep in ['missed','false_positive']:
        eligible=d.art_map.lt(65) if ep=='missed' else d.art_map.ge(65)
        event=d.cuff_map.ge(65) if ep=='missed' else d.cuff_map.lt(65)
        for lab in ['all','known','unknown',*LABELS]:
            mask=pd.Series(True,index=d.index) if lab=='all' else d.nee.notna() if lab=='known' else d.nee.isna() if lab=='unknown' else d.dose_stratum.eq(lab)
            rr=rates[rates.database.eq(db)&rates.endpoint.eq(ep)&rates.stratum.eq(lab)].iloc[0]
            passed=(int((eligible&event&mask).sum())==rr.numerator and int((eligible&mask).sum())==rr.denominator)
            checks.append(dict(check=f'raw_aggregate_{ep}_{lab}',database=db,passed=bool(passed),details='direct_parquet_recount_not_source_audit'))
        q=rates[rates.database.eq(db)&rates.endpoint.eq(ep)]
        allrow=q[q.stratum.eq('all')].iloc[0];part=q[q.stratum.isin([*LABELS,'unknown'])]
        checks.append(dict(check=f'conservation_{ep}',database=db,passed=bool(part.denominator.sum()==allrow.denominator and part.numerator.sum()==allrow.numerator),details='five strata plus unknown equals total'))
    if db=='mimic':d.assign(delta=d.cuff_map-d.art_map)[['patient','delta']].to_csv(RES/'ba_validation_mimic.csv',index=False)
for path in OUT.glob('phase1_*.csv'):
    z=pd.read_csv(path)
    if {'estimate','denominator','numerator'}.issubset(z.columns) and 'endpoint'in z.columns:
        ok=(z.numerator<=z.denominator).all() and z.estimate.dropna().between(0,1).all()
        checks.append(dict(check='valid_proportions',database=path.name,passed=bool(ok),details='counts and probabilities in range'))
    if {'ci_lower','ci_upper'}.issubset(z.columns):
        good=z[['ci_lower','ci_upper']].notna().all(axis=1)
        checks.append(dict(check='ordered_ci',database=path.name,passed=bool((z.loc[good,'ci_lower']<=z.loc[good,'ci_upper']).all()),details='NA intervals kept as NA'))
cv=pd.read_csv(OUT/'phase1_pair_coverage.csv')
checks.append(dict(check='coverage_not_over_100pct',database='all',passed=bool((cv.numerator<=cv.denominator).all()),details='shared observed span'))
save(pd.DataFrame(checks),'phase1_validation.csv')
if not all(x['passed'] for x in checks): raise RuntimeError('impossible_count_or_consistency_failure')
print(json.dumps({'aggregate_checks':len(checks),'all_passed':True}))
