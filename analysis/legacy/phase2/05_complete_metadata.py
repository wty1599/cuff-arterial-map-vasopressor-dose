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
s=pd.read_csv(OUT/'phase2_C_standardized.csv');sev=pd.read_csv(OUT/'phase2_C_arterial_severity.csv')
for i,r in s.iterrows():s.loc[i,'denominator_patients']=int(sev[sev.database.eq(r.database)&sev.scenario.eq(r.scenario)&sev.stratum.eq(r.stratum)].iloc[0].patients)
save(s,'phase2_C_standardized.csv')
d=read_pairs('eicu');r=pd.read_csv(INP/'eicu_recovered_pairs.csv',usecols=['cid','recovered_nee']);d=d.merge(r,on='cid',validate='one_to_one')
restore=d.nee_raw.isna()&d.recovered_nee.notna()&~d.dose_weight_flag&d.recovered_nee.between(0,5)
d.loc[restore,'nee']=d.loc[restore,'recovered_nee'];d['dose_stratum']=dose_strata(d.nee)
rows=[]
for lab in ['all','>1.0']:
 z=d[(np.ones(len(d),bool) if lab=='all' else d.dose_stratum.eq(lab))&d.art_map.lt(65)&d.cuff_map.ge(65)];delta=z.cuff_map-z.art_map
 rows.append(dict(database='eicu',scenario='recovered_dose',stratum=lab,missed_pairs=len(z),missed_patients=z.patient.nunique(),difference_median=delta.median(),difference_q1=delta.quantile(.25),difference_q3=delta.quantile(.75)))
spec=pd.read_csv(RES/'ratio_specs_eicu_recovered_dose.csv');dist=spec[spec.domain.eq('A_distribution')].merge(pd.DataFrame(rows),on=['database','scenario','stratum'],validate='many_to_one')
save(dist,'phase2_sensitivity_recovered_miss_distribution.csv')
save(pd.DataFrame([dict(item='EVERDAC_norepinephrine_tartrate_threshold',value=2.5,unit='ug/kg/min',source='https://doi.org/10.1056/NEJMoa2502136',public_full_text='https://internetbookofemergencymedicine.com/wp-content/uploads/2026/04/art-line.pdf',boundary='Protocol sums norepinephrine tartrate plus epinephrine'),dict(item='approximate_pure_norepinephrine_base_equivalent',value=1.25,unit='ug/kg/min',source='https://doi.org/10.1056/NEJMoa2502136',public_full_text='https://internetbookofemergencymedicine.com/wp-content/uploads/2026/04/art-line.pdf',boundary='Pure norepinephrine conversion only; not generic multivasopressor NE-eq')]),'phase2_external_reference.csv')
print('METADATA_COMPLETE_NO_RATE_CHANGES')
