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

import os as _cuff_env
"""Synthetic branch checks, explicitly separate from source-data verification."""
from pathlib import Path
import json,csv,datetime
import audit_pg_source as a
import audit_file_source as f
from audit_quality import anthro
ROOT=Path('@CUFF_WORKSPACE@')
p=json.loads((ROOT/'sql/phase0_parameters.json').read_text(encoding='utf-8'))
p['weight_quality_range']=[30,300]
checks=[]
def check(name,passed):
    checks.append(dict(check=name,passed=bool(passed),evidence_type='synthetic_only_not_source_data',at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
check('perkg_mass_rate_does_not_require_weight',a.nee_from_rate('norepinephrine',.5,'mcg/kg/min',1,p)[0]==.5)
check('absolute_mass_rate_invalid_weight_unknown',a.nee_from_rate('norepinephrine',7,'mcg/min',1,p)[0] is None)
check('vasopressin_u_min_does_not_require_weight',a.nee_from_rate('vasopressin',.04,'units/min',1,p)[0]==.1)
check('vasopressin_u_kg_min_invalid_weight_unknown',a.nee_from_rate('vasopressin',.04,'units/kg/min',1,p)[0] is None)
check('file_perkg_mass_rate_kept',f.convert(.5,'mcg','min',True,1,'norepinephrine',quality=True)==.5)
check('file_absolute_mass_invalid_weight_unknown',f.convert(7,'mcg','min',False,1,'norepinephrine',quality=True) is None)
check('anthropometry_lower_bounds_included',anthro(30,120)['source_clean_bmi'] is not None)
check('invalid_weight_invalidates_bmi',anthro(29,180)['source_clean_bmi'] is None)
check('recomputed_bmi_above80_unknown',anthro(300,120)['source_clean_bmi'] is None)
s={1:{'unitdischargeoffset':'100','admissionweight':'10'}}
def row(name,rate):return {'stay_id':'1','infusionoffset':'0','drugname':name,'drugrate':str(rate),'infusionrate':'','drug':'norepinephrine'}
rs,_=a.eicu_intervals(s,[row('Norepinephrine (mcg/kg/min)',.1),row('Norepinephrine (ml/hr)',1)],p)
check('unconvertible_alias_does_not_invalidate_perkg',rs[1][0]['nee']==.1)
rs,_=a.eicu_intervals(s,[row('Norepinephrine (mcg/kg/min)',.1),row('Norepinephrine (mcg/min)',1)],p)
check('invalid_weight_participating_alias_unknown',rs[1][0]['nee'] is None)
with (ROOT/'phase1/audit_synthetic_checks.csv').open('w',encoding='utf-8-sig',newline='') as h:
    w=csv.DictWriter(h,fieldnames=list(checks[0]));w.writeheader();w.writerows(checks)
assert all(x['passed'] for x in checks)
print(json.dumps({'synthetic_checks':len(checks),'passed':sum(x['passed'] for x in checks),'source_queries_executed':False}))
