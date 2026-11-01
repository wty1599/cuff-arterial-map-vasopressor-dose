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
"""Independent audit-only anthropometry cleaning; no production imports."""
from pathlib import Path
import json, math, datetime, statistics
ROOT=Path('@CUFF_WORKSPACE@')
SCOPE=json.loads((Path(__file__).parent/'phase2_parameters.json').read_text(encoding='utf-8'))
assert SCOPE['phase1_definitions_unchanged'] is True
P=json.loads((ROOT/'sql'/'phase1'/'phase1_parameters.json').read_text(encoding='utf-8'))
def number(x):
    try:
        v=float(x)
        return v if math.isfinite(v) else None
    except (TypeError,ValueError):return None
def checked(value,bounds):
    v=number(value)
    return v if v is not None and bounds[0]<=v<=bounds[1] else None
def anthro(weight,height):
    w=checked(weight,P['weight_kg_range']);h=checked(height,P['height_cm_range'])
    b=checked(w/(h/100)**2 if w is not None and h is not None else None,P['bmi_range'])
    return {'source_weight':number(weight),'source_height':number(height),'source_clean_weight':w,'source_clean_height':h,'source_clean_bmi':b}
def mimic_anthro(stays,measurements):
    result={}
    epoch=datetime.datetime(1970,1,1)
    for key,stay in stays.items():
        start=(datetime.datetime.fromisoformat(stay['intime'])-epoch).total_seconds()
        weights=[];heights=[]
        for r in measurements:
            if int(float(r['stay_id']))!=key:continue
            t=number(r.get('source_t'));v=number(r.get('value'));label=r.get('label','')
            if t is None or v is None or v<=0:continue
            if 'Weight' in label and t>=start:
                if 'lbs' in label:v*=.45359237
                weights.append((t,v))
            if 'Height' in label:
                if label=='Height':v*=2.54
                heights.append((abs(t-start),t,label,v))
        weight=None;height=None
        if weights:
            earliest=min(x[0] for x in weights)
            weight=statistics.median(x[1] for x in weights if x[0]==earliest)
        if heights:
            first=min(x[:3] for x in heights)
            height=statistics.median(x[3] for x in heights if x[:3]==first)
        result[key]=anthro(weight,height)
    return result
