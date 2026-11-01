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
"""Independent source audit. Does not import production conversion or parser."""
from pathlib import Path
import pandas as pd, numpy as np, re, json
BASE=Path('@CUFF_WORKSPACE@/phase1');R=BASE/'_restricted'
def n(x):
    try:return float(x)
    except (TypeError,ValueError):return np.nan
def canonical(x):
    t=x.upper()
    if 'NOREPI' in t or 'NORADREN' in t or 'LEVOPHED' in t:return 'norepinephrine'
    if 'EPINEPH' in t or 'ADRENALIN' in t:return 'epinephrine'
    if 'PHENYLEPH' in t or re.search('NEO.?SYNEPH',t):return 'phenylephrine'
    if 'DOPAMIN' in t:return 'dopamine'
    if 'VASOPRESS' in t or 'PITRESSIN' in t:return 'vasopressin'
    return ''
def audit_conc(label,kind):
    s=label.upper().replace(',','');opts=[]
    for a in re.finditer(r'(\d+(?:\.\d+)?)\s*(MCG|MG|G|UNITS?|U)\b',s):
        after=s[a.end():]
        m=re.search(r'(\d+(?:\.\d+)?)\s*ML\b',after)
        if m is None or re.search(r'\d',after[:m.start()]):continue
        bridge=after[:m.start()].strip()
        if len(bridge)>50:continue
        amount=n(a[1]);volume=n(m[1]);u=a[2]
        if amount<=0 or volume<50:continue
        if kind=='vasopressin' and u.startswith('U'):opts.append(amount/volume)
        elif kind!='vasopressin' and u in ['MG','MCG','G']:opts.append(amount*{'MG':1000,'MCG':1,'G':1e6}[u]/volume)
    return opts[0] if len(opts)==1 else np.nan
ex=pd.read_csv(R/'recovery_audit_expected.csv',dtype={'infusiondrugid':str,'sid':str})
src=pd.read_csv(R/'recovery_audit_infusion_source.csv',dtype=str,keep_default_na=False).set_index('infusiondrugid')
med=pd.read_csv(R/'recovery_audit_medication_source.csv',dtype=str,keep_default_na=False)
med['drug']=med.drugname.map(canonical)
answer=[]
for e in ex.itertuples(index=False):
    r=src.loc[e.infusiondrugid];d=canonical(r.drugname)
    aw=n(r.admissionweight);iw=n(r.patientweight);w=aw if 30<=aw<=300 else iw if 30<=iw<=300 else np.nan
    tm=n(r.infusionoffset);rate=n(r.drugrate);pump=n(r.infusionrate)
    par=[v for v in re.findall(r'\(([^()]*)\)',r.drugname) if '/' in v]
    u=(par[-1] if par else '').lower().replace(' ','').replace('units','u').replace('unit','u').replace('hour','hr').replace('nanograms','ng')
    result=np.nan;reason='unresolved';concs=[]
    if u=='ml/hr':
        if not (rate>0 and pump>0 and abs(rate-pump)>max(1e-6,rate*.01)):
            flow=rate if rate>0 else pump
            mm=med[(med.patientunitstayid==r.patientunitstayid)&(med.drug==d)]
            for m in mm.itertuples(index=False):
                if m.drugordercancelled.lower()=='yes' or not re.search('IV|INTRAVEN|CENTRAL|PERIPH',m.routeadmin.upper()):continue
                if not (n(m.drugstartoffset)<=tm<n(m.drugstopoffset)):continue
                v=audit_conc(m.drugname,d)
                if np.isfinite(v):concs.append(v)
            if concs and max(concs)-min(concs)<=max(1e-8,max(concs)*.01) and flow>0:
                co=float(np.median(concs))
                if d=='vasopressin':result=flow*co/60*2.5
                elif 30<=w<=300:result=flow*co/60/w*{'norepinephrine':1,'epinephrine':1,'phenylephrine':.1,'dopamine':.01}[d]
    elif rate>0 and (u.endswith('/min') or u.endswith('/hr')):
        time_factor=1/60 if u.endswith('/hr') else 1;kg='/kg/' in u
        if d=='vasopressin' and u.startswith('u/') and (not kg or 30<=w<=300):result=rate*time_factor*(w if kg else 1)*2.5
        elif d!='vasopressin':
            mass={'mcg':1,'mg':1000,'ng':.001,'g':1e6}.get(u.split('/')[0])
            if mass is not None and (kg or 30<=w<=300):result=rate*mass*time_factor/(1 if kg else w)*{'norepinephrine':1,'epinephrine':1,'phenylephrine':.1,'dopamine':.01}[d]
    expected=n(e.candidate_dose)
    equal=(np.isnan(expected) and np.isnan(result)) or (np.isfinite(expected) and np.isfinite(result) and np.isclose(expected,result,rtol=1e-8,atol=1e-10))
    same_source=(str(r.patientunitstayid)==str(e.sid) and n(r.infusionoffset)==n(e.infusionoffset) and r.drugname==e.drugname)
    answer.append(dict(infusiondrugid=e.infusiondrugid,expected=expected,independent=result,dose_matches=equal,source_matches=same_source,concentration_sources=len(concs)))
a=pd.DataFrame(answer);a.to_csv(R/'recovery_audit_comparison.csv',index=False)
s=pd.DataFrame([{'metric':'sampled_source_records','value':len(a)},{'metric':'sampled_recovered_records','value':int(a.expected.notna().sum())},{'metric':'sampled_unresolved_records','value':int(a.expected.isna().sum())},{'metric':'source_record_discrepancies','value':int((~a.source_matches).sum())},{'metric':'dose_discrepancies','value':int((~a.dose_matches).sum())}])
s.to_csv(BASE/'phase1_eicu_recovery_audit.csv',index=False)
print(json.dumps(dict(zip(s.metric,s.value))))
