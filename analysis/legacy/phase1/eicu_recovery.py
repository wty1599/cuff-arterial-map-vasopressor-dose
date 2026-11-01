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
import re, math, datetime

def num(v):
    try:
        x=float(v)
        return x if np.isfinite(x) else np.nan
    except (TypeError,ValueError): return np.nan

def drug(s):
    s=str(s).lower()
    if re.search(r'norepi|noradren|levophed',s): return 'norepinephrine'
    if re.search(r'epineph|adrenalin',s): return 'epinephrine'
    if re.search(r'phenyleph|neo.?syneph|neosyneph',s): return 'phenylephrine'
    if 'dopamin' in s: return 'dopamine'
    if re.search(r'vasopress|pitressin',s): return 'vasopressin'
    return None

def unit(name):
    vals=re.findall(r'\(([^()]*)\)',str(name))
    vals=[v for v in vals if re.search(r'/min|/hr|/hour',v,re.I)]
    return vals[-1].lower().replace(' ','').replace('units','u').replace('unit','u').replace('hour','hr').replace('nanograms','ng') if vals else ''

def concentration(name,drugname):
    # Match a named drug quantity followed by its explicit diluent amount.
    s=str(name).lower().replace(',','')
    matches=re.findall(r'(\d+(?:\.\d+)?)\s*(mcg|mg|g|units?|u)\b\s*(?:in\s*|/\s*|[^\d]{0,50}?)(\d+(?:\.\d+)?)\s*ml\b',s)
    out=[]
    for amount,u,vol in matches:
        amount=float(amount);vol=float(vol)
        if amount<=0 or vol<50: continue
        if drugname=='vasopressin':
            if not u.startswith('u'): continue
            out.append(amount/vol)
        else:
            if u not in ['mcg','mg','g']: continue
            out.append(amount*{'mcg':1,'mg':1000,'g':1e6}[u]/vol)
    if len(out)!=1: return np.nan
    return out[0]

def convert_direct(rate,u,w,d):
    if not np.isfinite(rate) or rate<=0: return np.nan
    if not (u.endswith('/min') or u.endswith('/hr')): return np.nan
    hour=60 if u.endswith('/hr') else 1
    kg='/kg/' in u
    if d=='vasopressin':
        if not u.startswith('u/'): return np.nan
        if kg and not (30<=w<=300): return np.nan
        return rate*(w if kg else 1)/hour*2.5
    fac={'mcg':1,'ng':.001,'mg':1000,'g':1e6}.get(u.split('/')[0])
    if fac is None or (not kg and not 30<=w<=300): return np.nan
    return rate*fac/hour/(1 if kg else w)*{'norepinephrine':1,'epinephrine':1,'phenylephrine':.1,'dopamine':.01}[d]

def run():
    c=con('eicu')
    c.execute(f"CREATE TEMP VIEW raw AS SELECT * FROM read_csv_auto('{(SRC/'eicu_drugs.csv').as_posix()}',all_varchar=true)")
    c.execute(f"CREATE TEMP VIEW drugmap AS SELECT * FROM read_csv_auto('{(ROOT/'phase0/map_vaso_eicu.csv').as_posix()}',all_varchar=true) WHERE selected='1'")
    c.execute("""CREATE TEMP TABLE eligible AS SELECT row_number() over(order by d.sid,d.drug,d.start,d.rawid) gkey,d.* FROM drugs d WHERE d.dose IS NULL AND EXISTS (SELECT 1 FROM pairs p WHERE p.is_primary AND p.nee IS NULL AND p.sid=d.sid AND p.t>=d.start AND p.t<d.end)""")
    rows=c.execute("""SELECT e.gkey,e.sid,e.drug,e.start,e.end,e.rawid,e.dose_reason,s.weight AS admission_weight,r.* FROM eligible e JOIN raw r0 ON e.rawid=r0.infusiondrugid JOIN raw r ON r.patientunitstayid=r0.patientunitstayid AND r.infusionoffset=r0.infusionoffset JOIN drugmap m ON m.itemid=r.drugname AND m.drug=e.drug JOIN stays s ON s.sid=e.sid""").df()
    med=pd.read_csv(RES/'recovery_medication_source.csv',dtype=str,keep_default_na=False)
    med['drug']=med.drugname.map(drug)
    med['concentration']=[concentration(n,d) for n,d in zip(med.drugname,med.drug)]
    med['a']=pd.to_numeric(med.drugstartoffset,errors='coerce')
    med['b']=pd.to_numeric(med.drugstopoffset,errors='coerce')
    med['route_ok']=med.routeadmin.str.contains(r'iv|intraven|central|periph',case=False,regex=True)
    med['valid']=med.concentration.notna()&med.route_ok&~med.drugordercancelled.str.lower().eq('yes')&(med.b>med.a)
    m=med[med.valid].copy()
    bymed={k:g for k,g in m.groupby(['patientunitstayid','drug'])}
    candidates=[]
    resume='--resume' in sys.argv
    for r in ([] if resume else rows.itertuples(index=False)):
        u=unit(r.drugname);rate=num(r.drugrate);ir=num(r.infusionrate);iw=num(r.patientweight);aw=num(r.admission_weight)
        weight=aw if 30<=aw<=300 else iw if 30<=iw<=300 else np.nan
        wsource='admission' if 30<=aw<=300 else 'same_infusion' if 30<=iw<=300 else 'unavailable'
        value=np.nan;method='unresolved';conc=np.nan;mid='';reason='rate_unit_or_concentration_unresolved'
        if r.dose_reason=='conflicting_simultaneous_rates':
            reason='existing_simultaneous_conflict'
        elif u!='ml/hr':
            value=convert_direct(rate,u,weight,r.drug)
            if np.isfinite(value): method='same_infusion_weight';reason='recovered_direct_rate_with_weight'
        else:
            flow=rate if rate>0 else ir
            if rate>0 and ir>0 and abs(rate-ir)>max(1e-6,abs(rate)*.01):
                reason='conflicting_pump_rates'
            elif not np.isfinite(flow) or flow<=0:
                reason='no_positive_pump_rate'
            else:
                ms=bymed.get((str(r.sid),r.drug))
                if ms is None: reason='no_explicit_diluted_order'
                else:
                    tm=num(r.infusionoffset);ms=ms[(ms.a<=tm)&(ms.b>tm)]
                    if len(ms)==0: reason='no_time_covering_diluted_order'
                    elif ms.concentration.max()-ms.concentration.min()>max(1e-8,ms.concentration.max()*.01): reason='conflicting_order_concentrations'
                    else:
                        conc=float(ms.concentration.median());mid='|'.join(sorted(set(ms.medicationid)))
                        if r.drug=='vasopressin': value=flow*conc/60*2.5
                        elif 30<=weight<=300: value=flow*conc/60/weight*{'norepinephrine':1,'epinephrine':1,'phenylephrine':.1,'dopamine':.01}[r.drug]
                        if np.isfinite(value):method='time_matched_explicit_diluted_medication';reason='recovered_from_concentration'
                        else: reason='weight_unknown_or_out_of_range'
        candidates.append(dict(gkey=r.gkey,sid=r.sid,drug=r.drug,start=r.start,end=r.end,rawid=r.rawid,infusiondrugid=r.infusiondrugid,infusionoffset=r.infusionoffset,drugname=r.drugname,unit=u,drugrate=rate,infusionrate=ir,weight=weight,weight_source=wsource,concentration=conc,medicationids=mid,candidate_dose=value,method=method,reason=reason))
    cand=pd.read_csv(RES/'recovery_source_candidates.csv',dtype={'sid':str,'rawid':str,'infusiondrugid':str,'medicationids':str}) if resume else pd.DataFrame(candidates)
    cand['medicationids']=cand.medicationids.fillna('').str.replace(r'^([0-9]+)\.0$',r'\1',regex=True)
    cand.to_csv(RES/'recovery_source_candidates.csv',index=False)
    segments=[]
    for gkey,g in ([] if resume else cand.groupby('gkey',sort=False)):
        v=g.candidate_dose.dropna();conflict=len(v)>1 and v.max()-v.min()>max(1e-6,v.max()*.01)
        dose=float(v.median()) if len(v)>0 and not conflict else np.nan
        first=g.iloc[0]
        segments.append(dict(sid=str(first.sid),drug=first.drug,start=first.start,end=first.end,rawid=str(first.rawid),new_dose=dose,reason='recovered' if np.isfinite(dose) else 'recovery_simultaneous_conflict' if conflict else '|'.join(sorted(set(g.reason))),source_count=len(g),source_record_ids='|'.join(g.infusiondrugid.astype(str)),methods='|'.join(sorted(set(g.loc[g.candidate_dose.notna(),'method'])))))
    seg=pd.read_csv(RES/'recovery_segments.csv',dtype={'sid':str,'rawid':str}) if resume else pd.DataFrame(segments)
    seg.to_csv(RES/'recovery_segments.csv',index=False)
    for col in seg:
        if isinstance(seg[col].dtype,pd.StringDtype):seg[col]=seg[col].astype(object).where(seg[col].notna(),None)
    c.register('recovered_segments_frame',seg)
    c.execute("""CREATE TEMP VIEW augmented AS SELECT d.*,coalesce(r.new_dose,d.dose) recovery_dose,CASE WHEN r.new_dose IS NOT NULL THEN r.reason ELSE d.dose_reason END recovery_reason FROM drugs d LEFT JOIN recovered_segments_frame r ON d.sid=r.sid AND d.drug=r.drug AND d.start=r.start AND d.end=r.end AND d.rawid=r.rawid""")
    p=c.execute("""WITH actual AS (SELECT p.cid,count(DISTINCT d.drug) n_active,sum(d.recovery_dose) total,count(*) FILTER(WHERE d.drug IS NOT NULL AND d.recovery_dose IS NULL) n_unknown,count(*) FILTER(WHERE d.recovery_reason='recovered') n_recovered FROM pairs p LEFT JOIN augmented d ON p.sid=d.sid AND p.t>=d.start AND p.t<d.end WHERE p.is_primary GROUP BY p.cid) SELECT p.cid,p.sid,p.pid,p.t,p.nee original_nee,p.nee_status original_status,e.n_expected,a.*,CASE WHEN p.nee IS NULL AND e.n_expected>0 AND e.n_expected=a.n_active AND a.n_unknown=0 THEN a.total END recovery_candidate FROM pairs p JOIN pair_expected e ON e.cid=p.cid JOIN actual a ON a.cid=p.cid WHERE p.is_primary""").df()
    if 'cid_1' in p:p=p.drop(columns='cid_1')
    p['recovered_nee']=p.recovery_candidate.where(p.recovery_candidate.le(P['nee_max']))
    p['recovery_status']=np.where(p.original_nee.notna(),'original_known_unchanged',np.where(p.recovery_candidate.gt(P['nee_max']),'recovery_exceeds_5_unknown',np.where(p.recovered_nee.notna(),'recovered','unresolved')))
    p['nee_with_recovery']=p.original_nee.fillna(p.recovered_nee)
    assert p.cid.is_unique
    assert np.all(p.loc[p.original_nee.notna(),'nee_with_recovery']==p.loc[p.original_nee.notna(),'original_nee'])
    p.to_csv(RES/'eicu_recovered_pairs.csv',index=False)
    output=[]
    def add(metric,value,denom=None):output.append(dict(metric=metric,value=value,denominator=denom))
    add('primary_pairs',len(p));add('primary_patients',p.pid.nunique());add('original_known_pairs',p.original_nee.notna().sum(),len(p));add('original_unknown_pairs',p.original_nee.isna().sum(),len(p))
    add('new_recovered_pairs',p.recovered_nee.notna().sum(),p.original_nee.isna().sum());add('patients_with_recovered_pairs',p.loc[p.recovered_nee.notna(),'pid'].nunique())
    add('known_pairs_with_recovery_raw',p.nee_with_recovery.notna().sum(),len(p));add('candidate_pairs_above5_set_unknown',p.recovery_candidate.gt(5).sum())
    add('unknown_pairs_remaining',p.nee_with_recovery.isna().sum(),len(p));add('changed_original_known_pairs',(p.original_nee.notna() & p.nee_with_recovery.ne(p.original_nee)).sum())
    add('source_candidate_segments',len(seg));add('segments_recovered',seg.new_dose.notna().sum());add('source_medication_rows',len(med));add('explicit_valid_medication_rows',med.valid.sum())
    add('candidate_raw_rows',len(cand));add('candidate_valid_weight_rows',cand.weight.notna().sum(),len(cand));add('candidate_unknown_weight_rows',cand.weight.isna().sum(),len(cand))
    add('candidate_weight_out_of_range_admission_rows',((pd.to_numeric(rows.admission_weight,errors='coerce')<30)|(pd.to_numeric(rows.admission_weight,errors='coerce')>300)).sum(),len(rows))
    add('row_count_conserved',int(len(p)==c.execute('SELECT count(*) FROM pairs WHERE is_primary').fetchone()[0]));add('all_original_known_unchanged',1)
    save(pd.DataFrame(output),'phase1_eicu_recovery_summary.csv')
    save(cand.groupby(['drug','method','reason']).agg(rows=('gkey','size'),segments=('gkey','nunique')).reset_index(),'phase1_eicu_recovery_reasons.csv')
    save(p.groupby('recovery_status').agg(pairs=('cid','size'),patients=('pid','nunique')).reset_index(),'phase1_eicu_recovery_status.csv')
    aud=pd.concat([cand[cand.candidate_dose.notna()].sample(n=min(50,cand.candidate_dose.notna().sum()),random_state=P['seed']),cand[cand.candidate_dose.isna()].sample(n=min(20,cand.candidate_dose.isna().sum()),random_state=P['seed'])]).copy()
    aud.to_csv(RES/'recovery_audit_expected.csv',index=False)
    aud[['infusiondrugid']].drop_duplicates().to_csv(RES/'recovery_audit_ids.csv',index=False)
    print(json.dumps({r['metric']:int(r['value']) for r in output},ensure_ascii=False))
    c.close()

if __name__=='__main__':run()
