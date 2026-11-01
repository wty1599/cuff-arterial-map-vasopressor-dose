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
import re, hashlib, datetime
from importlib import import_module
sys.path.append(str(ROOT/'sql'))
original_convert_rate=import_module('08_core').convert_rate

def num(x): return pd.to_numeric(x,errors='coerce')
def bad(x,limits): return x.notna() & ~x.between(*limits)
def read_raw(p): return pd.read_csv(p,dtype=str,keep_default_na=False)

def source_weight_flags(db,c,df,st):
    # Independent phase1 range screen of the weights actually required for dose conversion.
    if db=='sicdb':
        ids=st.loc[bad(st.weight,[30,300]),'sid']
        z=c.execute("SELECT sid,start,\"end\" stop FROM drugs WHERE drug<>'vasopressin' AND dose IS NOT NULL").df()
        z=z[z.sid.astype(str).isin(ids)].copy();z['sid']=z.sid.astype(str)
    elif db=='mimic':
        r=read_raw(SRC/'mimic_drugs.csv')
        r['sid']=r.stay_id.astype(str)
        fw=st.set_index('sid').dose_fallback_weight
        w=num(r.patientweight); w=w.where(w>0,r.sid.map(fw))
        u=r.rateuom.str.lower().str.replace(' ','',regex=False)
        mp=read_raw(ROOT/'phase0'/'map_vaso_mimic.csv');drug=r.itemid.map(mp.set_index('itemid').drug)
        perkg=u.str.contains('/kg',regex=False);vaso=drug.eq('vasopressin')
        needw=(~perkg&~vaso)|(perkg&vaso)
        mask=bad(w,[30,300])&needw&(num(r.rate)>0)&r.ordercategorydescription.eq('Continuous Med')&~r.statusdescription.isin(['Rewritten','Bolus'])
        cand=r.index[mask]
        if len(cand): mask.loc[cand]=np.array([np.isfinite(original_convert_rate(r.loc[i,'rate'],r.loc[i,'rateuom'],w.loc[i],drug.loc[i])[0]) for i in cand],dtype=bool)
        z=r[mask].copy();z['start']=pd.to_datetime(z.starttime).astype('datetime64[s]').astype('int64');z['stop']=pd.to_datetime(z.endtime).astype('datetime64[s]').astype('int64')
    elif db=='eicu':
        r=read_raw(SRC/'eicu_drugs.csv'); mp=read_raw(ROOT/'phase0'/'map_vaso_eicu.csv') if (ROOT/'phase0'/'map_vaso_eicu.csv').exists() else read_raw(ROOT/'phase0'/'map_vasopressor_eicu.csv')
        # Source mapping filename is detected explicitly below if needed.
        r=r[r.drugname.isin(mp.loc[mp.selected.eq('1'),'itemid'])].copy();r['sid']=r.patientunitstayid.astype(str)
        r['drug']=r.drugname.map(mp.set_index('itemid').drug)
        u=r.drugname.str.lower(); w=r.sid.map(st.set_index('sid').weight)
        perkg=u.str.contains('/kg',regex=False);vaso=r.drug.eq('vasopressin')
        needw=(~perkg&~vaso)|(perkg&vaso)
        r['bad']=bad(w,[30,300])&needw&(num(r.drugrate)>0)
        cand=r.index[r.bad]
        def source_unit(name):
            v=[s for s in re.findall(r'\(([^()]*)\)',name) if re.search(r'/min|/hr|/hour',s,re.I)]
            return v[-1] if v else ''
        if len(cand): r.loc[cand,'bad']=np.array([np.isfinite(original_convert_rate(r.loc[i,'drugrate'],source_unit(r.loc[i,'drugname']),w.loc[i],r.loc[i,'drug'])[0]) for i in cand],dtype=bool)
        r['start']=num(r.infusionoffset)*60
        z=r.groupby(['sid','drug','start'],as_index=False)['bad'].any().sort_values(['sid','drug','start'])
        z['stop']=np.minimum(z.groupby(['sid','drug']).start.shift(-1).fillna(np.inf),z.start+3600)
        z=z[z.bad].copy()
    else:
        r=read_raw(SRC/'amsterdam_drugs.csv.gz');r['sid']=r.admissionid.astype(str)
        w=r.sid.map(st.set_index('sid').weight)
        mp=read_raw(ROOT/'phase0'/'map_vaso_amsterdam.csv');drug=r.itemid.map(mp.set_index('itemid').drug)
        perkg=num(r.doserateperkg).eq(1);vaso=drug.eq('vasopressin')
        needw=(~perkg&~vaso)|(perkg&vaso)
        mask=bad(w,[30,300])&needw&(num(r.dose)>0)&num(r.iscontinuous).eq(1)
        cand=r.index[mask]
        if len(cand): mask.loc[cand]=np.array([np.isfinite(original_convert_rate(r.loc[i,'dose'],r.loc[i,'doseunit']+('/kg/' if perkg.loc[i] else '/')+r.loc[i,'doserateunit'],w.loc[i],drug.loc[i])[0]) for i in cand],dtype=bool)
        z=r[mask].copy()
        z['start']=num(z.start)/1000;z['stop']=num(z.stop)/1000
    if len(z)==0: return pd.Series(False,index=df.index)
    zi=z[['sid','start','stop']].copy();pi=df[['cid','sid','t']].copy()
    zi['sid']=zi.sid.astype(object);pi['sid']=pi.sid.astype(object)
    c.register('invalid_weight_intervals',zi);c.register('phase1_anchors',pi)
    ids=c.execute('SELECT DISTINCT p.cid FROM phase1_anchors p JOIN invalid_weight_intervals z ON p.sid=z.sid AND p.t>=z.start AND p.t<z.stop').df().cid
    c.unregister('invalid_weight_intervals');c.unregister('phase1_anchors')
    return df.cid.isin(ids)

allqc=[];counts=[];coverage=[];missing=[];prev=[];hashes=[]
for path in [Path(__file__).parent/'phase1_parameters.json',ROOT/'sql'/'phase0_parameters.json']:
    hashes.append({'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'snapshot_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
save(pd.DataFrame(hashes),'phase1_input_manifest.csv')
for db in DBS:
    c=con(db); st=c.execute('SELECT * FROM stays').df();st['sid']=st.sid.astype(str)
    d=c.execute('SELECT p.*,s.age,s.weight weight_raw,s.height height_raw FROM pairs p JOIN stays s USING(sid) WHERE p.is_primary').df();d['sid']=d.sid.astype(str);d['pid']=d.pid.astype(str)
    assert d.cid.is_unique
    d['database']=db;d['patient']=db+':'+d.pid;d['bmi_raw']=d.bmi;d['nee_raw']=d.nee
    for var,lims in [('weight',[30,300]),('height',[120,230]),('bmi',[12,80])]:
        d[var+'_flag']=bad(d[var+'_raw'],lims)
        d[var]=d[var+'_raw'].mask(d[var+'_flag'])
    d['bmi']=d.weight/(d.height/100)**2;d['bmi']=d.bmi.mask(bad(d.bmi,[12,80]))
    d['nee_flag']=d.nee.notna() & ((d.nee>5)|(d.nee<0))
    d['dose_weight_flag']=source_weight_flags(db,c,d,st)
    d['nee']=d.nee.mask(d.nee_flag|d.dose_weight_flag)
    d['quality_any']=d[['weight_flag','height_flag','bmi_flag','nee_flag','dose_weight_flag']].any(axis=1)
    d['dose_stratum']=strata(d.nee)
    d['bmi40']=np.where(d.bmi.notna(),(d.bmi>40).astype(float),np.nan)
    if db in ['mimic','eicu','amsterdam']:
        rr=read_raw(SRC/f'{db}_stays.csv');key={'mimic':'stay_id','eicu':'patientunitstayid','amsterdam':'admissionid'}[db]
        gender=rr.set_index(key).gender.str.lower().str.strip()
        d['sex']=d.sid.map(gender).map({'m':'M','male':'M','man':'M','f':'F','female':'F','vrouw':'F'})
    else: d['sex']=pd.NA
    af=pd.to_numeric(d.af.astype(int),errors='coerce')
    if db in ['mimic','eicu']:
        dx=read_raw(RES/f'{db}_diagnosis_coverage.csv');covered=set(dx.loc[dx.diagnosis_covered.isin(['t','true','True','1']),'sid'])
        d['af_observed']=d.af_source.eq('rhythm')|d.af_source.eq('diagnosis')|d.sid.isin(covered)
        d['af_model']=af.where(d.af_observed)
    else: d['af_observed']=d.af_source.eq('rhythm');d['af_model']=af.where(d.af_observed)
    d['af_subgroup']=af.where(d.af_source.eq('rhythm')) if db in ['mimic','amsterdam'] else d.af_model if db=='eicu' else np.nan
    d.to_parquet(RES/f'pairs_{db}.parquet',index=False)
    # Window sensitivities reuse the same cuff anchors and dose assignment.
    w2=c.execute('SELECT cid,art_map,art_values,nearest_diff FROM pairs_windows p JOIN stays s USING(sid) WHERE s.is_primary AND win=2').df()
    w2.to_parquet(RES/f'window2_{db}.parquet',index=False)
    sm=c.execute('''SELECT cf.cid,median(a.map_value) art_map,count(*) art_values,min(abs(a.t-cf.t)) nearest_diff
        FROM cuff60 cf JOIN stays s USING(sid) JOIN arterial a ON a.sid=cf.sid AND a.episode_start=cf.episode_start
        AND floor(a.t/60)=floor(cf.t/60) WHERE s.is_primary GROUP BY cf.cid''').df()
    sm.to_parquet(RES/f'same_minute_{db}.parquet',index=False)
    cv=c.execute('''WITH cs AS (SELECT cf.sid,episode_start,min(t) lo,max(t) hi FROM cuff60 cf JOIN stays s USING(sid) WHERE s.is_primary GROUP BY cf.sid,episode_start),
        ar AS (SELECT a.sid,episode_start,min(t) lo,max(t) hi FROM arterial a JOIN stays s USING(sid) WHERE s.is_primary GROUP BY a.sid,episode_start),
        spans AS (SELECT cs.sid,cs.episode_start,greatest(cs.lo,ar.lo) lo,least(cs.hi,ar.hi) hi FROM cs JOIN ar USING(sid,episode_start)),
        anchors AS (SELECT cf.cid,cf.sid FROM cuff60 cf JOIN spans s USING(sid,episode_start) WHERE cf.t BETWEEN s.lo AND s.hi)
        SELECT count(*) denominator,count(p.cid) numerator,count(DISTINCT a.sid) stays FROM anchors a LEFT JOIN pairs p USING(cid)''').df().iloc[0]
    coverage.append(dict(database=db,**cv.to_dict(),coverage=cv.numerator/cv.denominator if cv.denominator else np.nan,definition='shared_observed_span'))
    for v in ['weight','height','bmi','nee']:
        fl=d[v+'_flag'];clean=d[v]
        allqc.append(dict(database=db,variable=v,pairs=len(d),nonmissing_before=int(d[v+'_raw'].notna().sum()),out_of_range_pairs=int(fl.sum()),out_of_range_patients=d.loc[fl,'pid'].nunique(),nonmissing_after=int(clean.notna().sum()),missing_after=int(clean.isna().sum())))
    allqc.append(dict(database=db,variable='dose_conversion_weight',pairs=len(d),out_of_range_pairs=int(d.dose_weight_flag.sum()),out_of_range_patients=d.loc[d.dose_weight_flag,'pid'].nunique()))
    for level in [*LABELS,'unknown']:
        z=d[d.dose_stratum.eq(level)];counts.append(dict(database=db,stratum=level,pairs=len(z),patients=z.pid.nunique()))
    assert sum(x['pairs'] for x in counts if x['database']==db)==len(d)
    patient=d.drop_duplicates('pid')
    for v in ['age','sex','weight','height','bmi','bmi40','af_model','af_subgroup']:
        a=patient[v]; av=pd.to_numeric(a,errors='coerce') if v!='sex' else pd.Series(dtype=float)
        missing.append(dict(database=db,variable=v,unit='paired_patients',n=len(patient),missing=int(a.isna().sum()),missing_fraction=a.isna().mean(),median=av.median() if len(av.dropna()) else np.nan,q1=av.quantile(.25) if len(av.dropna()) else np.nan,q3=av.quantile(.75) if len(av.dropna()) else np.nan))
    for v in ['bmi40','af_model','af_subgroup']:
        a=patient[v];prev.append(dict(database=db,variable=v,positive=int(a.eq(1).sum()),negative=int(a.eq(0).sum()),unknown=int(a.isna().sum()),patients=len(patient)))
    print(json.dumps(dict(database=db,main_pairs=len(d),patients=d.pid.nunique(),window2_pairs=len(w2),same_minute_pairs=len(sm),known_nee=int(d.nee.notna().sum()),quality_flagged_pairs=int(d.quality_any.sum()))),flush=True)
    c.close()
save(pd.DataFrame(allqc),'phase1_quality_ranges.csv');save(pd.DataFrame(counts),'phase1_flow.csv');save(pd.DataFrame(coverage),'phase1_pair_coverage.csv');save(pd.DataFrame(missing),'phase1_covariate_missingness.csv');save(pd.DataFrame(prev),'phase1_covariate_prevalence.csv')
print('PREPARATION_COMPLETE',flush=True)
