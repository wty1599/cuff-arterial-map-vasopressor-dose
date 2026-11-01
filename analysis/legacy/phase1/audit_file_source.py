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
"""Independent raw-file audit. Does not import or query production transformations."""
from pathlib import Path
import sys,csv,json,math,struct,re,traceback,datetime,hashlib

import duckdb,pandas as pd
from audit_quality import P as QUALITY, anthro
ROOT=Path('@CUFF_WORKSPACE@');OUT=ROOT/'phase1';R=OUT/'_restricted';MAPS=ROOT/'phase0';SIC=Path('@SICDB_ROOT@')
def reg(c,n,d):
    d=d.copy()
    for k in d:
        if isinstance(d[k].dtype,pd.StringDtype):d[k]=d[k].astype(object).where(d[k].notna(),None)
    c.register('_input_frame',d);c.execute(f'CREATE OR REPLACE TABLE {n} AS SELECT * FROM _input_frame');c.unregister('_input_frame')
def rd(p,latin=False):return f"read_csv_auto('{Path(p).as_posix()}',all_varchar=true"+(",encoding='latin-1'" if latin else '')+")"
def numeric(v):
    try:z=float(v);return z if math.isfinite(z) else None
    except (TypeError,ValueError):return None
def mid(v):
    t=re.fullmatch(r'\s*(\d+(?:\.\d+)?)-(\d+(?:\.\d+)?)\s*',str(v));return (float(t[1])+float(t[2]))/2 if t else None
def regions(rows):
    ans=[]
    for left,right in sorted(rows):
        if ans and left<=ans[-1][1]+3600:ans[-1]=(ans[-1][0],max(ans[-1][1],right))
        else:ans.append((left,right))
    return ans
def convert(v,mass,timeunit,perkg,weight,agent,quality=False):
    if quality and (weight is None or not QUALITY["weight_kg_range"][0]<=weight<=QUALITY["weight_kg_range"][1]):weight=None
    v=numeric(v)
    if v is None or v<=0:return None
    if timeunit in ('min','minute'):scale=1
    elif timeunit in ('uur','hour','hr'):scale=60
    else:return None
    if agent=='vasopressin':
        if mass.lower() not in ('ie','u','unit','units'):return None
        if perkg:
            if weight is None or weight<=0:return None
            v*=weight
        return v/scale*2.5
    fact={'g':1e6,'mg':1000,'µg':1,'μg':1,'mcg':1,'ng':.001}.get(mass)
    if fact is None:return None
    if not perkg:
        if weight is None or weight<=0:return None
        v/=weight
    drugfactor={'norepinephrine':1,'epinephrine':1,'phenylephrine':.1,'dopamine':.01,'angiotensin_ii':10}[agent]
    return v*fact/scale*drugfactor
def main(db):
    started_at=datetime.datetime.now(datetime.timezone.utc).isoformat()
    c=duckdb.connect(str(R/f'audit_file_{db}.duckdb'));c.execute("SET memory_limit='700MB';SET threads=1;")
    c.execute(f"SET temp_directory='{R.as_posix()}/audit_{db}_tmp';")
    inp=pd.read_csv(R/f'audit_phase1_{db}.csv',dtype={'sid':str});reg(c,'checks',inp)
    bm=pd.read_csv(MAPS/f'map_bp_{db}.csv',dtype=str).fillna('');reg(c,'concepts',bm[bm.selected.eq('1')]);vm=pd.read_csv(MAPS/f'map_vaso_{db}.csv',dtype=str).fillna('');vm=vm[vm.selected.eq('1')];reg(c,'vaso_concepts',vm)
    saved_sql=[]
    def query(q):saved_sql.append(q+';');return c.execute(q)
    if db=='amsterdam':
        query(f"CREATE OR REPLACE TABLE source_bp AS SELECT x.admissionid sid,CAST(x.measuredat AS DOUBLE)/1000 ts,m.classification modality,m.component,m.family,TRY_CAST(x.value AS DOUBLE) measurement FROM {rd(MAPS/'_restricted'/'amsterdam_numeric_source.csv',True)} x JOIN concepts m ON x.itemid=m.itemid JOIN (SELECT sid,min(t)-300 lo,max(t)+300 hi FROM checks GROUP BY sid) z ON x.admissionid=z.sid AND CAST(x.measuredat AS DOUBLE)/1000 BETWEEN z.lo AND z.hi")
        meds=query(f"SELECT x.*,m.drug FROM {rd('@AMSTERDAM_ROOT@/drugitems.csv',True)} x JOIN vaso_concepts m ON x.itemid=m.itemid WHERE x.admissionid IN (SELECT sid FROM checks) AND x.iscontinuous='1' AND lower(x.ordercategory) LIKE '%spuitpompen%'").df()
        admissions=query(f"SELECT admissionid sid,admittedat,dischargedat,weightgroup,heightgroup FROM {rd('@AMSTERDAM_ROOT@/admissions.csv',True)} WHERE admissionid IN (SELECT sid FROM checks)").df();starts={str(x.sid):float(x.admittedat)/1000 for x in admissions.itertuples()};ends={str(x.sid):float(x.dischargedat)/1000 for x in admissions.itertuples()};weights={str(x.sid):mid(x.weightgroup) for x in admissions.itertuples()};heights={str(x.sid):mid(x.heightgroup) for x in admissions.itertuples()}
        raw=[];clean_raw=[]
        for x in meds.to_dict('records'):
            sid=str(x['admissionid']);left=numeric(x['start']);right=numeric(x['stop']);r=numeric(x['dose'])
            if left is None or right is None or r is None or r<=0:continue
            left=max(left/1000,starts[sid]);right=min(right/1000,ends[sid])
            if right<=left:continue
            value=convert(r,x['doseunit'],x['doserateunit'],str(x['doserateperkg'])=='1',weights[sid],x['drug']);raw.append((sid,x['drug'],left,right,value));clean_value=convert(r,x['doseunit'],x['doserateunit'],str(x['doserateperkg'])=='1',weights[sid],x['drug'],quality=True);clean_raw.append((sid,x['drug'],left,right,clean_value))
    else:
        if '--resume-blocks' in sys.argv:
            blocks=c.execute('SELECT * FROM audit_source_blocks').df()
        else:
            blocks=query(f"SELECT DISTINCT x.id,x.CaseID sid,x.DataID,x.Offset,x.rawdata,m.classification modality,m.component,m.family FROM {rd(SIC/'data_float_h.csv.gz')} x JOIN concepts m ON x.DataID=m.itemid JOIN checks z ON x.CaseID=z.sid AND TRY_CAST(x.Offset AS DOUBLE)<=z.t+300 AND TRY_CAST(x.Offset AS DOUBLE)+3600>z.t-300").df()
            reg(c,'audit_source_blocks',blocks)
        decoded=[]
        for x in blocks.to_dict('records'):
            h=x['rawdata']
            if not isinstance(h,str) or len(h)!=482:continue
            b=bytes.fromhex(h[2:])
            for minute in range(60):
                piece=b[4*minute:4*minute+4]
                if piece==b'\x00\x00\x00\x00':continue
                v=struct.unpack('<f',piece)[0]
                if math.isfinite(v):decoded.append((str(x['sid']),float(x['Offset'])+minute*60,x['modality'],x['component'],x['family'],v))
        reg(c,'source_bp',pd.DataFrame(decoded,columns=['sid','ts','modality','component','family','measurement']))
        meds=query(f"SELECT x.*,m.drug,m.unitname FROM {rd(SIC/'medication.csv.gz')} x JOIN vaso_concepts m ON x.DrugID=m.itemid WHERE x.CaseID IN (SELECT sid FROM checks) AND x.IsSingleDose='0'").df()
        admissions=query(f"SELECT CaseID sid,ICUOffset,TimeOfStay,TRY_CAST(WeightOnAdmission AS DOUBLE)/1000 AS weight_kg,HeightOnAdmission FROM {rd(SIC/'cases.csv.gz')} WHERE CaseID IN (SELECT sid FROM checks)").df();starts={str(x.sid):numeric(x.ICUOffset) for x in admissions.itertuples()};ends={str(x.sid):float(x.TimeOfStay) for x in admissions.itertuples()};weights={str(x.sid):numeric(x.weight_kg) for x in admissions.itertuples()};heights={str(x.sid):numeric(x.HeightOnAdmission) for x in admissions.itertuples()}
        raw=[];clean_raw=[]
        for x in meds.to_dict('records'):
            sid=str(x['CaseID']);left=numeric(x['Offset']);right=numeric(x['OffsetDrugEnd']);r=numeric(x['AmountPerMinute'])
            if left is None or right is None or r is None or r<=0 or starts[sid] is None:continue
            left=max(left,starts[sid]);right=min(right,ends[sid])
            if right<=left:continue
            value=convert(r,x['unitname'],'min',False,weights[sid],x['drug']);raw.append((sid,x['drug'],left,right,value));clean_value=convert(r,x['unitname'],'min',False,weights[sid],x['drug'],quality=True);clean_raw.append((sid,x['drug'],left,right,clean_value))
    meds.to_csv(R/f'audit_source_{db}_medication.csv',index=False)
    admissions.to_csv(R/f'audit_source_{db}_admissions.csv',index=False)
    reg(c,'source_drugs',pd.DataFrame(raw,columns=['sid','drug','begin','finish','dose']))
    reg(c,'source_clean_drugs',pd.DataFrame(clean_raw,columns=['sid','drug','begin','finish','dose']))
    clean_components={}
    for sid,drug,left,right,dose in clean_raw:clean_components.setdefault((sid,drug),[]).append((left,right,dose))
    intervals={};components={}
    for sid,drug,left,right,dose in raw:
        intervals.setdefault(sid,[]).append((left,right));components.setdefault((sid,drug),[]).append((left,right,dose))
    exposure={s:regions(v) for s,v in intervals.items()};comp_exposure={k:regions([(a,b) for a,b,_ in v]) for k,v in components.items()}
    reconstructed=[]
    for x in inp.to_dict('records'):
        sid=x['sid'];t=float(x['t']);episode=next((z for z in exposure.get(sid,[]) if z[0]<=t<z[1]),None);expected=[drug for (s,drug),v in comp_exposure.items() if s==sid and any(a<=t<b for a,b in v)];active=[(drug,dose) for (s,drug),v in components.items() if s==sid for a,b,dose in v if a<=t<b]
        ne=sum(v for _,v in active) if expected and set(expected)=={d for d,_ in active} and all(v is not None for _,v in active) else None
        clean_active=[(drug,dose) for (s,drug),v in clean_components.items() if s==sid for a,b,dose in v if a<=t<b]
        clean_ne=sum(v for _,v in clean_active) if expected and set(expected)=={d for d,_ in clean_active} and all(v is not None for _,v in clean_active) else None
        if clean_ne is not None and clean_ne>QUALITY['nee_max']:clean_ne=None
        av=anthro(weights.get(sid),heights.get(sid))
        reconstructed.append((x['audit_id'],sid,t,episode[0] if episode else None,episode[1] if episode else None,ne,clean_ne,*av.values()))
    reg(c,'reconstructed',pd.DataFrame(reconstructed,columns=['audit_id','sid','t','source_start','source_stop','source_nee','source_clean_nee','source_weight','source_height','source_clean_weight','source_clean_height','source_clean_bmi']))
    query("CREATE OR REPLACE TABLE source_triplets AS SELECT sid,ts,modality,family,quantile_cont(measurement,.5) FILTER(WHERE component='map') pressure,quantile_cont(measurement,.5) FILTER(WHERE component='sbp') systolic,quantile_cont(measurement,.5) FILTER(WHERE component='dbp') diastolic FROM source_bp GROUP BY sid,ts,modality,family")
    query("CREATE OR REPLACE TABLE source_clean AS SELECT sid,ts,modality,quantile_cont(pressure,.5) pressure FROM source_triplets WHERE pressure>=20 AND pressure<=200 AND (systolic IS NULL OR systolic<=300) AND (systolic IS NULL OR diastolic IS NULL OR (systolic>diastolic AND systolic-diastolic>=10)) GROUP BY sid,ts,modality")
    actual=query("SELECT r.audit_id,r.source_start,r.source_stop,r.source_nee,r.source_clean_nee,r.source_weight,r.source_height,r.source_clean_weight,r.source_clean_height,r.source_clean_bmi,(SELECT quantile_cont(b.pressure,.5) FROM source_clean b WHERE b.sid=r.sid AND b.modality='cuff' AND abs(b.ts-r.t)<.000001) source_cuff,(SELECT quantile_cont(b.pressure,.5) FROM source_clean b WHERE b.sid=r.sid AND b.modality='arterial' AND b.ts BETWEEN r.t-300 AND r.t+300 AND b.ts>=r.source_start AND b.ts<r.source_stop) source_art,(SELECT min(abs(b.ts-r.t)) FROM source_clean b WHERE b.sid=r.sid AND b.modality='arterial' AND b.ts BETWEEN r.t-300 AND r.t+300 AND b.ts>=r.source_start AND b.ts<r.source_stop) source_diff FROM reconstructed r").df()
    result=inp.merge(actual,on='audit_id');checks={'cuff_map':'source_cuff','art_map':'source_art','nee':'source_nee','nearest_diff':'source_diff','episode_start':'source_start','episode_stop':'source_stop','expected_phase0_weight':'source_weight','expected_phase0_height':'source_height'}
    for old,new in checks.items():
        a=pd.to_numeric(result[old],errors='coerce');b=pd.to_numeric(result[new],errors='coerce');tol=1e-4 if old in ('nearest_diff','episode_start','episode_stop') else 1e-6*(1+a.abs());result[old+'_match']=(a.isna()&b.isna())|((a-b).abs()<=tol)
    matchcols=[k+'_match' for k in checks];result['all_match']=result[matchcols].all(axis=1);result.to_csv(R/f'audit_file_comparison_{db}.csv',index=False)
    summary=[{'metric':'sampled_pairs','value':len(result)},{'metric':'any_discrepancy_pairs','value':int((~result.all_match).sum())}]+[{'metric':key+'_discrepancies','value':int((~result[key+'_match']).sum())} for key in checks]
    pd.DataFrame(summary).to_csv(OUT/f'phase1_source_audit_raw_{db}.csv',index=False)
    (ROOT/'sql'/'phase1'/f'audit_{db}_source_queries.sql').write_text('\n\n'.join(saved_sql),encoding='utf-8')
    manifest={'database':db,'sample_pairs':len(result),'discrepancy_pairs':int((~result.all_match).sum()),'source_queries_executed':True,'source_clean_unknown_pairs':int(result.source_clean_nee.isna().sum()),'source_query_started_at':started_at,'source_query_completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'timing':'阶段0结果之后；新50对样本固定后执行。阶段1主分层结果相对时序见 audit_execution_timeline.csv','phase1_clean_output_comparison':'pending','sample_sha256':hashlib.sha256((R/f'audit_phase1_{db}.csv').read_bytes()).hexdigest(),'audit_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (R/f'audit_summary_{db}.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(manifest,ensure_ascii=False));c.close()
if __name__=='__main__':
    try:main(sys.argv[1])
    except Exception as e:
        (R/f'audit_file_error_{sys.argv[1]}.txt').write_text(traceback.format_exc(),encoding='utf-8');print(json.dumps({'db':sys.argv[1],'error_class':type(e).__name__,'stage':'independent_file_audit'}));sys.exit(1)
