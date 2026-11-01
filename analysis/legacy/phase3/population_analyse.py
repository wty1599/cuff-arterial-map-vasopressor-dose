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

from pathlib import Path
import os,sys,json,datetime,hashlib
os.environ.setdefault('OPENBLAS_NUM_THREADS','2');os.environ.setdefault('OMP_NUM_THREADS','2')
ROOT=Path('@CUFF_WORKSPACE@');OUT=ROOT/'phase3';RES=OUT/'_restricted';CODE=ROOT/'sql/phase3'

import numpy as np,pandas as pd
LABELS=['<=0.1','>0.1-0.25','>0.25-0.5','>0.5-1.0','>1.0']
UNITS={'mimic':['Cardiac Vascular Intensive Care Unit (CVICU)'],'eicu':['CSICU','CTICU','CCU-CTICU']}
P={
 'written_at':datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
 'specified_before_S1_S2_results':True,'definitions_from_phase1_phase2_unchanged':True,
 'S1':{'name':'Exclude cardiac/cardiovascular surgical index ICU','source':{'mimic':'mimiciv_icu.icustays.first_careunit','eicu':'public.patient.unittype'},'excluded_index_units':UNITS,'missing_index_unit':'unknown; not retained in S1','supplementary_databases':{'amsterdam':'not_performed: admissions.location identifies IC/MC, not cardiac-surgical unit; specialty is not an ICU unit and is not substituted','sicdb':'not_performed: source unit dictionary distinguishes care intensity, not cardiac-surgical ICU; surgical diagnosis/admission type is not a unit and is not substituted'}},
 'S2':{'name':'Restrict to source-confirmed sepsis population','mimic':'At least one existing mimiciv_derived.sepsis3 row with sepsis3=true for index stay, without imposing a new onset-time window','mimic_absence':'Existing derived table contains positive rows only. No matching row is not source-confirmed non-sepsis; report as not listed in positive-derived table and exclude from the restriction. Null/false rows, if present, remain separately classified.','eicu':'Case-insensitive substring sepsis in nonempty public.patient.apacheadmissiondx; source field is APACHE admission diagnosis','eicu_missing':'Null/empty/unknown diagnosis is unascertainable; nonempty non-sepsis diagnosis is recorded other APACHE diagnosis, not a claim that infection is absent','supplementary_databases':'not_requested'},
 'scope':'Main endpoint overall and five original-known dose strata; deep missed hypotension arterial<55 overall and >1.0; arterial<65 given cuff>=65 overall and >1.0. Add layer arterial<65 prevalence, counts, source ascertainment, flow checks. No additional tests, models, pooled estimates, outcomes, or phenotype analyses.',
 'units':{'dose':'micrograms/kg/min norepinephrine equivalent','MAP':'mmHg'},
 'dose_strata':LABELS,'overall':'All selected pairs including unknown dose','known_dose':'Original phase1 range-checked known doses; no recovery',
 'bootstrap':{'replicates':2000,'seed':20261027,'unit':'patient within selected database and sensitivity population','statistic':'pair-weighted numerator/denominator','interval':'percentile 95%','implementation':'Phase1 summarize and phase2 ratio_engine resampling mechanics retained: sorted patient index, multinomial whole-patient weights in batches of50, no pair resampling','boundary':'Keep event and denominator counts; mark degenerate zero/all-event CI as boundary_not_informative; empty denominator not_estimable','precision':'denominator<30 or CI width>0.30 marked insufficient, never combine/delete'},
 'audit':'20 patients per main database, deterministically sampled over joint S1/S2 states. Independently query source careunit and diagnosis/Sepsis-3 rows with separate SQL; row-level records retained in _restricted.',
 'earlier_work':'Phase2 baseline/range/pair audits reused for unchanged variables; no claim of new full extraction verification.'
}
parameter_path=OUT/'population_parameters.json'
if parameter_path.exists():
    old=json.loads(parameter_path.read_text(encoding='utf-8'));P['written_at']=old['written_at'];assert P==old,'Parameter drift: do not silently rewrite pre-result rules'
else:parameter_path.write_text(json.dumps(P,ensure_ascii=False,indent=2),encoding='utf-8')
(CODE/'population_parameters.json').write_text(json.dumps(P,ensure_ascii=False,indent=2),encoding='utf-8')
PARAM_HASH=hashlib.sha256(parameter_path.read_bytes()).hexdigest()

def save(d,name):d.to_csv(OUT/name,index=False,encoding='utf-8-sig')
def masks(d):
    return {'all':np.ones(len(d),bool),'known':d.nee.notna().to_numpy(),'unknown':d.nee.isna().to_numpy(),**{lab:d.dose_stratum.eq(lab).to_numpy() for lab in LABELS}}

def ratios(d,db,scenario):
    pats=pd.Index(sorted(d.patient.unique()));codes=pats.get_indexer(d.patient);N=len(pats)
    mm=masks(d);art=d.art_map.to_numpy();cuff=d.cuff_map.to_numpy();lo=art<65;normal=cuff>=65;specs=[];vectors=[]
    for metric,threshold,labels in [('missed_hypotension',65,['all',*LABELS]),('deep_missed_hypotension',55,['all','>1.0']),('art_lt65_given_cuff_ge65',None,['all','>1.0'])]:
        for lab in labels:
            mask=mm[lab];eligible=mask&((art<threshold) if threshold else normal);event=normal if threshold else lo
            dc=np.bincount(codes,weights=eligible,minlength=N);ec=np.bincount(codes,weights=eligible&event,minlength=N)
            specs.append(dict(database=db,sensitivity=scenario,metric=metric,stratum=lab,numerator=int((eligible&event).sum()),denominator=int(eligible.sum()),denominator_patients=int((dc>0).sum()),layer_pairs=int(mask.sum()),layer_patients=int(np.unique(codes[mask]).size),layer_low_arterial_pairs=int((mask&lo).sum()),layer_low_arterial_prevalence=float((mask&lo).sum()/mask.sum()) if mask.any() else np.nan,selected_pairs=len(d),selected_patients=N))
            vectors.extend([ec,dc])
    matrix=np.array(vectors).T;chunks=[];rng=np.random.default_rng(P['bootstrap']['seed'])
    if N:
        for b in range(0,2000,50):
            w=rng.multinomial(N,np.full(N,1/N),size=min(50,2000-b));chunks.append(w@matrix)
        count=np.concatenate(chunks);ne=count[:,0::2];de=count[:,1::2];br=np.divide(ne,de,out=np.full_like(ne,np.nan),where=de>0)
    else:br=np.full((2000,len(specs)),np.nan);ne=de=np.zeros_like(br)
    for i,row in enumerate(specs):
        good=np.isfinite(br[:,i]);ll,uu=np.quantile(br[good,i],[.025,.975]) if good.any() else (np.nan,np.nan)
        n=row['numerator'];den=row['denominator'];degenerate=bool(den and n in [0,den] and ll==uu)
        row.update(estimate=n/den if den else np.nan,ci_lower=ll,ci_upper=uu,bootstrap_valid=int(good.sum()),small_cell=den<30,wide_ci=bool(uu-ll>.30),boundary=degenerate,precision='not_estimable' if not den else 'boundary_not_informative' if degenerate else 'insufficient' if den<30 or uu-ll>.30 else 'adequate',ci_status='boundary_degenerate' if degenerate else 'ok' if den else 'no_denominator',ci_method='patient_cluster_percentile_2000',parameter_sha256=PARAM_HASH)
    np.savez_compressed(RES/f'population_bootstrap_{db}_{scenario}.npz',numerator=ne,denominator=de)
    return pd.DataFrame(specs)

results=[];flows=[];ascertain=[];checks=[];manifests=[];states=[]
for db in ['mimic','eicu']:
    path=ROOT/f'phase1/_restricted/pairs_{db}.parquet';d=pd.read_parquet(path);d.sid=d.sid.astype(str)
    pt=d[['sid','patient']].drop_duplicates();assert pt.sid.is_unique and pt.patient.is_unique
    rawpath=RES/f'population_source_{db}.csv';x=pd.read_csv(rawpath,dtype=str,keep_default_na=False);assert x.sid.is_unique
    assert set(pt.sid)==set(x.sid),'Source join lost/added paired patients'
    pt=pt.merge(x,on='sid',validate='one_to_one')
    unitcol='first_careunit' if db=='mimic' else 'unittype';unit=pt[unitcol].str.strip()
    pt['s1_status']=np.where(unit.eq(''),'unknown_index_unit',np.where(unit.isin(UNITS[db]),'excluded_cardiac_surgical_index_unit','retained_other_recorded_index_unit'))
    pt['S1']=pt.s1_status.eq('retained_other_recorded_index_unit')
    if db=='mimic':
        pos=pd.to_numeric(pt.sepsis_positive_rows);neg=pd.to_numeric(pt.sepsis_negative_rows);rows=pd.to_numeric(pt.sepsis_rows)
        pt['s2_status']=np.select([pos.gt(0),neg.gt(0),rows.gt(0)],['confirmed_sepsis3_positive','recorded_false_in_existing_sepsis3_table','null_only_in_existing_sepsis3_table'],default='not_listed_in_positive_sepsis3_table')
        pt['S2']=pos.gt(0)
    else:
        dx=pt.apacheadmissiondx.str.strip();unknown=dx.str.lower().isin(['','unknown','none','null','nan','other/unknown'])
        pt['s2_status']=np.select([~unknown&dx.str.contains('sepsis',case=False,regex=False),unknown],['recorded_sepsis_APACHE_diagnosis','unknown_APACHE_admission_diagnosis'],default='recorded_other_APACHE_diagnosis')
        pt['S2']=pt.s2_status.eq('recorded_sepsis_APACHE_diagnosis')
    pt.to_csv(RES/f'population_labels_{db}.csv',index=False)
    for scenario in ['S1','S2']:
        statuscol=scenario.lower()+'_status';allstates=pt.groupby(statuscol).patient.agg(['size'])
        for state,g in pt.groupby(statuscol):
            dd=d[d.sid.isin(g.sid)];ascertain.append(dict(database=db,sensitivity=scenario,status=state,patients=len(g),patient_fraction=len(g)/len(pt),pairs=len(dd),pair_fraction=len(dd)/len(d)))
        chosen=pt[pt[scenario]];s=d[d.sid.isin(chosen.sid)].copy()
        save(ratios(s,db,scenario),f'population_result_{db}_{scenario}.csv');results.append(pd.read_csv(OUT/f'population_result_{db}_{scenario}.csv'))
        for lab,mask in masks(s).items():
            z=s[mask];flows.append(dict(database=db,sensitivity=scenario,stratum=lab,source_patients=len(pt),source_pairs=len(d),retained_patients=z.patient.nunique(),retained_pairs=len(z),low_arterial_pairs=int(z.art_map.lt(65).sum()),deep_low_arterial_pairs=int(z.art_map.lt(55).sum()),normal_cuff_pairs=int(z.cuff_map.ge(65).sum())))
        checks.extend([dict(database=db,sensitivity=scenario,check='retained_plus_excluded_source_pairs',passed=len(s)+len(d[~d.sid.isin(chosen.sid)])==len(d)),dict(database=db,sensitivity=scenario,check='five_known_plus_unknown_equals_retained',passed=sum(s.dose_stratum.eq(v).sum() for v in LABELS)+s.nee.isna().sum()==len(s)),dict(database=db,sensitivity=scenario,check='labels_not_missing',passed=pt[scenario].notna().all()),dict(database=db,sensitivity=scenario,check='patient_to_stay_unique',passed=s.groupby('patient').sid.nunique().le(1).all())])
    # Equal initial allocation across observed joint eligibility states then deterministic fill to20.
    pt['joint_state']=pt.s1_status+' | '+pt.s2_status;sample=[]
    for _,g in pt.groupby('joint_state',sort=True):sample.append(g.sample(min(4,len(g)),random_state=20261027))
    aud=pd.concat(sample).drop_duplicates('sid')
    if len(aud)>20:aud=aud.sample(20,random_state=20261027)
    elif len(aud)<20:aud=pd.concat([aud,pt[~pt.sid.isin(aud.sid)].sample(min(20-len(aud),len(pt)-len(aud)),random_state=20261027)])
    aud.sort_values('sid').to_csv(RES/f'population_audit_expected_{db}.csv',index=False)
    aud[['sid']].sort_values('sid').to_csv(RES/f'population_audit_ids_{db}.csv',index=False)
    for file in [path,rawpath,OUT/'population_parameters.json',ROOT/'sql/phase1/02_rates.py',ROOT/'sql/phase2/01_analyses.py']:
        manifests.append(dict(database=db,file=file.name,sha256=hashlib.sha256(file.read_bytes()).hexdigest(),role='read_only_input_or_parameter'))
    print(json.dumps({'database':db,'population_analysis':'complete','S1_patients':int(pt.S1.sum()),'S2_patients':int(pt.S2.sum()),'source_patients':len(pt)}),flush=True)
allres=pd.concat(results,ignore_index=True)
assert len(allres)==40 and not allres.duplicated(['database','sensitivity','metric','stratum']).any()
for r in allres.itertuples():
    checks.append(dict(database=r.database,sensitivity=r.sensitivity,check=f'{r.metric}_{r.stratum}_numerator_denominator_order',passed=0<=r.numerator<=r.denominator<=r.layer_pairs))
    checks.append(dict(database=r.database,sensitivity=r.sensitivity,check=f'{r.metric}_{r.stratum}_CI_bounds',passed=not r.denominator or (0<=r.ci_lower<=r.ci_upper<=1 and 0<=r.estimate<=1)))
assert all(x['passed'] for x in checks)
save(allres,'phase3_S1_S2_results.csv');save(pd.DataFrame(flows),'population_flow.csv');save(pd.DataFrame(ascertain),'population_ascertainment.csv');save(pd.DataFrame(checks),'population_checks.csv');save(pd.DataFrame(manifests),'population_input_manifest.csv')
save(pd.DataFrame([dict(database=db,sensitivity=s,status='not_performed',reason=(P['S1']['supplementary_databases'][db] if s=='S1' else 'Not requested for supplementary databases')) for db in ['amsterdam','sicdb'] for s in ['S1','S2']]),'population_not_performed.csv')
print('S1_S2_COMPLETE_WITH_SOURCE_AUDIT_PENDING',flush=True)
