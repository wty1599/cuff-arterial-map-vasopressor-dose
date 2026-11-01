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

from importlib import import_module
g=import_module('08_core');globals().update({k:v for k,v in vars(g).items() if not k.startswith('_')})
def build(db):
    c=connect(db);dm=mapping(db,'vaso');dm=dm[dm.selected=='1'];dmap=dm.set_index('itemid')['drug'].to_dict()
    if db=='mimic':
        a=read(R/'mimic_stays.csv',dtype=str);s=pd.DataFrame({'sid':a.stay_id,'pid':a.subject_id,'age':num(a.age),'t0':epoch(a.intime),'t1':epoch(a.outtime),'first':num(a.patient_icu_rank)==1,'first_unknown':False,'care_eligible':True,'care_source':a.first_careunit,'hospital':'single','hadm':a.hadm_id})
        raw=read(R/'mimic_drugs.csv',dtype=str)
        # Admission measurements are independent of weights attached to infusions.
        cov=mapping(db,'covariates');save(c,'cov_map',cov)
        cq=csvq(R/'mimic_bp_cov.csv');sq=pd.DataFrame({'sid':s.sid,'t0':s.t0});save(c,'stay_times',sq)
        cv=c.execute(f"SELECT x.stay_id sid,epoch(try_cast(x.charttime AS TIMESTAMP)) t,m.label,try_cast(x.valuenum AS DOUBLE) val FROM {cq} x JOIN cov_map m ON x.itemid=m.itemid JOIN stay_times s ON x.stay_id=s.sid WHERE try_cast(x.valuenum AS DOUBLE)>0 AND epoch(try_cast(x.charttime AS TIMESTAMP))>=s.t0").df()
        cv['scaled']=cv.val*np.where(cv.label.eq('Admission Weight (lbs.)'),.45359237,np.where(cv.label.eq('Height'),2.54,1))
        cv['kind']=np.where(cv.label.str.contains('Weight'),'weight',np.where(cv.label.str.contains('Height'),'height','other'))
        for kind in ['weight']:
            q=cv[cv.kind.eq(kind)].sort_values(['sid','t']).drop_duplicates('sid');s[kind]=s.sid.map(q.set_index('sid').scaled)
        # Final D007 rule: linked admission height may be charted before ICU intake.
        hs=c.execute(f"""WITH x AS (
          SELECT x.stay_id sid,epoch(try_cast(x.charttime AS TIMESTAMP)) t,s.t0,m.label,
          try_cast(x.valuenum AS DOUBLE)*CASE WHEN m.label='Height' THEN 2.54 ELSE 1 END height
          FROM {cq} x JOIN cov_map m ON x.itemid=m.itemid JOIN stay_times s ON x.stay_id=s.sid
          WHERE m.label ILIKE '%height%' AND try_cast(x.valuenum AS DOUBLE)>0)
          SELECT sid,arg_min(height,struct_pack(distance := abs(t-t0),source_t := t,label_name := label)) height FROM x GROUP BY sid""").df()
        s['height']=s.sid.map(hs.set_index('sid').height)
        # Weight fallback: latest OMR weight on/before ICU admission, then earliest ICU weight.
        omr=read(R/'mimic_omr.csv',dtype=str);omr=omr[omr.result_name.str.contains('weight',case=False)]
        omr['wt']=num(omr.result_value)*np.where(omr.result_name.str.contains('lbs',case=False),.45359237,1)
        omr['ot']=epoch(omr.chartdate);j=omr.merge(s[['sid','pid','t0']],left_on='subject_id',right_on='pid');j=j[(j.ot<=j.t0)&(j.wt>0)].sort_values(['sid','ot']).drop_duplicates('sid',keep='last');s['dose_fallback_weight']=s.sid.map(j.set_index('sid').wt).fillna(s.weight)
        rhids=set(cov.loc[cov.label.eq('Heart Rhythm'),'itemid']);rid=','.join("'"+x+"'" for x in rhids)
        rh=c.execute(f"SELECT stay_id sid,bool_or(regexp_matches(lower(value),'atrial fib|a[- ]?fib|afib')) af,count(*) n FROM {cq} WHERE itemid IN ({rid}) AND nullif(trim(value),'') IS NOT NULL GROUP BY stay_id").df()
        dx=set(read(R/'mimic_af_dx.csv',dtype=str).hadm_id);afm=rh.set_index('sid').af;s['af']=s.sid.map(afm).fillna(s.hadm.isin(dx)).astype(bool);s['af_source']=np.where(s.sid.isin(afm.index),'rhythm',np.where(s.hadm.isin(dx),'diagnosis','not_observed'))
        raw['sid']=raw.stay_id;raw['drug']=raw.itemid.map(dmap);raw['start']=epoch(raw.starttime);raw['end']=epoch(raw.endtime);raw['rate_num']=num(raw.rate);raw['unit']=raw.rateuom
        raw['dose_weight']=num(raw.patientweight).where(num(raw.patientweight)>0,raw.sid.map(s.set_index('sid').dose_fallback_weight))
        raw['continuous']=raw.ordercategorydescription.eq('Continuous Med') & ~raw.statusdescription.isin(['Rewritten','Bolus']);raw['rawid']=raw.orderid+'|'+raw.itemid
        raw['source_unit']=raw.unit;raw['valid_rate_evidence']=raw.rate_num>0
    elif db=='eicu':
        a=read(R/'eicu_stays.csv',dtype=str);a['year']=num(a.hospitaldischargeyear);a['firstyear']=a.year.eq(a.groupby('uniquepid').year.transform('min'));early=a[a.firstyear];h=early.groupby('uniquepid').patienthealthsystemstayid.nunique();ambig=set(h[h>1].index)
        a['hoff']=num(a.hospitaladmitoffset);a['visit']=num(a.unitvisitnumber);firstids=set(a[a.firstyear&~a.uniquepid.isin(ambig)].sort_values(['uniquepid','hoff','visit','patientunitstayid'],ascending=[True,False,True,True]).drop_duplicates('uniquepid').patientunitstayid)
        s=pd.DataFrame({'sid':a.patientunitstayid,'pid':a.uniquepid,'age':num(a.age.str.replace('>','',regex=False)),'t0':0.,'t1':num(a.unitdischargeoffset)*60,'first':a.patientunitstayid.isin(firstids),'first_unknown':a.uniquepid.isin(ambig),'care_eligible':True,'care_source':a.unittype,'hospital':a.hospitalid,'weight':num(a.admissionweight),'height':num(a.admissionheight)})
        dx=set(read(R/'eicu_af_dx.csv',dtype=str).patientunitstayid);s['af']=s.sid.isin(dx);s['af_source']=np.where(s.af,'diagnosis','not_observed')
        raw=read(R/'eicu_drugs.csv',dtype=str);raw=raw[raw.drugname.isin(dmap)].copy();raw['sid']=raw.patientunitstayid;raw['drug']=raw.drugname.map(dmap);raw['start']=num(raw.infusionoffset)*60;raw['rate_num']=num(raw.drugrate);raw['dose_weight']=raw.sid.map(s.set_index('sid').weight);raw['rawid']=raw.infusiondrugid
        def unit(name):
            vals=re.findall(r'\(([^()]*)\)',name);vals=[v for v in vals if re.search(r'/min|/hr|/hour',v,re.I)];return vals[-1] if vals else ''
        raw['unit']=raw.drugname.map(unit);raw['source_unit']=raw.unit;raw['continuous']=True;raw['valid_rate_evidence']=(raw.rate_num>0)|((num(raw.infusionrate)>0)&raw.unit.str.contains('ml/hr',case=False))
        # Explicit rate is used; volumes and undocumented concentration are not converted.
    elif db=='amsterdam':
        a=read(R/'amsterdam_stays.csv',dtype=str);a['at']=num(a.admittedat)/1000;a['icu_possible']=a.location.str.contains('IC',regex=False);firstids=set(a[a.icu_possible].sort_values(['patientid','at','admissionid']).drop_duplicates('patientid').admissionid)
        s=pd.DataFrame({'sid':a.admissionid,'pid':a.patientid,'age':a.agegroup.map(lambda v:number(re.findall(r'\d+',v)[0]) if re.findall(r'\d+',v) else np.nan),'t0':a['at'],'t1':num(a.dischargedat)/1000,'first':a.admissionid.isin(firstids)&a.location.eq('IC'),'first_unknown':a.admissionid.isin(firstids)&~a.location.eq('IC'),'care_eligible':a.location.eq('IC'),'care_source':a.location,'hospital':'single','weight':a.weightgroup.map(midpoint),'height':a.heightgroup.map(midpoint)})
        wb=a.weightgroup.map(bounds);hb=a.heightgroup.map(bounds);s['bmi_lower']=[w[0]/(h[1]/100)**2 if h[1]>0 else np.nan for w,h in zip(wb,hb)];s['bmi_upper']=[w[1]/(h[0]/100)**2 if h[0]>0 else np.inf for w,h in zip(wb,hb)]
        rh=read(R/'amsterdam_rhythm.csv.gz',dtype=str);rh['isaf']=rh.value.str.contains(r'atr\s*fib|atrial fib',case=False,regex=True);afm=rh.groupby('admissionid').isaf.any();s['af']=s.sid.map(afm).fillna(False).astype(bool);s['af_source']=np.where(s.sid.isin(afm.index),'rhythm','not_observed')
        raw=read(R/'amsterdam_drugs.csv.gz',dtype=str);raw=raw[raw.itemid.isin(dmap)].copy();raw['sid']=raw.admissionid;raw['drug']=raw.itemid.map(dmap);raw['start']=num(raw.start)/1000;raw['end']=num(raw.stop)/1000;raw['rate_num']=num(raw.dose);raw['unit']=raw.doseunit+np.where(num(raw.doserateperkg).eq(1),'/kg/','/')+raw.doserateunit;raw['dose_weight']=raw.sid.map(s.set_index('sid').weight);raw['continuous']=num(raw.iscontinuous).eq(1)&raw.ordercategory.str.contains('Spuitpompen',case=False);raw['valid_rate_evidence']=raw.rate_num>0;raw['source_unit']=raw.unit;raw['rawid']=raw.orderid+'|'+raw.itemid+'|'+raw.start.astype(str)
    elif db=='sicdb':
        a=read(SIC/'cases.csv.gz',dtype=str);a['firstoffset']=num(a.OffsetAfterFirstAdmission)+num(a.ICUOffset);firstids=set(a.sort_values(['PatientID','firstoffset','CaseID']).drop_duplicates('PatientID').CaseID)
        refs=read(SIC/'d_references.csv.gz',dtype=str);units=refs.set_index('ReferenceGlobalID').ReferenceValue.to_dict();s=pd.DataFrame({'sid':a.CaseID,'pid':a.PatientID,'age':num(a.AgeOnAdmission),'t0':num(a.ICUOffset),'t1':num(a.TimeOfStay),'first':a.CaseID.isin(firstids),'first_unknown':a.firstoffset.isna(),'care_eligible':True,'care_source':a.HospitalUnit.map(units).fillna(a.HospitalUnit),'hospital':'single','weight':num(a.WeightOnAdmission),'height':num(a.HeightOnAdmission)})
        unitrefs=refs[refs.ReferenceName.eq('HospitalUnit')].copy();unitrefs['classification']=np.where(unitrefs.ReferenceDescription.str.contains('intensive care',case=False),'ICU',np.where(unitrefs.ReferenceDescription.str.contains('intermediate|mixed',case=False),'MC','unknown'));unitrefs[['ReferenceGlobalID','ReferenceValue','ReferenceDescription','classification']].to_csv(OUT/'map_careunit_sicdb.csv',index=False)
        ulog=read(SIC/'unitlog.csv.gz',dtype=str);umap=unitrefs.set_index('ReferenceGlobalID').classification.to_dict();ulog['class']=ulog.HospitalUnit.map(umap).fillna('unknown');ug=ulog.groupby('CaseID')['class'].agg(lambda v:set(v));pure=set(ug[ug.map(lambda v:v=={'ICU'})].index);possible=set(ug[ug.map(lambda v:'ICU' in v)].index)
        firsticu=set(a[a.CaseID.isin(possible)].sort_values(['PatientID','firstoffset','CaseID']).drop_duplicates('PatientID').CaseID);s['care_eligible']=s.sid.isin(pure);s['first']=s.sid.isin(firsticu)&s.care_eligible;s['first_unknown']=(s.sid.isin(firsticu)&~s.care_eligible)|a.firstoffset.isna();s['care_evidence']=s.sid.map(ug.map(lambda v:'+'.join(sorted(v)))).fillna('no_unitlog')
        # SICdb released admission mass is in grams; the documentation's 5 kg is rounding width.
        s['weight']=s.weight/1000.0
        code=a.ICD10Main.str.upper().str.replace('.', '', regex=False).str.strip()
        explicit_af=code.str.match(r'^I48[012](?:[0-9]*)$')
        text_af=a.ICD10MainText.str.contains(r'atrial fibr|vorhofflim',case=False,regex=True)
        text_flutter=a.ICD10MainText.str.contains(r'flutter|vorhofflatt',case=False,regex=True)
        s['af']=explicit_af|(text_af&~text_flutter&~code.str.match(r'^I48[349]'))
        s['af_source']=np.where(s.af,'main_diagnosis_confirmed_af','not_observed')
        raw=read(SIC/'medication.csv.gz',dtype=str);raw=raw[raw.DrugID.isin(dmap)].copy();raw['sid']=raw.CaseID;raw['drug']=raw.DrugID.map(dmap);raw['start']=num(raw.Offset);raw['end']=num(raw.OffsetDrugEnd);raw['rate_num']=num(raw.AmountPerMinute);um=dm.set_index('itemid').unitname.to_dict();raw['unit']=raw.DrugID.map(um)+'/min';raw['dose_weight']=raw.sid.map(s.set_index('sid').weight);raw['continuous']=num(raw.IsSingleDose).eq(0);raw['valid_rate_evidence']=raw.rate_num>0;raw['source_unit']=raw.unit;raw['rawid']=raw.id
    s['weight']=s.weight.where(s.weight>0);s['height']=s.height.where(s.height>0);s['bmi']=s.weight/(s.height/100)**2;s['adult']=s.age>=18;s['is_primary']=s.adult&s['first']&s.care_eligible&~s.first_unknown
    save(c,'stays',s)
    conv=[convert_rate(v,u,w,d) for v,u,w,d in zip(raw.rate_num,raw.unit,raw.dose_weight,raw.drug)];raw['dose']=[x[0] for x in conv];raw['dose_reason']=[x[1] for x in conv]
    if db=='eicu':
        # Collapse simultaneous aliases of the same drug, never sum duplicate channels.
        gr=raw.groupby(['sid','drug','start'],sort=False)
        d=gr.agg(dose=('dose','median'),dose_min=('dose','min'),dose_max=('dose','max'),ongoing=('valid_rate_evidence','max'),rawid=('rawid','first'),unknown=('dose',lambda x:x.isna().all())).reset_index()
        conflict=(d.dose_max-d.dose_min)>np.maximum(1e-6,d.dose_max.abs()*.01);d.loc[conflict,'dose']=np.nan
        d['dose_reason']=np.where(conflict,'conflicting_simultaneous_rates',np.where(d.dose.isna(),'unit_weight_or_rate_unknown','known'))
        d=d.sort_values(['sid','drug','start']);d['next']=d.groupby(['sid','drug']).start.shift(-1);d['end']=np.minimum(d['next'].fillna(np.inf),d.start+3600);d['end30']=np.minimum(d['next'].fillna(np.inf),d.start+1800);d=d[d.ongoing].copy()
    else:
        d=raw[raw.continuous&raw.valid_rate_evidence].copy();d['end30']=d['end']
    d=d.merge(s[['sid','t0','t1','adult']],on='sid');d['start']=np.maximum(d.start,d.t0);d['end']=np.minimum(d.end,d.t1);d['end30']=np.minimum(d.end30,d.t1);d=d[d.adult&(d.end>d.start)]
    cols=['sid','drug','start','end','end30','dose','dose_reason','rawid'];save(c,'drugs',d[cols])
    q=raw.groupby(['drug','source_unit','dose_reason'],dropna=False).agg(n_rows=('sid','size'),patients_or_stays=('sid','nunique')).reset_index();q.to_csv(OUT/f'phase0_dose_conversion_{db}.csv',index=False)
    raw.groupby('drug').agg(source_records=('sid','size'),source_stays=('sid','nunique')).to_csv(OUT/f'phase0_drug_usage_{db}.csv')
    exported(c,f'phase0_careunit_{db}.csv','SELECT care_source,count(*) n_stays,sum(is_primary::INT) primary_stays,sum(first_unknown::INT) first_unknown_stays FROM stays GROUP BY care_source')
    print(json.dumps({'db':db,'stage':'stays_drugs_ready','source_stays':len(s),'continuous_segments':len(d)}));c.close()
if __name__=='__main__':run_guard(sys.argv[1],build)
