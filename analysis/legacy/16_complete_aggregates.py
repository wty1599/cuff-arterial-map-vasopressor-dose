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
def finish(db):
    c=connect(db);f=OUT/f'phase0_feasibility_{db}.csv';old=pd.read_csv(f);rows=[]
    def add(section,metric,q,den=None):rows.append(dict(section=section,metric=metric,value=c.execute(q).fetchone()[0],unit='count',denominator=den))
    add('flow','adult_source_stays','SELECT count(*) FROM stays WHERE adult')
    add('flow','adult_care_eligible_stays','SELECT count(*) FROM stays WHERE adult AND care_eligible')
    add('flow','first_order_unknown_patients','SELECT count(DISTINCT pid) FROM stays WHERE adult AND first_unknown')
    for metric,q in [
      ('paired5_stays','SELECT count(DISTINCT sid) FROM pairs p JOIN stays s USING(sid) WHERE s.adult AND s.care_eligible'),
      ('paired5_patients','SELECT count(DISTINCT p.pid) FROM pairs p JOIN stays s USING(sid) WHERE s.adult AND s.care_eligible'),
      ('paired5_at_least10_stays','SELECT count(*) FROM (SELECT p.sid FROM pairs p JOIN stays s USING(sid) WHERE s.adult AND s.care_eligible GROUP BY p.sid HAVING count(*)>=10)')]:add('all_adult_stays',metric,q)
    add('dose_records','eligible_continuous_segments','SELECT count(*) FROM drugs')
    add('dose_records','unknown_nee_segments','SELECT count(*) FROM drugs WHERE dose IS NULL')
    add('dose_records','weight_missing_segments','SELECT count(*) FROM drugs WHERE dose_reason=\'weight_missing\'')
    # Direct arterial-line records when available; BP evidence is retained as a distinct proxy.
    if db=='mimic':
        z=csvq(R/'mimic_arterial_line.csv');add('prevalence','direct_arterial_line_records_patients',f'SELECT count(DISTINCT s.pid) FROM stays s JOIN {z} l ON l.stay_id=s.sid WHERE s.is_primary')
    elif db=='sicdb':
        refs=read(SIC/'d_references.csv.gz',dtype=str);line=refs[refs.ReferenceName.eq('Line')&refs.ReferenceValue.str.contains('arterial',case=False)];save(c,'line_dictionary',line);z=csvq(SIC/'data_range.csv.gz');add('prevalence','direct_arterial_line_records_patients',f'SELECT count(DISTINCT s.pid) FROM stays s JOIN {z} l ON l.CaseID=s.sid JOIN line_dictionary m ON l.DataID=m.ReferenceGlobalID WHERE s.is_primary')
    elif db=='amsterdam':
        # Midpoint BMI is undefined for open categories; bounds are a separate availability probe.
        a=read(R/'amsterdam_stays.csv',dtype=str);computed=[]
        for x in a.itertuples():
            w=bounds(x.weightgroup);h=bounds(x.heightgroup);valid=all(not np.isnan(v) for v in (*w,*h))
            low=w[0]/(h[1]/100)**2 if valid and h[1]>0 else np.nan;high=(w[1]/(h[0]/100)**2 if h[0]>0 else np.inf) if valid else np.nan
            computed.append((x.admissionid,low,high))
        save(c,'bmi_bounds_audit',pd.DataFrame(computed,columns=['sid','lower_bmi','upper_bmi']))
        add('interval_availability','bmi_definite_gt40_patients','SELECT count(DISTINCT pid) FROM pairs JOIN bmi_bounds_audit USING(sid) WHERE is_primary AND lower_bmi>40')
        add('interval_availability','bmi_possible_gt40_patients','SELECT count(DISTINCT pid) FROM pairs JOIN bmi_bounds_audit USING(sid) WHERE is_primary AND upper_bmi>40')
    keys={(r['section'],r['metric']) for r in rows};keep=old.loc[[tuple(x) not in keys for x in old[['section','metric']].itertuples(index=False,name=None)]];pd.concat([keep,pd.DataFrame(rows)],ignore_index=True).to_csv(f,index=False)
    exported(c,f'phase0_drug_primary_patients_{db}.csv','SELECT d.drug,count(DISTINCT s.pid) primary_patients,count(*) eligible_segments,count(*) FILTER(WHERE dose IS NULL) unknown_dose_segments FROM drugs d JOIN stays s USING(sid) WHERE is_primary GROUP BY d.drug')
    # Prove the sequential flow is monotone and quantify each exclusion.
    updated=pd.read_csv(f);vals=dict(zip(updated[updated.section.eq('flow')].metric,updated[updated.section.eq('flow')].value));sequence=['all_source_stays','adult_source_stays','adult_care_eligible_stays','adult_first_icu','vasopressor_patients','arterial_during_therapy','both_bp_during_therapy','paired5_patients','paired5_at_least10_patients'];flow=[]
    for left,right in zip(sequence,sequence[1:]):
        assert vals[left]>=vals[right],(left,right)
        flow.append({'from_stage':left,'to_stage':right,'before':int(vals[left]),'after':int(vals[right]),'excluded':int(vals[left]-vals[right]),'conserved':True})
    pd.DataFrame(flow).to_csv(OUT/f'phase0_flow_conservation_{db}.csv',index=False)
    c.close();print(json.dumps({'db':db,'aggregates_completed':True}))
if __name__=='__main__':run_guard(sys.argv[1],finish)
