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
c=duckdb.connect();c.execute("SET memory_limit='512MB';SET threads=2;")
refs=read(SIC/'d_references.csv.gz',dtype=str);cr=refs[refs.ReferenceName.eq('Laboratory')&refs.ReferenceValue.str.contains('creatin|kreatin',case=False)];save(c,'lab_dict',cr)
z=csvq(SIC/'laboratory.csv.gz');exported(c,'sicdb_creatinine_source_presence.csv',f'SELECT d.ReferenceValue AS concept_label,d.ReferenceUnit AS unitname,count(*) AS n_records,count(try_cast(x.Offset AS DOUBLE)) AS timed_records FROM {z} x JOIN lab_dict d ON x.LaboratoryID=d.ReferenceGlobalID GROUP BY d.ReferenceValue,d.ReferenceUnit')
refs[refs.apply(lambda r:bool(re.search(r'atrial|fibrill|vorhof|rhythm|rhytm|arrhyth|\bafib\b|\bvhf\b',str(r.to_dict()),re.I)),axis=1)].to_csv(OUT/'sicdb_rhythm_dictionary_search.csv',index=False)
for db,kind,field,value in [('mimic','bp','linksto','chartevents'),('mimic','vaso','linksto','inputevents'),('sicdb','bp','category','SignalFloat')]:
    p=OUT/f'map_{kind}_{db}.csv';m=read(p,dtype=str);m['count_status']=np.where(m[field].eq(value),'source_counted','not_applicable_to_target_table');m.loc[~m[field].eq(value),'rows']='';m.to_csv(p,index=False)
for db in ['mimic','eicu']:
    j=json.loads((R/f'audit_summary_{db}.json').read_text(encoding='utf-8'));pd.DataFrame([{'db':db,'audit_type':'positive_pairs','started_at':j['source_query_started_at'],'completed_at':j['completed_at'],'timing':j['timing'],'sample_size':j['sample_pairs'],'discrepancies':j['discrepancy_pairs']}]).to_csv(OUT/f'audit_execution_record_{db}.csv',index=False)
print('Final metadata checks complete.')
