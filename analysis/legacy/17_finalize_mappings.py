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
def update(db,kind,counts,source):
    p=OUT/f'map_{kind}_{db}.csv';m=read(p,dtype=str);m['rows']=m.itemid.map(counts).fillna(0).astype('int64');m['count_source']=source;m.to_csv(p,index=False)
mi=read(OUT/'mimic_bp_dictionary_source_counts.csv',dtype={'itemid':str});update('mimic','bp',mi.set_index('itemid').n_rows.to_dict(),'chartevents; non-chartevents items not applicable')
mi=read(R/'mimic_drugs.csv',usecols=['itemid'],dtype=str).itemid.value_counts().to_dict();update('mimic','vaso',mi,'inputevents; non-inputevents concepts not applicable')
for kind,fn in [('bp','amsterdam_numeric_item_census.csv'),('vaso','amsterdam_drug_item_census.csv')]:
    d=read(OUT/fn,dtype={'itemid':str});update('amsterdam',kind,d.groupby('itemid').rows.sum().to_dict(),fn)
d=read(OUT/'sicdb_numeric_item_census.csv',dtype={'itemid':str});update('sicdb','bp',d.set_index('itemid').rows.to_dict(),'data_float_h hourly blocks; not minute count')
d=read(SIC/'medication.csv.gz',usecols=['DrugID'],dtype=str).DrugID.value_counts().to_dict();update('sicdb','vaso',d,'medication source rows')
c=duckdb.connect();c.execute("SET memory_limit='512MB';SET threads=2;");counts={}
for mod,fn,table,cols in [('arterial','eicu_art.csv','vitalperiodic',['systemicmean','systemicsystolic','systemicdiastolic']),('cuff','eicu_cuff.csv','vitalaperiodic',['noninvasivemean','noninvasivesystolic','noninvasivediastolic'])]:
    q=csvq(R/fn);v=c.execute('SELECT '+','.join(f'count(try_cast({x} AS DOUBLE))' for x in cols)+' FROM '+q).fetchone();counts.update({table+'.'+k:n for k,n in zip(cols,v)})
update('eicu','bp',counts,'non-null source field records')
# Retain entire drug-name catalog counts, including explicitly excluded volume channels.
v=read(OUT/'eicu_drug_dictionary_all.csv',dtype=str);counts=v.set_index('drugname').rows.to_dict();update('eicu','vaso',counts,'complete infusiondrug drugname census')
for db in ['mimic','eicu']:
    j=json.loads((R/f'audit_summary_{db}.json').read_text(encoding='utf-8'));rows=[{'metric':'sampled_pairs','value':j['sample_pairs']},{'metric':'any_discrepancy_pairs','value':j['discrepancy_pairs']},{'metric':'source_queries_executed','value':int(j['source_queries_executed'])},{'metric':'unknown_nee_in_sample','value':sum(j['source_unknown_dose_reasons'].values())}];pd.DataFrame(rows).to_csv(OUT/f'phase0_source_audit_{db}.csv',index=False)
refs=read(SIC/'d_references.csv.gz',dtype=str);cr=refs[refs.ReferenceName.eq('Laboratory')&refs.ReferenceValue.str.contains('creatin|kreatin',case=False)];cr.to_csv(OUT/'sicdb_creatinine_dictionary.csv',index=False)
reg=cr.copy()
for col in reg:
    if isinstance(reg[col].dtype,pd.StringDtype):reg[col]=reg[col].astype(object)
c.register('lab_dict',reg)
z=csvq(SIC/'laboratory.csv.gz');r=c.execute(f'SELECT d.ReferenceValue AS concept_label,d.ReferenceUnit AS unitname,count(*) AS n_records,count(try_cast(x.Offset AS DOUBLE)) AS timed_records FROM {z} x JOIN lab_dict d ON x.LaboratoryID=d.ReferenceGlobalID GROUP BY d.ReferenceValue,d.ReferenceUnit').df();r.to_csv(OUT/'sicdb_creatinine_source_presence.csv',index=False)
refs[refs.apply(lambda r:bool(re.search('atrial|fibrill|vorhof|rhythm|rhytm|arrhyth|\\bafib\\b|\\bvhf\\b',str(r.to_dict()),re.I)),axis=1)].to_csv(OUT/'sicdb_rhythm_dictionary_search.csv',index=False)
print('Dictionary counts and aggregate audit summaries completed.')
