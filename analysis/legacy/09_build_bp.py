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
import struct,gzip
def build(db):
    c=connect(db);bm=mapping(db,'bp');bm=bm[bm.selected.eq('1')];save(c,'bp_map',bm)
    c.execute('CREATE OR REPLACE TABLE vaso_stays AS SELECT DISTINCT sid FROM drugs')
    resume='--resume-wide' in sys.argv
    if resume:
        assert db!='eicu' and c.execute('SELECT count(*) FROM bp_long').fetchone()[0]>0
    elif db=='mimic':
        q=csvq(R/'mimic_bp_cov.csv')
        c.execute(f"CREATE OR REPLACE TABLE bp_long AS SELECT x.stay_id sid,epoch(try_cast(x.charttime AS TIMESTAMP)) t,m.classification modality,m.component,m.family,try_cast(x.valuenum AS DOUBLE) val,x.itemid FROM {q} x JOIN bp_map m ON x.itemid=m.itemid")
        exported(c,'mimic_bp_item_counts.csv','SELECT itemid,count(*) n_rows FROM bp_long GROUP BY itemid')
    elif db=='eicu':
        art=csvq(R/'eicu_art.csv');cf=csvq(R/'eicu_cuff.csv')
        c.execute(f"CREATE OR REPLACE TABLE bp_wide_all AS SELECT patientunitstayid sid,try_cast(observationoffset AS DOUBLE)*60 t,'arterial' modality,try_cast(systemicmean AS DOUBLE) map_value,try_cast(systemicsystolic AS DOUBLE) sbp,try_cast(systemicdiastolic AS DOUBLE) dbp FROM {art} UNION ALL SELECT patientunitstayid sid,try_cast(observationoffset AS DOUBLE)*60 t,'cuff' modality,try_cast(noninvasivemean AS DOUBLE) map_value,try_cast(noninvasivesystolic AS DOUBLE) sbp,try_cast(noninvasivediastolic AS DOUBLE) dbp FROM {cf}")
        c.execute('CREATE OR REPLACE TABLE bp_presence_all AS SELECT sid,bool_or(modality=\'arterial\' AND map_value IS NOT NULL) art,bool_or(modality=\'cuff\' AND map_value IS NOT NULL) cuff FROM bp_wide_all GROUP BY sid')
        c.execute("CREATE OR REPLACE TABLE bp_wide AS SELECT sid,t,modality,median(map_value) map_value,median(sbp) sbp,median(dbp) dbp,count(*) raw_duplicates FROM bp_wide_all JOIN vaso_stays USING(sid) GROUP BY sid,t,modality")
        c.execute('DROP TABLE bp_wide_all')
    elif db=='amsterdam':
        q=csvq(R/'amsterdam_numeric_source.csv',encoding='latin-1')
        c.execute(f"CREATE OR REPLACE TABLE bp_long AS SELECT x.admissionid sid,try_cast(x.measuredat AS DOUBLE)/1000 t,m.classification modality,m.component,m.family,try_cast(x.value AS DOUBLE) val,x.itemid FROM {q} x JOIN bp_map m ON x.itemid=m.itemid")
        exported(c,'amsterdam_bp_item_counts.csv','SELECT itemid,count(*) n_rows FROM bp_long GROUP BY itemid')
    elif db=='sicdb':
        q=csvq(R/'sicdb_bp_hourly.csv.gz');c.execute(f"CREATE OR REPLACE TABLE source_bp_hourly AS SELECT * FROM {q}")
        c.execute("CREATE OR REPLACE TABLE bp_presence_all AS SELECT x.CaseID sid,bool_or(m.classification='arterial' AND try_cast(x.cnt AS INT)>0) art,bool_or(m.classification='cuff' AND try_cast(x.cnt AS INT)>0) cuff FROM source_bp_hourly x JOIN bp_map m ON x.DataID=m.itemid GROUP BY x.CaseID")
        c.execute('CREATE OR REPLACE TABLE bp_long(sid VARCHAR,t DOUBLE,modality VARCHAR,component VARCHAR,family VARCHAR,val DOUBLE,itemid VARCHAR)')
        # Bounded decoding; preserve positions of null slots rather than compressing time.
        reader=c.cursor();reader.execute('SELECT x.CaseID,x.DataID,x.Offset,x.rawdata,m.classification,m.component,m.family FROM source_bp_hourly x JOIN bp_map m ON x.DataID=m.itemid JOIN vaso_stays v ON x.CaseID=v.sid')
        def blocks():
            while True:
                part=reader.fetchmany(5000)
                if not part:break
                yield from part
        batch=[];decoded=0;invalid=0;considered=0
        for sid,item,t,raw,mod,component,family in blocks():
            considered+=1
            if not raw or len(raw)!=482 or not raw.startswith('0x'):invalid+=1;continue
            values=np.frombuffer(bytes.fromhex(raw[2:]),dtype='<f4');ind=np.flatnonzero((values!=0)&np.isfinite(values));times=float(t)+ind*60
            batch.extend(zip([sid]*len(ind),times,[mod]*len(ind),[component]*len(ind),[family]*len(ind),values[ind].astype(float),[item]*len(ind)))
            if len(batch)>=500000:
                df=pd.DataFrame(batch,columns=['sid','t','modality','component','family','val','itemid']);save(c,'_decoded_batch',df);c.execute('INSERT INTO bp_long SELECT * FROM _decoded_batch');decoded+=len(batch);batch=[]
        if batch:
            save(c,'_decoded_batch',pd.DataFrame(batch,columns=['sid','t','modality','component','family','val','itemid']));c.execute('INSERT INTO bp_long SELECT * FROM _decoded_batch');decoded+=len(batch)
        reader.close()
        pd.DataFrame([{'decoded_minute_values':decoded,'missing_or_invalid_raw_blocks':invalid,'blocks_considered':considered}]).to_csv(OUT/'sicdb_raw_decode_qc.csv',index=False)
        c.execute('DROP TABLE _decoded_batch')
    if db!='eicu':
        if db!='sicdb':c.execute("CREATE OR REPLACE TABLE bp_presence_all AS SELECT sid,bool_or(modality='arterial' AND component='map' AND val IS NOT NULL) art,bool_or(modality='cuff' AND component='map' AND val IS NOT NULL) cuff FROM bp_long GROUP BY sid")
        c.execute("CREATE OR REPLACE TABLE bp_wide AS SELECT sid,t,modality,family,median(val) FILTER(WHERE component='map') map_value,median(val) FILTER(WHERE component='sbp') sbp,median(val) FILTER(WHERE component='dbp') dbp,count(*) raw_component_rows FROM bp_long JOIN vaso_stays USING(sid) GROUP BY sid,t,modality,family")
    c.execute("CREATE OR REPLACE TABLE bp_quality AS SELECT *,map_value BETWEEN 20 AND 200 AND (sbp IS NULL OR sbp<=300) AND (sbp IS NULL OR dbp IS NULL OR (sbp>dbp AND sbp-dbp>=10)) AS valid, sbp IS NULL OR dbp IS NULL AS incomplete_artifact_components FROM bp_wide")
    c.execute("CREATE OR REPLACE TABLE bp AS SELECT sid,t,modality,median(map_value) map_value,bool_or(incomplete_artifact_components) incomplete_artifact_components FROM bp_quality WHERE valid GROUP BY sid,t,modality")
    exported(c,f'phase0_bp_quality_{db}.csv',"SELECT modality,count(*) candidate_time_family_rows,count(map_value) map_nonmissing,sum((map_value IS NULL)::INT) missing_map,sum((map_value IS NOT NULL AND NOT valid)::INT) observed_artifact_exclusions,sum((map_value IS NOT NULL AND incomplete_artifact_components)::INT) incomplete_artifact_components,sum(valid::INT) valid_rows FROM bp_quality GROUP BY modality")
    print(json.dumps({'db':db,'stage':'bp_ready','valid_timepoints':c.execute('SELECT count(*) FROM bp').fetchone()[0]}));c.close()
if __name__=='__main__':run_guard(sys.argv[1],build)
