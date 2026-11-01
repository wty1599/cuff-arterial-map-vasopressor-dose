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
"""New, bounded 30-pair audit directly against original SICdb source files.

No production extraction functions, models, or bootstrap code are imported.
All case-level records and source offsets remain in the restricted directory.
"""
from pathlib import Path
import sys, csv, json, math, struct, random, hashlib, datetime, traceback

import duckdb
import pandas as pd

ROOT = Path(r'@CUFF_WORKSPACE@')
OUT = ROOT / 'phase3/revision8a'
RESTRICTED = ROOT / 'restricted/rev8a_source_audit'
RAW = Path(r'@SICDB_ROOT@')
SEED = 20261027
RESTRICTED.mkdir(parents=True, exist_ok=True)

def register(con, name, frame):
    frame = frame.copy()
    for col in frame:
        if isinstance(frame[col].dtype, pd.StringDtype):
            frame[col] = frame[col].astype(object).where(frame[col].notna(), None)
    con.register('_incoming', frame)
    con.execute(f'CREATE OR REPLACE TABLE {name} AS SELECT * FROM _incoming')
    con.unregister('_incoming')

def source(path):
    return f"read_csv_auto('{path.as_posix()}', all_varchar=true)"

def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None

def merged_intervals(values):
    result = []
    for a, b in sorted(values):
        if result and a <= result[-1][1] + 3600:
            result[-1] = (result[-1][0], max(b, result[-1][1]))
        else:
            result.append((a, b))
    return result

def main():
    start = datetime.datetime.now(datetime.timezone.utc).isoformat()
    con = duckdb.connect(str(RESTRICTED / 'new_30_pair_audit.duckdb'))
    con.execute("SET memory_limit='1GB'")
    con.execute('SET threads=2')
    con.execute(f"SET temp_directory='{(RESTRICTED / 'tmp').as_posix()}'")
    pair_file = ROOT / 'phase1/_restricted/pairs_sicdb.parquet'
    frame = con.execute(f"SELECT cid,sid,pid,t,cuff_map,art_map,art_values,nearest_diff,nearest_art_t,art_min_t,art_max_t FROM read_parquet('{pair_file.as_posix()}') WHERE cuff_map<65 AND art_map>=65 ORDER BY cid").df()
    frozen=pd.read_csv(ROOT/'phase3/revision4/data/S2_missed_and_false_alarm.csv')
    expected=frozen.loc[(frozen.database=='SICdb') & (frozen.scenario=='Primary five-minute pairing') &
                        (frozen.endpoint=='False-alarm rate') & (frozen.stratum=='all'),'numerator']
    if len(expected)!=1:raise ValueError('Expected one frozen SICdb false-alarm aggregate')
    if len(frame)!=int(expected.iloc[0]):
        raise ValueError('False-alarm sampling frame differs from the frozen aggregate')
    sample = frame.iloc[sorted(random.Random(SEED).sample(range(len(frame)), 30))].copy()
    sample.insert(0, 'audit_id', [f'R8A-{i:02d}' for i in range(1, 31)])
    sample.to_csv(RESTRICTED / 'new_sample_30.csv', index=False)
    register(con, 'checks', sample)

    # Independently specify concept IDs after checking the original dictionary.
    concepts = pd.DataFrame([
        ('701','arterial','sbp','art1'), ('702','arterial','dbp','art1'),
        ('703','arterial','map','art1'), ('704','cuff','sbp','cuff'),
        ('705','cuff','dbp','cuff'), ('706','cuff','map','cuff')
    ], columns=['itemid','modality','component','family'])
    register(con, 'concepts', concepts)
    dictionary = con.execute(f"SELECT r.* FROM {source(RAW / 'd_references.csv.gz')} r JOIN concepts m ON r.ReferenceGlobalID=m.itemid").df()
    dictionary.to_csv(RESTRICTED / 'original_bp_dictionary.csv', index=False)
    if len(dictionary) != 6 or not dictionary.ReferenceUnit.eq('mmHg').all():
        raise ValueError('Blood pressure dictionary or units differ from expected concepts')
    query = f"""CREATE OR REPLACE TABLE original_hour_blocks AS
    SELECT DISTINCT x.id,x.CaseID sid,x.DataID,x.Offset,x.Val,x.cnt,x.rawdata,
           m.modality,m.component,m.family
    FROM {source(RAW / 'data_float_h.csv.gz')} x
    JOIN concepts m ON x.DataID=m.itemid
    JOIN checks z ON x.CaseID=z.sid
      AND TRY_CAST(x.Offset AS DOUBLE)<=z.t+300
      AND TRY_CAST(x.Offset AS DOUBLE)+3600>z.t-300"""
    (RESTRICTED / 'independent_source_queries.sql').write_text(query+';\n', encoding='utf-8')
    (RESTRICTED / 'sampling_manifest.json').write_text(json.dumps({
        'seed':SEED,'sampling_unit':'false-alarm pair','frame_pairs':len(frame),
        'sample_pairs':len(sample),'sample_patients':int(sample.pid.nunique()),
        'sampling_algorithm':'Python random.Random(seed).sample over cid-sorted pairs',
        'started_at_utc':start,'timing':'After original study results; before this new source audit',
        'sample_sha256':hashlib.sha256((RESTRICTED/'new_sample_30.csv').read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'source_files':{p.name:{'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns} for p in [RAW/'data_float_h.csv.gz', RAW/'cases.csv.gz', RAW/'medication.csv.gz',RAW/'d_references.csv.gz']}
    }, indent=2), encoding='utf-8')
    print(json.dumps({'stage':'resume_original_blocks_from_this_new_audit' if '--resume-blocks' in sys.argv else 'new_sample_fixed_original_gzip_scan_started','frame_pairs':len(frame),'sample_pairs':30,'sample_patients':int(sample.pid.nunique())}),flush=True)
    if '--resume-blocks' not in sys.argv:
        con.execute(query)
    blocks = con.execute('SELECT * FROM original_hour_blocks').df()
    decoded = []
    malformed = 0
    no_payload = 0
    for row in blocks.to_dict('records'):
        encoded = row['rawdata']
        if not isinstance(encoded,str) or encoded=='':
            no_payload += 1
            continue
        if not isinstance(encoded, str) or len(encoded)!=482 or encoded[:2] not in ('0x','\\x'):
            malformed += 1
            continue
        try:
            data = bytes.fromhex(encoded[2:])
        except ValueError:
            malformed += 1
            continue
        for minute in range(60):
            part = data[minute*4:(minute+1)*4]
            if part == b'\x00\x00\x00\x00':
                continue
            value = struct.unpack('<f',part)[0]
            if math.isfinite(value):
                decoded.append((row['id'],str(row['sid']),float(row['Offset'])+60*minute,row['modality'],row['component'],row['family'],value))
    register(con, 'decoded_minute_source',pd.DataFrame(decoded,columns=['source_record','sid','ts','modality','component','family','value']))
    con.execute("CREATE OR REPLACE TABLE minute_triplets AS SELECT sid,ts,modality,family,median(value) FILTER(WHERE component='map') pressure,median(value) FILTER(WHERE component='sbp') sbp,median(value) FILTER(WHERE component='dbp') dbp FROM decoded_minute_source GROUP BY sid,ts,modality,family")
    con.execute("CREATE OR REPLACE TABLE clean_minutes AS SELECT sid,ts,modality,median(pressure) pressure FROM minute_triplets WHERE pressure BETWEEN 20 AND 200 AND (sbp IS NULL OR sbp<=300) AND (sbp IS NULL OR dbp IS NULL OR (sbp>dbp AND sbp-dbp>=10)) GROUP BY sid,ts,modality")

    admissions = con.execute(f"SELECT CaseID sid,ICUOffset,TimeOfStay FROM {source(RAW/'cases.csv.gz')} WHERE CaseID IN (SELECT sid FROM checks)").df()
    meds = con.execute(f'''SELECT CaseID sid,DrugID,"Offset",OffsetDrugEnd,AmountPerMinute FROM {source(RAW/'medication.csv.gz')} WHERE CaseID IN (SELECT sid FROM checks) AND DrugID IN ('1502','1550','1562','1593','1618','1764') AND IsSingleDose='0' ''').df()
    admissions.to_csv(RESTRICTED/'original_case_bounds.csv',index=False)
    meds.to_csv(RESTRICTED/'original_medication_intervals.csv',index=False)
    bounds = {str(r.sid):(number(r.ICUOffset),number(r.TimeOfStay)) for r in admissions.itertuples()}
    intervals = {}
    for row in meds.to_dict('records'):
        sid=str(row['sid']); left=number(row['Offset']);right=number(row['OffsetDrugEnd']);dose=number(row['AmountPerMinute'])
        icu,end=bounds.get(sid,(None,None))
        if None in (left,right,dose,icu,end) or dose<=0:
            continue
        left=max(left,icu);right=min(right,end)
        if right>left:
            intervals.setdefault(sid,[]).append((left,right))
    exposure={sid:merged_intervals(v) for sid,v in intervals.items()}
    episodes=[]
    for row in sample.itertuples():
        interval=next((v for v in exposure.get(str(row.sid),[]) if v[0]<=row.t<v[1]),(None,None))
        episodes.append((row.audit_id,str(row.sid),row.t,*interval))
    register(con,'audit_episodes',pd.DataFrame(episodes,columns=['audit_id','sid','t','start','stop']))
    con.execute("CREATE OR REPLACE TABLE arterial_windows AS SELECT z.audit_id,b.ts,b.pressure,abs(b.ts-z.t) distance FROM audit_episodes z JOIN clean_minutes b ON b.sid=z.sid AND b.modality='arterial' AND b.ts BETWEEN z.t-300 AND z.t+300 AND b.ts>=z.start AND b.ts<z.stop")
    actual=con.execute("""SELECT z.audit_id,z.start source_episode_start,z.stop source_episode_stop,
    (SELECT median(pressure) FROM clean_minutes b WHERE b.sid=z.sid AND b.modality='cuff' AND abs(b.ts-z.t)<0.000001) source_cuff,
    (SELECT median(pressure) FROM arterial_windows b WHERE b.audit_id=z.audit_id) source_art,
    (SELECT min(distance) FROM arterial_windows b WHERE b.audit_id=z.audit_id) source_nearest_diff,
    (SELECT count(*) FROM arterial_windows b WHERE b.audit_id=z.audit_id) source_art_values,
    (SELECT min(ts) FROM arterial_windows b WHERE b.audit_id=z.audit_id) source_art_min_t,
    (SELECT max(ts) FROM arterial_windows b WHERE b.audit_id=z.audit_id) source_art_max_t,
    (SELECT median(pressure) FROM arterial_windows b WHERE b.audit_id=z.audit_id AND distance<0.000001) source_same_minute_art
    FROM audit_episodes z""").df()
    ties=con.execute("SELECT audit_id,to_json(list(ts) FILTER(WHERE distance=nearest)) nearest_ties FROM (SELECT *,min(distance) OVER(PARTITION BY audit_id) nearest FROM arterial_windows) GROUP BY audit_id").df()
    result=sample.merge(actual,on='audit_id').merge(ties,on='audit_id',how='left')
    matches=[]
    for expected,observed in [('cuff_map','source_cuff'),('art_map','source_art'),('nearest_diff','source_nearest_diff'),('art_values','source_art_values'),('art_min_t','source_art_min_t'),('art_max_t','source_art_max_t')]:
        tol=1e-4 if expected in ('nearest_diff','art_min_t','art_max_t') else 1e-6*(1+result[expected].abs())
        col=expected+'_match'; result[col]=(result[expected]-result[observed]).abs()<=tol;matches.append(col)
    result['nearest_art_t_match']=result.apply(lambda r:any(abs(float(r.nearest_art_t)-float(v))<1e-4 for v in json.loads(r.nearest_ties)) if isinstance(r.nearest_ties,str) else False,axis=1)
    matches.append('nearest_art_t_match')
    result['all_source_fields_match']=result[matches].all(axis=1)
    result['source_false_alarm']=(result.source_cuff<65)&(result.source_art>=65)
    result['same_minute_present']=result.source_same_minute_art.notna()
    result['same_minute_still_false_alarm']=result.same_minute_present&(result.source_cuff<65)&(result.source_same_minute_art>=65)
    result.to_csv(RESTRICTED/'source_comparison_30.csv',index=False)
    summary=[]
    def add(metric,value,denominator='',note=''):
        summary.append({'metric':metric,'value':value,'denominator':denominator,'note':note})
    add('sampling_frame_false_alarm_pairs',449)
    add('seed',SEED)
    add('sampled_pairs',len(result))
    add('sampled_unique_patients',int(result.pid.nunique()))
    add('source_hour_blocks',len(blocks),note='Original data_float_h.csv.gz; 60 little-endian float32 minute slots per hourly rawdata block')
    add('malformed_source_blocks',malformed,len(blocks))
    add('source_blocks_without_minute_payload',no_payload,len(blocks),'Source rawdata was null; all five such blocks had cnt=1; hourly summary values were not substituted for minute observations')
    add('source_blocks_decoded',len(blocks)-malformed-no_payload,len(blocks))
    add('bp_dictionary_concepts_mmHg',len(dictionary),6)
    add('source_episode_reconstructed_pairs',int(result.source_episode_start.notna().sum()),30)
    add('all_source_fields_match_pairs',int(result.all_source_fields_match.sum()),30)
    for col in matches:
        add(col+'_pairs',int(result[col].sum()),30)
    add('false_alarm_reproduced_pairs',int(result.source_false_alarm.sum()),30)
    add('same_minute_arterial_available_pairs',int(result.same_minute_present.sum()),30)
    add('same_minute_false_alarm_persists_pairs',int(result.same_minute_still_false_alarm.sum()),int(result.same_minute_present.sum()),'Descriptive audit tally in the 30 sampled pairs; not a new population sensitivity analysis')
    add('same_minute_false_alarm_resolves_pairs',int((result.same_minute_present&~result.same_minute_still_false_alarm).sum()),int(result.same_minute_present.sum()))
    add('nearest_time_difference_seconds_median',float(result.source_nearest_diff.median()))
    add('nearest_time_difference_seconds_min',float(result.source_nearest_diff.min()))
    add('nearest_time_difference_seconds_max',float(result.source_nearest_diff.max()))
    add('arterial_readings_per_pair_median',float(result.source_art_values.median()))
    add('arterial_readings_per_pair_min',int(result.source_art_values.min()))
    add('arterial_readings_per_pair_max',int(result.source_art_values.max()))
    add('source_cuff_values_match_and_unit_mmHg',int(result.cuff_map_match.sum()),30)
    add('source_arterial_values_match_and_unit_mmHg',int(result.art_map_match.sum()),30)
    add('waveform_quality_or_device_cause_verified',0,30,'This source audit cannot adjudicate waveform quality, damping, cuff technique, or why paired measurements were obtained')
    pd.DataFrame(summary).to_csv(OUT/'data/SICdb_false_alarm_source_audit_summary.csv',index=False)
    receipt={'completed_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sample_pairs':30,'source_field_discrepancy_pairs':int((~result.all_source_fields_match).sum()),'false_alarm_reproduced_pairs':int(result.source_false_alarm.sum()),'same_minute_available':int(result.same_minute_present.sum()),'same_minute_false_alarm_persists':int(result.same_minute_still_false_alarm.sum()),'malformed_source_blocks':malformed,'source_blocks_without_minute_payload':no_payload,'models_or_bootstrap_rerun':False,'patient_records_deliverable':False,'execution_note':'One SQL reserved-word parsing error was fixed; raw hourly blocks from this new audit were retained and independently decoded after repair'}
    (OUT/'qa/source_audit_receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    print(json.dumps(receipt),flush=True)
    con.close()

if __name__=='__main__':
    try:
        main()
    except Exception as exc:
        (RESTRICTED/'error.txt').write_text(traceback.format_exc(),encoding='utf-8')
        print(json.dumps({'status':'failed','error_class':type(exc).__name__,'details':'Restricted error log; no case-level records printed'}),flush=True)
        sys.exit(1)
