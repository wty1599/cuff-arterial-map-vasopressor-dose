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
ep=import_module('10_pair_and_counts');import time
c=connect('sicdb');ep.episodes(c,'decode_episodes',False,60,'end')
care=set(c.execute('SELECT sid FROM stays WHERE adult AND care_eligible').fetchnumpy()['sid']);intervals={}
for sid,start,stop in c.execute('SELECT sid,start,stop FROM decode_episodes').fetchall():
    if sid in care:intervals.setdefault(sid,[]).append((start,stop))
bm=mapping('sicdb','bp');bm=bm[bm.selected.eq('1')];mods=bm.set_index('itemid').classification.to_dict();comps=bm.set_index('itemid').component.to_dict();fams=bm.set_index('itemid').family.to_dict()
c.execute('CREATE OR REPLACE TABLE bp_long(sid VARCHAR,t DOUBLE,modality VARCHAR,component VARCHAR,family VARCHAR,val DOUBLE,itemid VARCHAR)')
scanned=0;considered=0;decoded=0;invalid=0;t0=time.monotonic();prefixes={}
for chunk in pd.read_csv(R/'sicdb_bp_hourly.csv.gz',dtype=str,keep_default_na=False,chunksize=25000):
    scanned+=len(chunk);sub=chunk[chunk.CaseID.isin(intervals)].copy()
    if sub.empty:continue
    sub['seconds']=pd.to_numeric(sub.Offset,errors='coerce');keep=[any(t<b and t+3600>a for a,b in intervals[s]) for s,t in zip(sub.CaseID,sub.seconds)];sub=sub.loc[keep]
    if sub.empty:continue
    for prefix,n in sub.rawdata.str[:2].value_counts().items():prefixes[prefix]=prefixes.get(prefix,0)+int(n)
    good=sub.rawdata.str.len().eq(482)&sub.rawdata.str[:2].isin(['0x','\\x']);invalid+=int((~good).sum());sub=sub[good].reset_index(drop=True);considered+=len(sub)
    if sub.empty:continue
    values=np.vstack([np.frombuffer(bytes.fromhex(v[2:]),dtype='<f4') for v in sub.rawdata]);mask=(values!=0)&np.isfinite(values);ri,mi=np.nonzero(mask);times=sub.seconds.to_numpy()[ri]+mi*60;sids=sub.CaseID.to_numpy()[ri]
    # Only treatment-time values are needed; all-stay raw availability remains in source_bp_hourly.
    in_treatment=np.zeros(len(ri),dtype=bool)
    for sid in set(sids):
        loc=np.flatnonzero(sids==sid);tt=times[loc];ok=np.zeros(len(loc),dtype=bool)
        for a,b in intervals[sid]:ok|=(tt>=a)&(tt<b)
        in_treatment[loc]=ok
    ri=ri[in_treatment];mi=mi[in_treatment];times=times[in_treatment];sids=sids[in_treatment]
    if not len(ri):continue
    ids=sub.DataID.to_numpy()[ri];frame=pd.DataFrame({'sid':sids,'t':times,'modality':[mods[v] for v in ids],'component':[comps[v] for v in ids],'family':[fams[v] for v in ids],'val':values[ri,mi].astype(float),'itemid':ids})
    save(c,'decode_batch',frame);c.execute('INSERT INTO bp_long SELECT * FROM decode_batch');decoded+=len(frame)
    (OUT/'_logs/sicdb_decode_progress.json').write_text(json.dumps({'hourly_rows_scanned':scanned,'decoded_treatment_values':decoded,'elapsed_seconds':round(time.monotonic()-t0),'prefixes':prefixes}))
pd.DataFrame([{'decoded_minute_values':decoded,'missing_or_invalid_raw_blocks':invalid,'blocks_considered':considered,'hourly_rows_scanned':scanned,'scope':'adult confirmed ICU stays, treatment periods; source availability separately full'}]).to_csv(OUT/'sicdb_raw_decode_qc.csv',index=False)
c.execute('DROP TABLE IF EXISTS decode_batch');c.close();print(json.dumps({'sicdb_decoder_complete':True,'decoded_values':decoded,'invalid_blocks':invalid,'prefixes':prefixes}))
