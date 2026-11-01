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
def episodes(c,name,per_drug,gap,endcol):
    part='sid,drug' if per_drug else 'sid';group=part+',grp'
    endcol='"'+endcol+'"'
    c.execute(f"""CREATE OR REPLACE TABLE {name} AS WITH a AS (
      SELECT sid,drug,rawid,start,{endcol} stop,max({endcol}) OVER(PARTITION BY {part} ORDER BY start,{endcol},rawid ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) prev
      FROM drugs WHERE {endcol}>start), b AS (
      SELECT *,sum(CASE WHEN prev IS NULL OR start>prev+{gap*60} THEN 1 ELSE 0 END) OVER(PARTITION BY {part} ORDER BY start,stop,rawid ROWS UNBOUNDED PRECEDING) grp FROM a)
      SELECT {part},min(start) AS "start",max(stop) AS "stop" FROM b GROUP BY {group}""")
def build(db):
    c=connect(db)
    for gap in [60,30]:
        endcol='end' if gap==60 else 'end30';episodes(c,f'episodes{gap}',False,gap,endcol);episodes(c,f'drug_episodes{gap}',True,gap,endcol)
        c.execute(f"CREATE OR REPLACE TABLE cuff{gap} AS SELECT row_number() OVER(ORDER BY b.sid,b.t) cid,b.sid,b.t,b.map_value cuff_map,b.incomplete_artifact_components,e.start episode_start,e.stop episode_stop FROM bp b JOIN episodes{gap} e ON b.sid=e.sid AND b.t>=e.start AND b.t<e.stop WHERE b.modality='cuff'")
    c.execute("CREATE OR REPLACE TABLE arterial AS SELECT b.*,e.start episode_start FROM bp b JOIN episodes60 e ON b.sid=e.sid AND b.t>=e.start AND b.t<e.stop WHERE b.modality='arterial'")
    c.execute("""CREATE OR REPLACE TABLE pairs_windows AS SELECT cf.cid,cf.sid,cf.t,cf.cuff_map,w.win,
    median(a.map_value) art_map,count(*) art_values,min(abs(a.t-cf.t)) nearest_diff,
    arg_min(a.t,struct_pack(delta_s := abs(a.t-cf.t),source_t := a.t)) nearest_art_t,
    avg(CASE WHEN abs(a.map_value/5-round(a.map_value/5))<0.0000001 THEN 1.0 ELSE 0 END) art_ends0or5_fraction,
    min(a.t) art_min_t,max(a.t) art_max_t
    FROM cuff60 cf CROSS JOIN (VALUES(2),(5),(10)) w(win) JOIN arterial a ON a.sid=cf.sid AND a.episode_start=cf.episode_start AND abs(a.t-cf.t)<=w.win*60
    GROUP BY cf.cid,cf.sid,cf.t,cf.cuff_map,w.win""")
    c.execute("CREATE OR REPLACE TABLE pair_base AS SELECT * FROM pairs_windows WHERE win=5")
    c.execute("""CREATE OR REPLACE TABLE pair_expected AS SELECT p.cid,count(DISTINCT e.drug) n_expected FROM pair_base p LEFT JOIN drug_episodes60 e ON p.sid=e.sid AND p.t>=e.start AND p.t<e.stop GROUP BY p.cid""")
    c.execute("""CREATE OR REPLACE TABLE pair_actual AS SELECT p.cid,count(DISTINCT d.drug) n_active,sum(d.dose) sum_dose,count(*) FILTER(WHERE d.drug IS NOT NULL AND d.dose IS NULL) n_unknown FROM pair_base p LEFT JOIN drugs d ON p.sid=d.sid AND p.t>=d.start AND p.t<d.end GROUP BY p.cid""")
    c.execute("""CREATE OR REPLACE TABLE pairs AS SELECT p.*,s.pid,s.is_primary,s.adult,s.hospital,s.bmi,s.af,s.af_source,CASE WHEN e.n_expected>0 AND e.n_expected=a.n_active AND a.n_unknown=0 THEN a.sum_dose END nee,
    CASE WHEN e.n_expected=0 OR e.n_expected<>a.n_active THEN 'undocumented_gap_component' WHEN a.n_unknown>0 THEN 'unit_weight_or_rate_unknown' ELSE 'known' END nee_status
    FROM pair_base p JOIN stays s USING(sid) JOIN pair_expected e USING(cid) JOIN pair_actual a USING(cid)""")
    records=[]
    def add(section,metric,value,unit='count',denom=None):records.append(dict(section=section,metric=metric,value=value,unit=unit,denominator=denom))
    def scalar(q):return c.execute(q).fetchone()[0]
    flow={
      'all_source_stays':'SELECT count(*) FROM stays',
      'adult_first_icu':'SELECT count(*) FROM stays WHERE is_primary',
      'adult_first_order_unknown':'SELECT count(*) FROM stays WHERE adult AND first_unknown',
      'vasopressor_patients':'SELECT count(*) FROM stays s WHERE is_primary AND EXISTS(SELECT 1 FROM drugs d WHERE d.sid=s.sid)',
      'arterial_during_therapy':'SELECT count(*) FROM stays s WHERE is_primary AND EXISTS(SELECT 1 FROM arterial a WHERE a.sid=s.sid)',
      'both_bp_during_therapy':'SELECT count(*) FROM stays s WHERE is_primary AND EXISTS(SELECT 1 FROM arterial a WHERE a.sid=s.sid) AND EXISTS(SELECT 1 FROM cuff60 cf WHERE cf.sid=s.sid)',
      'paired5_patients':'SELECT count(DISTINCT pid) FROM pairs WHERE is_primary',
      'paired5_at_least10_patients':'SELECT count(*) FROM (SELECT pid FROM pairs WHERE is_primary GROUP BY pid HAVING count(*)>=10)',
    }
    for k,q in flow.items():add('flow',k,scalar(q))
    for scope,cond in [('primary','s.is_primary'),('all_adult_stays','s.adult AND s.care_eligible')]:
        total=scalar(f'SELECT count(*) FROM cuff60 cf JOIN stays s USING(sid) WHERE {cond}');add(scope,'eligible_cuff_readings',total)
        for win in [2,5,10]:
            n=scalar(f'SELECT count(*) FROM pairs_windows p JOIN stays s USING(sid) WHERE {cond} AND win={win}');add(scope,f'paired{win}_readings',n,denom=total)
            add(scope,f'paired{win}_fraction',n/total if total else None,'fraction',total)
        n30=scalar(f'SELECT count(*) FROM cuff30 cf JOIN stays s USING(sid) WHERE {cond}');add(scope,'gap30_eligible_cuff_readings',n30)
        # Only count sensitivity, never agreement calculations.
        n30pairs=scalar(f"SELECT count(*) FROM cuff30 cf JOIN stays s USING(sid) WHERE {cond} AND EXISTS(SELECT 1 FROM bp a WHERE a.sid=cf.sid AND a.modality='arterial' AND a.t>=cf.episode_start AND a.t<cf.episode_stop AND abs(a.t-cf.t)<=300)");add(scope,'gap30_paired5_readings',n30pairs,denom=n30)
    strata=[('le0.1','nee<=0.1'),('gt0.1_le0.25','nee>0.1 AND nee<=0.25'),('gt0.25_le0.5','nee>0.25 AND nee<=0.5'),('gt0.5_le1','nee>0.5 AND nee<=1'),('gt1','nee>1'),('gt0.5_total','nee>0.5'),('unknown','nee IS NULL'),('bmi_gt40','bmi>40'),('af_any','af')]
    for label,cond in strata:
        p,n=c.execute(f'SELECT count(DISTINCT pid),count(*) FROM pairs WHERE is_primary AND ({cond})').fetchone();add('stratum_patients',label,p);add('stratum_pairs',label,n)
    totalpairs=scalar('SELECT count(*) FROM pairs WHERE is_primary');known=scalar('SELECT count(*) FROM pairs WHERE is_primary AND nee IS NOT NULL');add('ascertainability','nee_known_pairs_fraction',known/totalpairs if totalpairs else None,'fraction',totalpairs)
    # Missingness uses all treated primary patients, before conditioning on paired values.
    c.execute("CREATE OR REPLACE TABLE treated_primary AS SELECT * FROM stays s WHERE s.is_primary AND EXISTS(SELECT 1 FROM drugs d WHERE d.sid=s.sid)")
    nt=scalar('SELECT count(*) FROM treated_primary')
    for var in ['weight','height']:
        n=scalar(f'SELECT count(*) FROM treated_primary WHERE {var} IS NULL');add('missingness',var,n,denom=nt);add('missingness',var+'_fraction',n/nt if nt else None,'fraction',nt)
    for modality in ['arterial','cuff']:
        n=scalar(f"SELECT count(*) FROM treated_primary s WHERE NOT EXISTS(SELECT 1 FROM bp_quality b JOIN episodes60 e ON b.sid=e.sid AND b.t>=e.start AND b.t<e.stop WHERE b.sid=s.sid AND b.modality='{modality}' AND b.map_value IS NOT NULL)");add('missingness',modality+'_map_raw',n,denom=nt);add('missingness',modality+'_map_raw_fraction',n/nt if nt else None,'fraction',nt)
    for metric,cond in [('arterial_record_proxy','b.art'),('both_bp_records','b.art AND b.cuff')]:
        n=scalar(f'SELECT count(*) FROM stays s JOIN bp_presence_all b USING(sid) WHERE s.is_primary AND ({cond})');den=scalar('SELECT count(*) FROM stays WHERE is_primary');add('prevalence',metric,n,denom=den)
    # Recording-time and rounding summaries, using nearest source time for a median window.
    td=exported(c,f'phase0_timediff_{db}.csv',"""SELECT win,count(*) pairs,median(nearest_diff)/60 median_minutes,quantile_cont(nearest_diff,.25)/60 q1_minutes,quantile_cont(nearest_diff,.75)/60 q3_minutes,avg((nearest_diff=0)::INT) exact_zero_fraction,avg((nearest_diff>0 AND nearest_diff<=120)::INT) gt0_le2_fraction,avg((nearest_diff>120 AND nearest_diff<=300)::INT) gt2_le5_fraction,avg((nearest_diff>300 AND nearest_diff<=600)::INT) gt5_le10_fraction,avg((floor(t/60)=floor(nearest_art_t/60))::INT) same_clock_minute_fraction,avg(art_ends0or5_fraction) art_ending0or5_fraction FROM pairs_windows JOIN stays USING(sid) WHERE is_primary GROUP BY win ORDER BY win""")
    for n in ['raw_primary_art_values','raw_primary_cuff_values']:
        mod='arterial' if '_art_' in n else 'cuff';add('recording',n,scalar(f"SELECT count(*) FROM bp b JOIN stays s USING(sid) WHERE s.is_primary AND modality='{mod}'"))
    pd.DataFrame(records).to_csv(OUT/f'phase0_feasibility_{db}.csv',index=False)
    exported(c,f'phase0_nee_unknown_reasons_{db}.csv','SELECT nee_status,count(*) pairs,count(DISTINCT pid) patients FROM pairs WHERE is_primary GROUP BY nee_status')
    exported(c,f'phase0_af_sources_{db}.csv','SELECT af_source,count(DISTINCT pid) patients,count(*) pairs FROM pairs WHERE is_primary AND af GROUP BY af_source')
    if db=='eicu':
        exported(c,'phase0_eicu_hospital_pairs.csv','SELECT hospital,count(*) pairs,count(DISTINCT pid) patients,count(*)*1.0/sum(count(*)) OVER() fraction FROM pairs WHERE is_primary GROUP BY hospital ORDER BY pairs DESC')
        hc=csvq(R/'eicu_hospital_infusion_coverage.csv');exported(c,'phase0_eicu_hospital_coverage.csv',f"SELECT s.hospital,count(*) all_stays,count(*) FILTER(WHERE b.art) systemicmean_stays,count(*) FILTER(WHERE d.sid IS NOT NULL) vaso_stays,count(*) FILTER(WHERE b.art AND d.sid IS NOT NULL) both_vaso_and_systemicmean_stays,max(try_cast(h.infusion_any_stays AS BIGINT)) infusion_any_stays FROM stays s LEFT JOIN bp_presence_all b USING(sid) LEFT JOIN (SELECT DISTINCT sid FROM drugs) d USING(sid) LEFT JOIN {hc} h ON h.hospitalid=s.hospital GROUP BY s.hospital")
    print(json.dumps({'db':db,'stage':'counts_ready','primary_pairs':totalpairs,'primary_paired_patients':scalar('SELECT count(DISTINCT pid) FROM pairs WHERE is_primary')}));c.close()
if __name__=='__main__':run_guard(sys.argv[1],build)
