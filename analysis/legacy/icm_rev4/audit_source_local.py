from __future__ import annotations

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
import sys as _cuff_sys
_cuff_sys.path.insert(0,str(_research_workspace/'sql/phase1'))
"""Independent source audit for 50 new low-MAP sampled MIMIC/eICU cuff-arterial pairs.

Does not import production cohort, drug, BP, or pairing modules. SQL files contain
paths and generic queries only. Every row-level artifact and detailed log stays
under phase1/_restricted. Source queries run only for the explicit `run` command.
"""


import argparse
import bisect
import collections
import csv
import datetime as dt
import hashlib
import json
import math
import re
import statistics
import subprocess
from pathlib import Path

ROOT = Path('@CUFF_WORKSPACE@')
SQLDIR = ROOT / 'sql' / 'icm_rev4'
OUT = ROOT / 'phase3' / 'icm_revision4' / 'analysis'
MAPS = ROOT / 'phase0'
PARAMFILE = ROOT / 'sql' / 'phase0_parameters.json'
import sys
sys.path.insert(0,str(ROOT/'sql/phase1'))
from audit_quality import P as QUALITY, anthro, mimic_anthro
RESTRICTED = OUT / '_restricted'
SAMPLE_COLUMNS = ['audit_id','sid','t','cuff_map','art_map','nee','nee_status',
                  'nearest_diff','nearest_art_t','art_values','episode_start','episode_stop']
EPOCH = dt.datetime(1970, 1, 1)


def read_csv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def number(x):
    try:
        v = float(str(x).strip())
        return v if math.isfinite(v) else None
    except (ValueError, TypeError):
        return None


def sid(x):
    return int(float(x))


def seconds(x):
    if not x:
        return None
    return (dt.datetime.fromisoformat(str(x)).replace(tzinfo=None) - EPOCH).total_seconds()


def median(values):
    v = [x for x in values if x is not None and math.isfinite(x)]
    return statistics.median(v) if v else None


def close(a, b, atol=1e-6, rtol=1e-6):
    if a is None or b is None:
        return a is None and b is None
    return abs(a-b) <= max(atol, rtol*max(abs(a), abs(b)))


def wsl_path(path):
    return Path(path).as_posix()


def sql_literal(value):
    return "'" + str(value).replace("'", "''") + "'"


def temp_import(table, path, expected=None):
    if Path(path).exists():
        with Path(path).open(encoding='utf-8-sig', newline='') as f:
            header = next(csv.reader(f))
    elif expected:
        header = expected
    else:
        raise FileNotFoundError('Required mapping file is unavailable')
    if not all(re.fullmatch(r'[A-Za-z_][A-Za-z_0-9]*', h) for h in header):
        raise ValueError('Unexpected CSV column name')
    cols = ','.join('"' + h + '" text' for h in header)
    return (f'CREATE TEMP TABLE {table} ({cols});\n'
            f'\\copy {table} FROM {sql_literal(wsl_path(path))} CSV HEADER\n')


def make_sql(db):
    sample = RESTRICTED / f'audit_sample_{db}.csv'
    text = '\\set ON_ERROR_STOP on\nSET client_min_messages=warning;\nSET TIME ZONE \'UTC\';\n'
    text += temp_import('audit_sample', sample, SAMPLE_COLUMNS)
    text += temp_import('audit_bpmap', MAPS / f'map_bp_{db}.csv')
    text += temp_import('audit_drugmap', MAPS / f'map_vaso_{db}.csv')
    if db == 'mimic':
        text += temp_import('audit_covmap', MAPS / 'map_covariates_mimic.csv')
    text += 'BEGIN READ ONLY;\nSET LOCAL statement_timeout=0;\n'
    queries = {}
    if db == 'mimic':
        queries['stays'] = '''SELECT i.stay_id,i.subject_id,i.hadm_id,i.intime,i.outtime
            FROM mimiciv_icu.icustays i JOIN (SELECT DISTINCT sid::bigint AS stay_id FROM audit_sample) a USING(stay_id)'''
        queries['bp'] = '''SELECT a.audit_id,c.stay_id,extract(epoch FROM c.charttime) AS source_t,
            c.itemid,m.classification AS modality,m.component,m.family,c.valuenum AS value,c.valueuom
            FROM audit_sample a JOIN mimiciv_icu.chartevents c ON c.stay_id=a.sid::bigint
            JOIN (SELECT DISTINCT itemid,classification,component,family FROM audit_bpmap WHERE selected='1') m
              ON c.itemid=m.itemid::bigint
            WHERE c.charttime >= timestamp '1970-01-01'+(a.t::double precision-300)*interval '1 second'
              AND c.charttime <= timestamp '1970-01-01'+(a.t::double precision+300)*interval '1 second'
              AND (m.classification='arterial' OR (m.classification='cuff' AND
                 c.charttime=timestamp '1970-01-01'+a.t::double precision*interval '1 second'))'''
        queries['drugs'] = '''SELECT x.stay_id,m.drug,x.itemid,x.starttime,x.endtime,x.rate,x.rateuom,
            x.patientweight,x.orderid,x.linkorderid,x.ordercategoryname,x.ordercategorydescription,x.statusdescription
            FROM mimiciv_icu.inputevents x JOIN (SELECT DISTINCT sid::bigint AS stay_id FROM audit_sample) a USING(stay_id)
            JOIN (SELECT DISTINCT itemid,drug FROM audit_drugmap WHERE selected='1') m ON x.itemid=m.itemid::bigint'''
        queries['weights'] = '''SELECT c.stay_id,extract(epoch FROM c.charttime) AS source_t,c.itemid,
            c.valuenum AS value,c.valueuom,m.label
            FROM mimiciv_icu.chartevents c JOIN (SELECT DISTINCT sid::bigint AS stay_id FROM audit_sample) a USING(stay_id)
            JOIN audit_covmap m ON c.itemid=m.itemid::bigint
            WHERE m.label IN ('Admission Weight (Kg)','Admission Weight (lbs.)','Daily Weight') OR lower(m.label) LIKE '%height%' '''
        queries['omr'] = '''SELECT o.subject_id,o.chartdate,o.seq_num,o.result_name,o.result_value
            FROM mimiciv_hosp.omr o JOIN (SELECT DISTINCT i.subject_id FROM mimiciv_icu.icustays i
                JOIN (SELECT DISTINCT sid::bigint AS stay_id FROM audit_sample) a USING(stay_id)) s USING(subject_id)
            WHERE o.result_name ~* 'weight' '''
    elif db == 'eicu':
        # Column names come from the reviewed eICU field mapping, not literal item IDs.
        bpmap = read_csv(MAPS / 'map_bp_eicu.csv')
        columns = {}
        for row in bpmap:
            if row.get('selected') == '1':
                table, col = row['itemid'].split('.', 1)
                if table not in {'vitalperiodic','vitalaperiodic'} or not re.fullmatch(r'[a-z]+',col):
                    raise ValueError('Unexpected eICU blood-pressure source mapping')
                columns[(row['classification'],row['component'])] = (table,col)
        bp_parts = []
        for modality in ['arterial','cuff']:
            source_table = columns[(modality,'map')][0]
            values = ','.join(f"({sql_literal(comp)},v.{columns[(modality,comp)][1]}::double precision)"
                              for comp in ['map','sbp','dbp'])
            time_filter = ('v.observationoffset::double precision*60 BETWEEN a.t::double precision-300 AND a.t::double precision+300'
                           if modality == 'arterial' else 'v.observationoffset::double precision*60=a.t::double precision')
            bp_parts.append(f'''SELECT a.audit_id,v.patientunitstayid AS stay_id,
                v.observationoffset::double precision*60 AS source_t,{sql_literal(modality)} AS modality,
                z.component,{sql_literal('art1' if modality == 'arterial' else 'cuff')} AS family,
                z.value,'mmHg' AS valueuom FROM audit_sample a
                JOIN public.{source_table} v ON v.patientunitstayid=a.sid::bigint
                CROSS JOIN LATERAL (VALUES {values}) z(component,value) WHERE {time_filter}''')
        queries['bp'] = ' UNION ALL '.join(bp_parts)
        queries['stays'] = '''SELECT p.patientunitstayid AS stay_id,p.unitdischargeoffset,p.admissionweight,p.admissionheight
            FROM public.patient p JOIN (SELECT DISTINCT sid::bigint AS patientunitstayid FROM audit_sample) a USING(patientunitstayid)'''
        queries['drugs'] = '''SELECT x.infusiondrugid,x.patientunitstayid AS stay_id,x.infusionoffset,
            x.drugname,x.drugrate,x.infusionrate,x.drugamount,x.volumeoffluid,m.drug
            FROM public.infusiondrug x JOIN (SELECT DISTINCT sid::bigint AS patientunitstayid FROM audit_sample) a USING(patientunitstayid)
            JOIN (SELECT DISTINCT itemid,drug FROM audit_drugmap WHERE selected='1') m ON x.drugname=m.itemid'''
    else:
        raise ValueError('Unknown audit database')
    for name, query in queries.items():
        # Only replace source-code line breaks; never split or collapse literal spaces.
        q = query.replace('\r', '').replace('\n', ' ')
        text += f"\\copy ({q}) TO {sql_literal(wsl_path(RESTRICTED / f'audit_source_{db}_{name}.csv'))} CSV HEADER\n"
    text += 'COMMIT;\n'
    path = SQLDIR / f'audit_source_{db}.sql'
    path.write_text(text, encoding='utf-8')
    return path


def normalize_unit(unit):
    u = str(unit or '').strip().lower().replace('µ','u').replace('μ','u')
    for old,new in [('micrograms','mcg'),('microgram','mcg'),('nanograms','ng'),('nanogram','ng'),
                    ('milligrams','mg'),('milligram','mg'),('minutes','min'),('minute','min'),
                    ('hours','hr'),('hour','hr'),('units','u'),('unit','u')]:
        u = u.replace(old,new)
    u = u.replace('ug','mcg').replace(' ','').replace('/h','/hr') if not u.endswith('/hr') else u.replace('ug','mcg').replace(' ','')
    return u


def source_unit_from_name(name):
    units = []
    for bracket in re.findall(r'\(([^()]*)\)', name):
        u = normalize_unit(bracket)
        if re.fullmatch(r'(?:mcg|mg|ng|g|u|ml)(?:/kg)?/(?:min|hr)',u):
            units.append(u)
    unique = set(units)
    return next(iter(unique)) if len(unique) == 1 else None


def nee_from_rate(drug, rate, unit, weight, params, name=''):
    if 'weight_quality_range' in params and (weight is None or not params['weight_quality_range'][0]<=weight<=params['weight_quality_range'][1]):
        weight=None
    if rate is None or rate <= 0:
        return None, 'not_positive_rate'
    unit = normalize_unit(unit)
    if unit in {'ml/hr','ml/min'}:
        # Only a stated concentration with explicit mass/activity unit is evidence.
        conc = re.search(r'(\d+(?:\.\d+)?)\s*(mcg|mg|ng|g|units?|iu|u)\b.*?(\d+(?:\.\d+)?)\s*ml\b',name,re.I)
        if not conc:
            return None,'volume_rate_without_verified_concentration'
        amount, amount_unit, volume = float(conc[1]),normalize_unit(conc[2]),float(conc[3])
        if volume <= 0:
            return None,'invalid_concentration'
        rate = rate*amount/volume
        unit = amount_unit+('/hr' if unit.endswith('/hr') else '/min')
    m = re.fullmatch(r'(mcg|mg|ng|g|u)(/kg)?/(min|hr)',unit)
    if not m:
        return None,'unknown_rate_unit'
    mass, perkg, period = m.groups()
    r = rate/(60 if period=='hr' else 1)
    if drug == 'vasopressin':
        if mass != 'u':
            return None,'vasopressin_without_activity_unit'
        if perkg:
            if weight is None or weight <= 0:
                return None,'weight_missing'
            r *= weight
        return r*params['vasopressin_units_min_factor'], 'known'
    if mass == 'u':
        return None,'activity_unit_for_mass_drug'
    r *= {'mcg':1.,'mg':1000.,'ng':.001,'g':1000000.}[mass]
    if not perkg:
        if weight is None or weight <= 0:
            return None,'weight_missing'
        r /= weight
    if drug == 'angiotensin_ii':
        return r*(1000*params['angiotensin_ii_ng_kg_min_factor']), 'known'
    factor = params['nee_factors_mcg_kg_min'].get(drug)
    return (r*factor,'known') if factor is not None else (None,'unknown_drug')


def weight_kg(value, unit_or_label):
    v = number(value)
    if v is None or v <= 0:
        return None
    unit = str(unit_or_label).lower()
    if 'lb' in unit or 'pound' in unit:
        return v*.45359237
    if 'kg' in unit or 'kilogram' in unit or 'daily weight' in unit:
        return v
    return None


def fallback_mimic_weights(stays, weight_rows, omr_rows):
    chart = collections.defaultdict(list)
    for x in weight_rows:
        s = sid(x['stay_id']); tm = number(x['source_t'])
        w = weight_kg(x['value'], x.get('valueuom','')+' '+x.get('label',''))
        if w is not None and tm is not None and tm >= seconds(stays[s]['intime']):
            chart[s].append((tm,w))
    omr = collections.defaultdict(list)
    for x in omr_rows:
        w = weight_kg(x['result_value'], x['result_name'])
        if w is not None:
            omr[sid(x['subject_id'])].append((x['chartdate'],number(x.get('seq_num')) or 0,w))
    result = {}
    for s,stay in stays.items():
        before = [x for x in omr[sid(stay['subject_id'])] if x[0] <= stay['intime'][:10]]
        if before:
            latest = max((x[0],x[1]) for x in before)
            result[s] = (median([x[2] for x in before if (x[0],x[1])==latest]),'omr_before_icu_date')
        elif chart[s]:
            earliest = min(x[0] for x in chart[s])
            result[s] = (median([x[1] for x in chart[s] if x[0]==earliest]),'first_post_icu_chart_weight')
        else:
            result[s] = (None,'weight_missing')
    return result


def mimic_intervals(stays, rows, fallback, params):
    records = collections.defaultdict(list)
    exclusion = collections.Counter()
    for x in rows:
        s=sid(x['stay_id']); rate=number(x.get('rate'))
        category=(x.get('ordercategoryname','')+' '+x.get('ordercategorydescription','')).lower()
        if x.get('ordercategorydescription')!='Continuous Med' or rate is None or rate<=0 or 'bolus' in category or x.get('statusdescription','').lower()=='rewritten':
            exclusion['not_eligible_direct_infusion'] += 1
            continue
        a=seconds(x['starttime']); b=seconds(x['endtime'])
        a=max(a,seconds(stays[s]['intime'])); b=min(b,seconds(stays[s]['outtime']))
        if b<=a:
            exclusion['nonpositive_or_outside_stay_interval'] += 1
            continue
        w=number(x.get('patientweight'))
        if w is None or w<=0:
            w,weight_source=fallback[s]
        else:
            weight_source='inputevents.patientweight'
        value,reason=nee_from_rate(x['drug'],rate,x.get('rateuom'),w,params)
        records[s].append({'drug':x['drug'],'start':a,'stop':b,'nee':value,'reason':reason,
                           'weight':w,'weight_source':weight_source,'source_key':x.get('orderid','')})
    return records,dict(exclusion)


def eicu_intervals(stays, rows, params):
    grouped=collections.defaultdict(list)
    exclusion=collections.Counter()
    for x in rows:
        if re.search(r'\bvolume\s*\(\s*ml\s*\)',x['drugname'],re.I):
            exclusion['cumulative_volume_channel']+=1
            continue
        s=sid(x['stay_id']); t=number(x['infusionoffset'])
        if t is None:
            exclusion['missing_source_time']+=1
            continue
        unit=source_unit_from_name(x['drugname']); rate=number(x.get('drugrate'))
        if unit=='ml/hr' and (rate is None or rate<=0):
            candidate=number(x.get('infusionrate'))
            if candidate is not None and candidate>0:
                rate=candidate
        weight=number(stays[s].get('admissionweight'))
        positive=rate is not None and rate>0
        value,reason=nee_from_rate(x['drug'],rate,unit,weight,params,x['drugname'])
        original_value,_=nee_from_rate(x['drug'],rate,unit,weight,
            {k:v for k,v in params.items() if k!='weight_quality_range'},x['drugname'])
        invalid_participating_weight=('weight_quality_range' in params and
            original_value is not None and value is None)
        grouped[(s,x['drug'],t*60)].append((positive,value,reason,weight,invalid_participating_weight))
    series=collections.defaultdict(list)
    for (s,drug,t),records in grouped.items():
        positive=any(r[0] for r in records)
        known=[r[1] for r in records if r[0] and r[1] is not None]
        if any(r[4] for r in records if r[0]):
            value=None; reason='invalid_weight_participated_in_original_alias_median'
        elif known and max(known)-min(known)>max(1e-6,max(abs(v) for v in known)*.01):
            value=None; reason='conflicting_simultaneous_alias_doses'
        elif known:
            value=statistics.median(known); reason='known'
        else:
            value=None; reason=';'.join(sorted({r[2] for r in records if r[0]})) or 'not_positive_rate'
        series[(s,drug)].append((t,positive,value,reason))
    intervals=collections.defaultdict(list)
    for (s,drug),events in series.items():
        events.sort(key=lambda x:x[0]); discharge=number(stays[s]['unitdischargeoffset'])*60
        for i,(t,positive,value,reason) in enumerate(events):
            if not positive:
                continue
            following=events[i+1][0] if i+1<len(events) else math.inf
            a=max(0.,t); b=min(following,t+60*params['gap_minutes_primary'],discharge)
            if b>a:
                intervals[s].append({'drug':drug,'start':a,'stop':b,'nee':value,'reason':reason,
                                     'weight':number(stays[s]['admissionweight']),
                                     'weight_source':'patient.admissionweight','source_key':''})
    return intervals,dict(exclusion)


def merge_intervals(intervals, gap):
    result=[]
    for a,b in sorted((x['start'],x['stop']) for x in intervals):
        if not result or a>result[-1][1]+gap:
            result.append([a,b])
        else:
            result[-1][1]=max(result[-1][1],b)
    return result


def dose_at(records, t, gap):
    bydrug=collections.defaultdict(list)
    for r in records:
        bydrug[r['drug']].append(r)
    details=[]; total=0.; unknown=False; expected=0
    for drug,rs in sorted(bydrug.items()):
        drug_episodes=merge_intervals(rs,gap)
        if not any(a<=t<b for a,b in drug_episodes):
            continue
        expected+=1
        active=[r for r in rs if r['start']<=t<r['stop']]
        if not active:
            unknown=True
            details.append({'drug':drug,'status':'gap_unknown','active_records':0})
        elif any(r['nee'] is None for r in active):
            unknown=True
            details.append({'drug':drug,'status':'dose_unknown','active_records':len(active),
                            'reasons':sorted({r['reason'] for r in active if r['nee'] is None})})
        else:
            amount=sum(r['nee'] for r in active); total+=amount
            details.append({'drug':drug,'status':'known','active_records':len(active),'nee':amount,
                            'weight_sources':sorted({r['weight_source'] for r in active})})
    if expected==0:
        return None,'no_active_component_in_bridged_global_gap',details
    return (None,'unknown',details) if unknown else (total,'known',details)


def clean_bp(rows, params):
    components=collections.defaultdict(lambda:collections.defaultdict(list))
    for r in rows:
        value=number(r['value']); t=number(r['source_t'])
        if value is not None and t is not None:
            components[(r['modality'],r['family'],t)][r['component']].append(value)
    valid=collections.defaultdict(list); rejected=[]
    for (modality,family,t),c in components.items():
        values={name:median(v) for name,v in c.items()}
        mp=values.get('map'); sbp=values.get('sbp'); dbp=values.get('dbp'); causes=[]
        if mp is None:
            causes.append('map_missing')
        elif not params['map_range_inclusive'][0]<=mp<=params['map_range_inclusive'][1]:
            causes.append('map_range')
        if sbp is not None and sbp>params['exclude_sbp_gt']:
            causes.append('sbp_gt_limit')
        if sbp is not None and dbp is not None:
            if sbp<=dbp:
                causes.append('sbp_le_dbp')
            if sbp-dbp<params['exclude_pulse_pressure_lt']:
                causes.append('pulse_pressure_lt_limit')
        if causes:
            rejected.append({'modality':modality,'family':family,'t':t,'causes':causes})
        else:
            valid[(modality,t)].append(mp)
    # At a given modality/time each source family contributes at most one value.
    return {(m,t):median(v) for (m,t),v in valid.items()},rejected


def write_results(db, rows, summary):
    detail=RESTRICTED/f'audit_comparison_{db}.csv'
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with detail.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
    summary_path=RESTRICTED/f'audit_summary_{db}.json'
    summary_path.write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    return summary


def file_digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_contract(db):
    paths=[RESTRICTED/f'audit_sample_{db}.csv',PARAMFILE,ROOT/'phase3/icm_revision4/analysis/icm_rev4_params.yaml',
           MAPS/f'map_bp_{db}.csv',MAPS/f'map_vaso_{db}.csv']
    if db=='mimic':
        paths.append(MAPS/'map_covariates_mimic.csv')
    return {str(p.relative_to(ROOT)):file_digest(p) for p in paths}


def compare(db):
    params=json.loads((PARAMFILE).read_text(encoding='utf-8'))
    sample_path=RESTRICTED/f'audit_sample_{db}.csv'
    manifest_path=RESTRICTED/f'audit_query_manifest_{db}.json'
    if not manifest_path.exists():
        raise RuntimeError('No completed independent source-query manifest')
    query_manifest=json.loads(manifest_path.read_text(encoding='utf-8'))
    if query_manifest.get('input_contract')!=source_contract(db):
        raise RuntimeError('Audit input, parameters or mappings changed since the source query')
    for rel,expected_hash in query_manifest.get('source_files',{}).items():
        if file_digest(ROOT/rel)!=expected_hash:
            raise RuntimeError('A source-audit export changed after source query completion')
    sample=read_csv(sample_path)
    if len({r['audit_id'] for r in sample})!=len(sample):
        raise ValueError('Audit IDs are not unique')
    if len(sample)!=50:
        raise ValueError('Expected exactly 50 presampled audit pairs')
    def source(name):
        return read_csv(RESTRICTED/f'audit_source_{db}_{name}.csv')
    stays={sid(x['stay_id']):x for x in source('stays')}
    if db=='mimic':
        fallback=fallback_mimic_weights(stays,source('weights'),source('omr'))
        all_intervals,exclusion=mimic_intervals(stays,source('drugs'),fallback,params)
    else:
        all_intervals,exclusion=eicu_intervals(stays,source('drugs'),params)
    clean_params={**params,'weight_quality_range':QUALITY['weight_kg_range']}
    if db=='mimic':
        clean_intervals,_=mimic_intervals(stays,source('drugs'),fallback,clean_params)
        source_anthro=mimic_anthro(stays,source('weights'))
    else:
        clean_intervals,_=eicu_intervals(stays,source('drugs'),clean_params)
        source_anthro={s:anthro(x.get('admissionweight'),x.get('admissionheight')) for s,x in stays.items()}
    bp=collections.defaultdict(list)
    for r in source('bp'):
        bp[r['audit_id']].append(r)
    results=[]; fields=collections.Counter(); reasons=collections.Counter(); unknown_reasons=collections.Counter()
    gap=params['gap_minutes_primary']*60
    for row in sample:
        s=sid(row['sid']); t=number(row['t']); audit_id=row['audit_id']; differences=[]
        records=all_intervals[s]
        episodes=merge_intervals(records,gap)
        containing=[(a,b) for a,b in episodes if a<=t<b]
        if len(containing)!=1:
            differences.append('therapy_episode_not_identified')
            ep_start=ep_stop=None
        else:
            ep_start,ep_stop=containing[0]
            if not close(ep_start,number(row['episode_start']),atol=1e-5,rtol=0):
                differences.append('episode_start')
            if not close(ep_stop,number(row['episode_stop']),atol=1e-5,rtol=0):
                differences.append('episode_stop')
        clean,rejections=clean_bp(bp[audit_id],params)
        cuff=clean.get(('cuff',t))
        if not close(cuff,number(row['cuff_map'])):
            differences.append('cuff_map_or_artifact')
        arterial=sorted((at,value) for (mod,at),value in clean.items() if mod=='arterial' and
                        abs(at-t)<=params['pair_window_primary']*60 and ep_start is not None and ep_start<=at<ep_stop)
        art=median([v for at,v in arterial])
        if not close(art,number(row['art_map'])):
            differences.append('art_median_or_artifact')
        nearest=min((abs(at-t) for at,v in arterial),default=None)
        ties=[at for at,v in arterial if close(abs(at-t),nearest,atol=1e-5,rtol=0)]
        if not close(nearest,number(row['nearest_diff']),atol=1e-5,rtol=0):
            differences.append('nearest_diff')
        original_nearest=number(row['nearest_art_t'])
        if original_nearest is None or not any(close(original_nearest,z,atol=1e-5,rtol=0) for z in ties):
            differences.append('nearest_art_t')
        if number(row['art_values'])!=len(arterial):
            differences.append('art_values')
        nee,status,components=dose_at(records,t,gap)
        clean_nee,clean_status,clean_components=dose_at(clean_intervals[s],t,gap)
        if clean_nee is not None and clean_nee>QUALITY['nee_max']:
            clean_nee=None;clean_status='quality_nee_above_5'
        av=source_anthro[s]
        for source_key,expected_key in [('source_weight','expected_phase0_weight'),('source_height','expected_phase0_height')]:
            if not close(av[source_key],number(row.get(expected_key))):
                differences.append(source_key)
        expected=number(row['nee'])
        if not close(nee,expected,atol=1e-6,rtol=1e-6):
            differences.append('nee')
        expected_known=str(row['nee_status']).strip().lower()=='known'
        if (status=='known')!=expected_known:
            differences.append('nee_known_status')
        for name in differences:
            fields[name]+=1
        if differences:
            reason=';'.join(differences);reasons[reason]+=1
        if status!='known':
            unknown_reasons[status]+=1
        results.append({**row,**av,'source_clean_nee':clean_nee,'source_clean_nee_status':clean_status,'clean_component_evidence':json.dumps(clean_components,ensure_ascii=False),'source_cuff_map':cuff,'source_art_map':art,'source_art_values':len(arterial),
                        'source_nearest_diff':nearest,'source_nearest_art_t_candidates':json.dumps(ties),
                        'source_nee':nee,'source_nee_status':status,'source_episode_start':ep_start,'source_episode_stop':ep_stop,
                        'component_evidence':json.dumps(components,ensure_ascii=False),
                        'artifact_evidence':json.dumps(rejections,ensure_ascii=False),
                        'discrepancy':int(bool(differences)),'discrepancy_fields':';'.join(differences)})
    summary={'database':db,'audit_kind':'independent_source_table_requery_and_separate_comparison',
             'sample_pairs':len(sample),'expected_sample_pairs':50,
             'discrepancy_pairs':sum(r['discrepancy'] for r in results),
             'discrepancy_fields':dict(fields),'discrepancy_reason_combinations':dict(reasons),
             'source_unknown_dose_reasons':dict(unknown_reasons),'source_record_exclusions':exclusion,
             'source_queries_executed':True,'scope':'source sample only; other four QA domains not performed by this audit',
             'timing':'ICM rev4主分析结果已知之后；新高剂量50对样本写定后；独立重查真实源数据库。',
             'source_query_started_at':query_manifest['started_at'],
             'source_query_completed_at':query_manifest['completed_at'],
             'completed_at':dt.datetime.now(dt.timezone.utc).isoformat(),
             'sample_sha256':hashlib.sha256(sample_path.read_bytes()).hexdigest(),
             'audit_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             'phase1_parameter_sha256':file_digest(ROOT/'phase3/icm_revision4/analysis/icm_rev4_params.yaml'),'source_clean_unknown_pairs':sum(x['source_clean_nee'] is None for x in results),'phase1_clean_output_comparison':'pending','parameter_sha256':hashlib.sha256((PARAMFILE).read_bytes()).hexdigest()}
    return write_results(db,results,summary)


def run(db):
    RESTRICTED.mkdir(parents=True,exist_ok=True)
    sample=RESTRICTED/f'audit_sample_{db}.csv'
    if not sample.exists():
        raise FileNotFoundError('Presampled restricted input is not yet available')
    sql=make_sql(db)
    contract=source_contract(db)
    started_at=dt.datetime.now(dt.timezone.utc).isoformat()
    database='mimiciv' if db=='mimic' else 'eicu'
    command=[_cuff_env.environ.get('CUFF_PSQL','psql'),'-X','-v',
             'ON_ERROR_STOP=1','-d',database,'-f',wsl_path(sql)]
    result=subprocess.run(command,capture_output=True,text=True,encoding='utf-8',errors='replace')
    (RESTRICTED/f'audit_psql_{db}.log').write_text(result.stdout+'\n'+result.stderr,encoding='utf-8')
    if result.returncode:
        raise RuntimeError('Source query failed; detailed diagnostics retained only in restricted log')
    if source_contract(db)!=contract:
        raise RuntimeError('Audit inputs or parameters changed during the source query')
    names=['stays','bp','drugs','weights','omr'] if db=='mimic' else ['stays','bp','drugs']
    files=[RESTRICTED/f'audit_source_{db}_{n}.csv' for n in names]
    manifest={'database':db,'started_at':started_at,
              'completed_at':dt.datetime.now(dt.timezone.utc).isoformat(),
              'input_contract':contract,'sql_sha256':file_digest(sql),
              'source_files':{str(p.relative_to(ROOT)):file_digest(p) for p in files}}
    (RESTRICTED/f'audit_query_manifest_{db}.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    return compare(db)


def selftest():
    p=json.loads((PARAMFILE).read_text(encoding='utf-8'))
    assert close(nee_from_rate('angiotensin_ii',20,'ng/kg/min',70,p)[0],.2)
    assert close(nee_from_rate('vasopressin',.04,'units/min',70,p)[0],.1)
    assert close(nee_from_rate('norepinephrine',7,'mcg/min',70,p)[0],.1)
    assert nee_from_rate('norepinephrine',3,'ml/hr',70,p)[0] is None
    assert source_unit_from_name('EPINEPHrine(Adrenalin)STD 4 mg 250 ml (mcg/min)')=='mcg/min'
    assert source_unit_from_name('Phenylephrine (mcg/kg/min) (mcg/kg/min)')=='mcg/kg/min'
    record=lambda drug,a,b,n:{'drug':drug,'start':a,'stop':b,'nee':n,'reason':'known','weight_source':'synthetic'}
    values=[record('norepinephrine',0,30,.1),record('norepinephrine',90,120,.2),record('epinephrine',0,120,.1)]
    assert dose_at(values,60,60)[0] is None
    assert close(dose_at(values,20,60)[0],.2)
    assert dose_at(values,120,60)[0] is None
    r=lambda c,v:{'modality':'arterial','family':'art1','source_t':'0','component':c,'value':str(v)}
    valid,rejected=clean_bp([r('map',65),r('sbp',80),r('dbp',75)],p)
    assert not valid and rejected
    valid,rejected=clean_bp([r('map',65)],p)
    assert valid[('arterial',0)]==65
    return {'synthetic_checks':'passed','source_queries_executed':False,'models_run':False}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('action',choices=['prepare','run','compare','selftest'])
    ap.add_argument('--db',choices=['mimic','eicu','both'],default='both')
    args=ap.parse_args()
    try:
        if args.action=='selftest':
            print(json.dumps(selftest(),ensure_ascii=False));return
        for db in (['mimic','eicu'] if args.db=='both' else [args.db]):
            if args.action=='prepare':
                path=make_sql(db)
                output={'database':db,'prepared_sql':str(path),'source_queries_executed':False}
            else:
                output=run(db) if args.action=='run' else compare(db)
            print(json.dumps(output,ensure_ascii=False))
    except Exception as error:
        # Exception text may contain a source value: keep it out of terminal/chat.
        if RESTRICTED.exists():
            import traceback
            (RESTRICTED/'audit_failure.log').write_text(traceback.format_exc(),encoding='utf-8')
        print(json.dumps({'status':'failed','error_type':type(error).__name__,
                          'diagnostics':'phase3/icm_revision4/analysis/_restricted/audit_failure.log','source_validation':'not_passed'},ensure_ascii=False))
        raise SystemExit(1)


if __name__=='__main__':
    main()
