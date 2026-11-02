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

"""Post hoc SICdb cluster-bootstrap diagnostic. Detailed outputs stay restricted."""
from pathlib import Path
import sys, csv, json, hashlib, datetime, platform
from types import SimpleNamespace

ROOT = Path('@CUFF_WORKSPACE@')
if '@' in str(ROOT) or not ROOT.is_absolute():
    raise RuntimeError('Prepare an absolute private workspace before executing this template.')
ROOT = ROOT.resolve()
for candidate in Path(__file__).resolve().parents:
    if (candidate/'analysis/prepare_workspace.py').is_file() and (candidate/'outputs/aggregate').is_dir():
        if ROOT == candidate or candidate in ROOT.parents:
            raise RuntimeError('Diagnostic outputs must remain outside the repository checkout.')
OUT = ROOT / 'restricted/rev8d_bootstrap'
DETAIL = OUT / 'replicate_detail'
for directory in [OUT, DETAIL, OUT/'data', OUT/'diagnostics']: directory.mkdir(parents=True,exist_ok=True)
import numpy as np
import pandas as pd
import scipy
from scipy.special import logit, expit
from scipy.stats import norm, skew, bootstrap
import duckdb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PARAM = json.loads((OUT / 'diagnostics/SICdb_bootstrap_parameters.json').read_text(encoding='utf-8'))
SEED, B = PARAM['seed'], PARAM['replicates']
LABELS = ['<=0.1', '>0.1-0.25', '>0.25-0.5', '>0.5-1.0', '>1.0']
checks = []

def write_csv(path, rows):
    # R1/R2: single-patient contribution and deletion details never leave
    # the restricted workspace. All other tabular exports pass the gate.
    from pathlib import Path
    import io, importlib.util
    path = Path(path)
    rows = list(rows)
    if path.name == 'SICdb_patient_contribution_shares.csv':
        rows = [r for r in rows if int(r.get('top_patients', 0)) >= 10]
    elif path.name == 'SICdb_high_contributor_exclusion.csv':
        path = DETAIL / path.name
    stream = io.StringIO(newline='')
    if rows:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    payload = stream.getvalue().encode('utf-8-sig')
    gate_path = next((parent / 'tools/small_cell_gate.py' for parent in Path(__file__).resolve().parents if (parent / 'tools/small_cell_gate.py').is_file()), ROOT / 'tools/small_cell_gate.py')
    if not gate_path.is_file():
        raise RuntimeError('data_gate_required')
    spec = importlib.util.spec_from_file_location('diagnostic_data_gate', gate_path)
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)
    hits, _ = gate.scan_bytes(payload, str(path))
    if hits:
        path = DETAIL / 'blocked_exports' / path.name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)

def record(item, observed, expected, scope='aggregate_reproduction', atol=1e-12):
    if isinstance(observed, (str, bool)) or isinstance(expected, (str, bool)):
        ok = observed == expected
    else:
        ok = bool(np.isclose(observed, expected, atol=atol, rtol=0, equal_nan=True))
    checks.append(dict(scope=scope, item=item, observed=observed, expected=expected, passed=bool(ok)))
    return ok

def read_csv(path):
    return list(csv.DictReader(path.open(encoding='utf-8-sig', newline='')))

def numeric(value):
    return float(value) if value not in ('', None) else np.nan

def ratio(num, den):
    return np.divide(num, den, out=np.full(np.shape(num), np.nan, dtype=float), where=np.asarray(den) > 0)

def bca_same_draws(ev, den, draws, estimate):
    """Cluster jackknife acceleration, midpoint tie rank for bootstrap bias."""
    jack = ratio(ev.sum() - ev, den.sum() - den)
    good = np.isfinite(draws)
    samples = draws[good]
    p_less = (np.sum(samples < estimate) + .5 * np.sum(samples == estimate)) / len(samples)
    z0 = norm.ppf(p_less)
    centered = np.mean(jack) - jack
    acc_den = 6 * np.sum(centered ** 2) ** 1.5
    acc = np.sum(centered ** 3) / acc_den if acc_den > 0 else np.nan
    za = norm.ppf([.025, .975])
    adjusted = norm.cdf(z0 + (z0 + za) / (1 - acc * (z0 + za)))
    ci = np.quantile(samples, adjusted, method='linear') if np.all(np.isfinite(adjusted)) and np.all(np.isfinite(jack)) else np.array([np.nan, np.nan])
    return ci, z0, acc, adjusted, jack

pairfile = ROOT / 'phase1/_restricted/pairs_sicdb.parquet'
con = duckdb.connect()
con.execute("SET threads=1")
df = con.execute('SELECT patient, cuff_map, art_map, nee, dose_stratum FROM read_parquet(?)', [str(pairfile)]).fetchdf()
pats = pd.Index(sorted(df.patient.unique()))
codes = pats.get_indexer(df.patient)
n = len(pats)
assert n > 1 and df[['patient', 'cuff_map', 'art_map']].notna().all().all()
assert ((df.art_map < 65) | (df.art_map >= 65)).all()
masks = {'all': np.ones(len(df), bool), 'known': df.nee.notna().to_numpy(), 'unknown': df.nee.isna().to_numpy()}
masks.update({label: df.dose_stratum.eq(label).to_numpy() for label in LABELS})
record('dose masks conserve all pairs', sum(masks[label].sum() for label in LABELS) + masks['unknown'].sum(), len(df), 'conservation')
record('known-dose strata conserve known pairs', sum(masks[label].sum() for label in LABELS), masks['known'].sum(), 'conservation')
specs, vectors = [], []
for endpoint in ['missed', 'false_positive']:
    eligible = (df.art_map < 65).to_numpy() if endpoint == 'missed' else (df.art_map >= 65).to_numpy()
    event = (df.cuff_map >= 65).to_numpy() if endpoint == 'missed' else (df.cuff_map < 65).to_numpy()
    for label, mask in masks.items():
        dm, em = eligible & mask, eligible & mask & event
        ev = np.bincount(codes, weights=em, minlength=n)
        den = np.bincount(codes, weights=dm, minlength=n)
        assert np.all(ev <= den)
        specs.append((endpoint, label, ev, den))
        vectors.extend([ev, den])
matrix = np.asarray(vectors).T
rng = np.random.default_rng(SEED)
weights = np.concatenate([rng.multinomial(n, np.full(n, 1/n), size=min(50, B-b)) for b in range(0, B, 50)])
assert weights.shape == (B, n) and np.all(weights.sum(axis=1) == n)
counts = weights @ matrix
num, denom = counts[:, 0::2], counts[:, 1::2]
br = ratio(num, denom)
assert np.all((br[np.isfinite(br)] >= 0) & (br[np.isfinite(br)] <= 1))
np.savez_compressed(DETAIL / 'sicdb_bootstrap_replicates.npz', weights=weights, numerator=num, denominator=denom, rates=br)

# Independent SQL recount of the existing paired file, not a new source audit.
sql = '''SELECT count(*) all_pairs, count(DISTINCT patient) all_patients,
 count(*) FILTER(WHERE art_map<65 AND cuff_map<65) tp,
 count(*) FILTER(WHERE art_map<65 AND cuff_map>=65) fn,
 count(*) FILTER(WHERE art_map>=65 AND cuff_map<65) fp,
 count(*) FILTER(WHERE art_map>=65 AND cuff_map>=65) tn,
 count(DISTINCT patient) FILTER(WHERE art_map>=65) fa_patients,
 count(DISTINCT patient) FILTER(WHERE art_map<65) missed_patients,
 count(*) FILTER(WHERE nee IS NOT NULL) known_pairs,
 count(*) FILTER(WHERE nee IS NULL) unknown_pairs
 FROM read_parquet(?)'''
cur = con.execute(sql, [str(pairfile)])
quad = dict(zip([c[0] for c in cur.description], cur.fetchone()))
record('four cells sum', sum(quad[k] for k in ['tp','fn','fp','tn']), quad['all_pairs'], 'conservation')
record('known plus unknown pairs', quad['known_pairs'] + quad['unknown_pairs'], quad['all_pairs'], 'conservation')

original = read_csv(ROOT / 'phase1/phase1_all_rates.csv')
publication = read_csv(ROOT / 'phase3/revision4/data/S2_missed_and_false_alarm.csv')
repro = []
for k, (ep, lab, ev, den) in enumerate(specs):
    src = next(r for r in original if r['database']=='sicdb' and r['endpoint']==ep and r['stratum']==lab and r['scenario']=='primary_5min')
    vals = br[np.isfinite(br[:,k]), k]
    ci = np.quantile(vals, [.025,.975]) if len(vals) else [np.nan,np.nan]
    boundary = (num[:,k]==0) | (num[:,k]==denom[:,k])
    transformed = logit(np.where(boundary, (num[:,k]+.5)/(denom[:,k]+1), br[:,k]))
    logit_se = np.std(transformed[np.isfinite(br[:,k])], ddof=1) if len(vals)>1 else np.nan
    est = float(ev.sum()/den.sum()) if den.sum() else np.nan
    values = {'numerator':int(ev.sum()),'denominator':int(den.sum()),'denominator_patients':int((den>0).sum()),'estimate':est,'ci_lower':ci[0],'ci_upper':ci[1],'bootstrap_valid':len(vals),'logit_se':logit_se}
    for key,value in values.items():record(f'{ep}/{lab}/{key}', value, numeric(src[key]))
    pubep = 'Missed hypotension rate' if ep=='missed' else 'False-alarm rate'
    pub = next(r for r in publication if r['database']=='SICdb' and r['endpoint']==pubep and r['stratum']==lab)
    for key in ['numerator','denominator','denominator_patients','estimate','ci_lower','ci_upper']:
        record(f'publication/{ep}/{lab}/{key}', values[key], numeric(pub[key]))
    repro.append(dict(database='SICdb',endpoint=ep,stratum=lab,**values,original_ci_lower=src['ci_lower'],original_ci_upper=src['ci_upper'],source_row_id=pub['row_id']))

summary, intervals, quantiles, concentration, deletion, histrows = [], [], [], [], [], []
private_clusters, private_jackknife = [], []
target_columns = {}
for ep in ['false_positive','missed']:
    k = next(i for i,s in enumerate(specs) if s[0]==ep and s[1]=='all')
    target_columns[ep]=k
    ev, den = specs[k][2:]
    estimate = ev.sum()/den.sum()
    samples = br[np.isfinite(br[:,k]), k]
    pct = np.quantile(samples, [.025,.975])
    bca_ci, z0, acceleration, adj, jack = bca_same_draws(ev, den, br[:,k], estimate)
    # SciPy checks the BCa calculation using these existing draws: zero new resamples.
    def statistic(event_counts, denominator_counts, axis=-1):
        return np.sum(event_counts, axis=axis)/np.sum(denominator_counts, axis=axis)
    independent = bootstrap((ev, den), statistic, paired=True, vectorized=True, method='BCa', n_resamples=0,
                            bootstrap_result=SimpleNamespace(bootstrap_distribution=samples), rng=np.random.default_rng(SEED))
    record(ep+'/BCa lower', bca_ci[0], independent.confidence_interval.low, 'BCa_independent_formula')
    record(ep+'/BCa upper', bca_ci[1], independent.confidence_interval.high, 'BCa_independent_formula')
    raw_boundary = (num[:,k]==0) | (num[:,k]==denom[:,k])
    m = np.where(raw_boundary, (num[:,k]+.5)/(denom[:,k]+1), br[:,k])
    sd_logit = np.std(logit(m[np.isfinite(br[:,k])]), ddof=1)
    wald = expit(logit(estimate) + np.array([-1,1])*norm.ppf(.975)*sd_logit)
    gaps = np.abs(bca_ci-pct)
    summary.append(dict(database='SICdb',endpoint=ep,stratum='all',numerator=int(ev.sum()),denominator=int(den.sum()),
        denominator_patients=int((den>0).sum()),resampling_patients=n,zero_denominator_patients=int((den==0).sum()),
        estimate=estimate,bootstrap_requested=B,bootstrap_valid=len(samples),zero_bootstrap_denominators=int((denom[:,k]==0).sum()),
        boundary_bootstrap_replicates=int(raw_boundary[np.isfinite(br[:,k])].sum()),bootstrap_mean=np.mean(samples),
        bootstrap_median=np.median(samples),bootstrap_sd=np.std(samples,ddof=1),bootstrap_mean_minus_estimate=np.mean(samples)-estimate,
        bootstrap_skewness=skew(samples,bias=False),logit_bootstrap_sd=sd_logit,bca_z0=z0,bca_acceleration=acceleration,
        bca_adjusted_lower_quantile=adj[0],bca_adjusted_upper_quantile=adj[1],bca_lower_difference_pp=(bca_ci[0]-pct[0])*100,
        bca_upper_difference_pp=(bca_ci[1]-pct[1])*100,bca_any_endpoint_difference_gt5pp=bool(np.any(gaps>.05))))
    for method, ci in [('Percentile',pct),('BCa',bca_ci),('Logit-Wald',wald)]:
        intervals.append(dict(database='SICdb',endpoint=ep,stratum='all',method=method,estimate=estimate,ci_lower=ci[0],ci_upper=ci[1],bootstrap_valid=len(samples),main_reported_interval=method=='Percentile'))
    for q in PARAM['quantiles']:
        quantiles.append(dict(database='SICdb',endpoint=ep,probability=q,estimate_quantile=np.quantile(samples,q),bootstrap_valid=len(samples)))
    # Ties follow original patient order without publishing any identifier or patient count.
    order = np.lexsort((np.arange(n), -ev))
    for top in [1,5,10]:
        chosen = order[:top]
        concentration.append(dict(database='SICdb',endpoint=ep,ranked_by='endpoint event count',top_patients=top,
            event_share=ev[chosen].sum()/ev.sum(),denominator_share=den[chosen].sum()/den.sum()))
    for top in [1,5]:
        chosen = order[:top]
        after = (ev.sum()-ev[chosen].sum())/(den.sum()-den[chosen].sum())
        deletion.append(dict(database='SICdb',endpoint=ep,excluded_top_event_contributors=top,original_estimate=estimate,
                             estimate_after_exclusion=after,difference_pp=(after-estimate)*100))
    for j in range(n):
        private_clusters.append(dict(endpoint=ep,patient=pats[j],events=int(ev[j]),eligible_pairs=int(den[j])))
        private_jackknife.append(dict(endpoint=ep,patient=pats[j],leave_one_out_estimate=jack[j]))
    hh, edges = np.histogram(samples,bins=np.linspace(0,1,41))
    for i,value in enumerate(hh):
        histrows.append(dict(database='SICdb',endpoint=ep,bin_lower=edges[i],bin_upper=edges[i+1],replicate_count=int(value)))
    record(ep+'/histogram counts', int(hh.sum()), B, 'diagnostic_conservation')

write_csv(DETAIL/'patient_contributions.csv',private_clusters)
write_csv(DETAIL/'patient_jackknife.csv',private_jackknife)
write_csv(OUT/'data/SICdb_bootstrap_reproduction.csv',repro)
write_csv(OUT/'data/SICdb_bootstrap_summary.csv',summary)
write_csv(OUT/'data/SICdb_bootstrap_intervals.csv',intervals)
write_csv(OUT/'data/SICdb_bootstrap_quantiles.csv',quantiles)
write_csv(OUT/'data/SICdb_patient_contribution_shares.csv',concentration)
write_csv(OUT/'data/SICdb_high_contributor_exclusion.csv',deletion)
write_csv(OUT/'data/SICdb_bootstrap_histogram_bins.csv',histrows)
write_csv(OUT/'diagnostics/SICdb_bootstrap_checks.csv',checks)

# Two diagnostic panels, using aggregate histogram bins only.
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,'axes.titlesize':11,
                     'svg.fonttype':'none','pdf.fonttype':42})
fig, axes = plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
for ax,ep,title,color in zip(axes,['false_positive','missed'],['A  False-alarm rate','B  Missed-hypotension rate'],['#3D6C91','#AE7048']):
    hh=[r for r in histrows if r['endpoint']==ep]
    ax.bar([r['bin_lower']*100 for r in hh],[r['replicate_count'] for r in hh],width=2.5,align='edge',color=color,alpha=.72,edgecolor='white',linewidth=.5)
    sr=next(r for r in summary if r['endpoint']==ep)
    ci=next(r for r in intervals if r['endpoint']==ep and r['method']=='Percentile')
    ax.axvline(sr['estimate']*100,color='#111111',linewidth=1.4,label='Original estimate')
    ax.axvline(ci['ci_lower']*100,color='#111111',linestyle='--',linewidth=1,label='Percentile 95% CI')
    ax.axvline(ci['ci_upper']*100,color='#111111',linestyle='--',linewidth=1)
    ax.set(xlabel='Pair-weighted rate (%)',ylabel='Bootstrap replicates',title=title,xlim=(0,75))
    ax.spines[['top','right']].set_visible(False)
    ax.text(.98,.96,f"Skewness: {sr['bootstrap_skewness']:.2f}\nValid replicates: {int(sr['bootstrap_valid']):,}",ha='right',va='top',transform=ax.transAxes,fontsize=9)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, frameon=False, fontsize=9, loc='outside lower center', ncol=2)
for ext in ['png','svg','pdf']:
    fig.savefig(OUT/'diagnostics'/f'SICdb_bootstrap_histogram.{ext}',dpi=300,facecolor='white')
plt.close(fig)

failed=[x for x in checks if not x['passed']]
receipt={'completed_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'parameters_recorded_at_utc':PARAM['recorded_at_utc'],
 'checks':len(checks),'failed':len(failed),'seed':SEED,'replicates':B,'resampling_patients':n,
 'rng':type(rng.bit_generator).__name__,'all_patient_sort':'lexicographic source patient string','batch_size':50,
 'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'duckdb':duckdb.__version__,'matplotlib':matplotlib.__version__},
 'source_sha256':{'paired_file':hashlib.sha256(pairfile.read_bytes()).hexdigest()},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'raw_source_requeried_this_round':False,'raw_source_audit_reused':'Existing study source-audit report; no new source query in this diagnostic',
 'main_result_files_modified':False,'bootstrap_detail_location':'restricted/rev8d_bootstrap','BCa_formula_checked_by_scipy_without_new_resampling':True}
(OUT/'diagnostics/SICdb_bootstrap_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
timeline_path=OUT/'diagnostics/SICdb_bootstrap_decision_timeline.json'
timeline=json.loads(timeline_path.read_text(encoding='utf-8')) if timeline_path.exists() else []
timeline.append({'recorded_at_utc':receipt['completed_at_utc'],'stage':'after diagnostic results','decision':'Production percentiles and logit SD verified; BCa and influence diagnostics retained separately without modifying manuscript results.','context':f"{len(checks)} checks, {len(failed)} discrepancies; author flags recorded in summary.",'numbering':'Root to assign D142 onward'})
(OUT/'diagnostics/SICdb_bootstrap_decision_timeline.json').write_text(json.dumps(timeline,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'status': 'completed', 'restricted_directory': str(DETAIL), 'check_rows': len(checks), 'patient_level_stdout': False}))
if failed:
    write_csv(DETAIL / 'failed_checks.csv', failed)
    print(json.dumps({'failed_count': len(failed), 'restricted_directory': str(DETAIL), 'clinical_values_returned': False}))
    raise RuntimeError('Aggregate diagnostic checks failed; inspect aggregate audit output.')
