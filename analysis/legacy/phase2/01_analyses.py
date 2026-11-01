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

from common import *
from scipy.special import logit
import hashlib,datetime

def ratio_engine(d,db,scenario,only=None,limited=False):
    pats=pd.Index(sorted(d.patient.unique()));codes=pats.get_indexer(d.patient);N=len(pats)
    ms=masks(d,groups=(db in MAIN and scenario=='primary_5min'))
    if only:ms={k:v for k,v in ms.items() if k in only}
    art=d.art_map.to_numpy();cuff=d.cuff_map.to_numpy();lo=art<65;normal=cuff>=65;miss=lo&normal
    specs=[];vectors=[];quant=[]
    def register(domain,metric,lab,cat,eligible,event,cell_n=None):
        dc=np.bincount(codes,weights=eligible,minlength=N);ec=np.bincount(codes,weights=eligible&event,minlength=N)
        de=int(eligible.sum());ev=int((eligible&event).sum())
        specs.append(dict(database=db,scenario=scenario,domain=domain,metric=metric,stratum=lab,category=cat,numerator=ev,denominator=de,denominator_patients=int((dc>0).sum()),cell_pairs=de if cell_n is None else int(cell_n)))
        vectors.extend([ec,dc])
    for lab,mask in ms.items():
        if limited:
            register('A','deep_miss_art_lt55',lab,'all',mask&(art<55),normal)
            if db in MAIN:register('B','art_lt65_given_cuff_ge65',lab,'all_normal_cuff',mask&normal,lo)
            continue
        if not lab.startswith('dosegroup:'):
            z=d[mask&miss];delta=z.cuff_map-z.art_map
            quant.append(dict(database=db,scenario=scenario,stratum=lab,missed_pairs=len(z),missed_patients=z.patient.nunique(),difference_median=delta.median(),difference_q1=delta.quantile(.25),difference_q3=delta.quantile(.75)))
            for cat,bm in [('60-64',(art>=60)&(art<65)),('55-59',(art>=55)&(art<60)),('50-54',(art>=50)&(art<55)),('<50',art<50)]:
                register('A_distribution','arterial_bin_given_missed',lab,cat,mask&miss,bm,cell_n=(mask&miss&bm).sum())
            for threshold in [55,60]:register('A',f'deep_miss_art_lt{threshold}',lab,'all',mask&(art<threshold),normal)
        if db in MAIN:
            register('B','arterial_lt65_prevalence',lab,'all_pairs',mask,lo)
            register('B','art_lt65_given_cuff_ge65',lab,'all_normal_cuff',mask&normal,lo)
            for cat,cm in [('65-69',(cuff>=65)&(cuff<70)),('70-74',(cuff>=70)&(cuff<75)),('75-79',(cuff>=75)&(cuff<80)),('>=80',cuff>=80)]:
                register('B','art_lt65_given_cuff_bin',lab,cat,mask&cm,lo)
            if not lab.startswith('dosegroup:'):
                register('B','art_lt60_given_cuff65_69',lab,'65-69',mask&(cuff>=65)&(cuff<70),art<60)
                register('B','art_lt65_given_cuff_lt65',lab,'low_cuff',mask&(cuff<65),lo)
                for cat,bm in [('60-64',(art>=60)&(art<65)),('55-59',(art>=55)&(art<60)),('<55',art<55)]:
                    register('C_cross','cuff_ge65_given_art_bin',lab,cat,mask&bm,normal)
                    register('C_bin_fraction','art_bin_given_art_lt65',lab,cat,mask&lo,bm)
    mat=np.array(vectors).T;chunks=[];rng=np.random.default_rng(P['seed'])
    for b in range(0,2000,50):
        w=rng.multinomial(N,np.full(N,1/N),size=min(50,2000-b));chunks.append(w@mat)
    count=np.concatenate(chunks);ne=count[:,0::2];de=count[:,1::2]
    br=np.divide(ne,de,out=np.full_like(ne,np.nan),where=de>0)
    for i,row in enumerate(specs):
        good=np.isfinite(br[:,i]);ll,uu=np.quantile(br[good,i],[.025,.975]) if good.any() else (np.nan,np.nan)
        n=row['numerator'];den=row['denominator'];p=n/den if den else np.nan
        degenerate=bool(den and n in [0,den] and ll==uu)
        row.update(estimate=p,ci_lower=ll,ci_upper=uu,bootstrap_valid=int(good.sum()),small_cell=row['cell_pairs']<30,wide_ci=bool(uu-ll>.30),boundary=degenerate,precision='boundary_not_informative' if degenerate else 'insufficient' if row['cell_pairs']<30 or uu-ll>.30 else 'adequate' if den else 'not_estimable',ci_status='boundary_degenerate' if degenerate else 'ok' if den else 'no_denominator')
        pp=np.where((ne[:,i]==0)|(ne[:,i]==de[:,i]),(ne[:,i]+.5)/(de[:,i]+1),br[:,i]);logs=logit(pp[good]);row['logit_estimate']=float(logit((n+.5)/(den+1) if n in [0,den] else p)) if den else np.nan;row['logit_se']=float(np.std(logs,ddof=1)) if len(logs)>1 else np.nan
    out=pd.DataFrame(specs)
    if db in MAIN and limited:
        for lab,mask in ms.items():
            ix=out.stratum.eq(lab);den=int(mask.sum());num=int((mask&lo).sum())
            out.loc[ix,'layer_low_arterial_pairs']=num;out.loc[ix,'layer_pairs']=den;out.loc[ix,'layer_low_arterial_prevalence']=num/den if den else np.nan
        out['condition_low_arterial_pairs']=np.where(out.domain.eq('B'),out.numerator,np.nan)
        out['condition_low_arterial_prevalence']=np.where(out.domain.eq('B'),out.estimate,np.nan)
    if db in MAIN and not limited:
        prev=out[out.metric.eq('arterial_lt65_prevalence')].set_index('stratum')
        for col,new in [('numerator','layer_low_arterial_pairs'),('denominator','layer_pairs'),('estimate','layer_low_arterial_prevalence'),('ci_lower','layer_prevalence_ci_lower'),('ci_upper','layer_prevalence_ci_upper')]:out[new]=out.stratum.map(prev[col])
        # For the <60 probability, also show the <65 prevalence within that cuff denominator.
        out['condition_low_arterial_pairs']=np.nan;out['condition_low_arterial_prevalence']=np.nan
        bm=out.domain.eq('B')
        for i,row in out[bm].iterrows():
            if row.metric=='art_lt60_given_cuff65_69':
                other=out[out.stratum.eq(row.stratum)&out.metric.eq('art_lt65_given_cuff_bin')&out.category.eq('65-69')].iloc[0]
                out.loc[i,'condition_low_arterial_pairs']=other.numerator;out.loc[i,'condition_low_arterial_prevalence']=other.estimate
            else:
                out.loc[i,'condition_low_arterial_pairs']=row.numerator;out.loc[i,'condition_low_arterial_prevalence']=row.estimate
    out['bootstrap_column']=np.arange(len(out))
    out.to_csv(RES/f'ratio_specs_{db}_{scenario}.csv',index=False)
    np.savez_compressed(RES/f'bootstrap_{db}_{scenario}.npz',numerator=ne,denominator=de)
    q=pd.DataFrame(quant)
    return out,q

def severity(d,db,scenario='primary_5min',only=None):
    rows=[]
    for lab,mask in masks(d).items():
        if only and lab not in only:continue
        z=d[mask&d.art_map.lt(65)];a=z.art_map
        rows.append(dict(database=db,scenario=scenario,stratum=lab,low_arterial_pairs=len(z),patients=z.patient.nunique(),arterial_map_median=a.median(),arterial_map_q1=a.quantile(.25),arterial_map_q3=a.quantile(.75),art_lt55_pairs=int(a.lt(55).sum()),art_lt55_fraction=a.lt(55).mean()))
    return pd.DataFrame(rows)

def main():
    manifests=[];results={};distributions=[];deep=[];B=[];C=[];sever=[];sensitivity=[]
    for db in DBS:
        p=INP/f'pairs_{db}.parquet';manifests.append(dict(database=db,file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
        d=read_pairs(db);out,q=ratio_engine(d,db,'primary_5min');results[db]=out
        distributions.append(out[out.domain.eq('A_distribution')].merge(q,on=['database','scenario','stratum'],validate='many_to_one'))
        deep.append(out[out.domain.eq('A')])
        if db in MAIN:
            B.append(out[out.domain.eq('B')]);C.append(out[out.domain.eq('C_cross')]);sever.append(severity(d,db))
        sm=pd.read_parquet(INP/f'same_minute_{db}.parquet');s=d.drop(columns=['art_map','art_values','nearest_diff']).merge(sm,on='cid',validate='one_to_one')
        so,_=ratio_engine(s,db,'same_minute',['all','>1.0'],True);sensitivity.append(so)
        print(json.dumps(dict(database=db,phase2_ratios='complete',pairs=len(d))),flush=True)
    d=read_pairs('eicu');r=pd.read_csv(INP/'eicu_recovered_pairs.csv',usecols=['cid','recovered_nee']);d=d.merge(r,on='cid',validate='one_to_one')
    restore=d.nee_raw.isna()&d.recovered_nee.notna()&~d.dose_weight_flag&d.recovered_nee.between(0,5)
    d.loc[restore,'nee']=d.loc[restore,'recovered_nee'];d['dose_stratum']=dose_strata(d.nee)
    er,eq=ratio_engine(d,'eicu','recovered_dose',['all','>1.0']);sensitivity.append(er[~er.domain.isin(['C_bin_fraction'])]);sever.append(severity(d,'eicu','recovered_dose',['all','>1.0']))
    save(er[er.domain.eq('A_distribution')].merge(eq,on=['database','scenario','stratum'],validate='many_to_one'),'phase2_sensitivity_recovered_miss_distribution.csv')
    save(pd.concat(distributions,ignore_index=True),'phase2_A_miss_distribution.csv');save(pd.concat(deep,ignore_index=True),'phase2_A_deep_miss.csv')
    save(pd.concat(B,ignore_index=True),'phase2_B_cuff_perspective.csv');save(pd.concat(C,ignore_index=True),'phase2_C_cross_table.csv');save(pd.concat(sever,ignore_index=True),'phase2_C_arterial_severity.csv')
    ss=pd.concat(sensitivity,ignore_index=True)
    for scenario,z in ss.groupby('scenario'):save(z,f'phase2_sensitivity_{scenario}.csv')
    from standardized import calculate
    calculate();save(pd.DataFrame(manifests),'phase2_input_manifest.csv')
    print('PHASE2_ANALYSES_COMPLETE',flush=True)
if __name__=='__main__':main()
