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
from scipy.optimize import minimize_scalar

def summarize(df,db,scenario,subgroups=False):
    # Patient clusters are the resampling unit; every row remains pair-weighted.
    pats=pd.Index(sorted(df.patient.unique()));codes=pats.get_indexer(df.patient);n=len(pats)
    masks={'all':np.ones(len(df),bool),'known':df.nee.notna().to_numpy(),'unknown':df.nee.isna().to_numpy()}
    masks.update({x:df.dose_stratum.eq(x).to_numpy() for x in LABELS})
    if subgroups:
        if db in ['mimic','eicu']:
            masks.update({'bmi40_yes':df.bmi40.eq(1).to_numpy(),'bmi40_no':df.bmi40.eq(0).to_numpy()})
        if db in ['mimic','eicu','amsterdam']:
            masks.update({'af_yes':df.af_subgroup.eq(1).to_numpy(),'af_no':df.af_subgroup.eq(0).to_numpy()})
    specs=[];arr=[]
    for endpoint in ['missed','false_positive']:
        eligible=(df.art_map<65).to_numpy() if endpoint=='missed' else (df.art_map>=65).to_numpy()
        event=(df.cuff_map>=65).to_numpy() if endpoint=='missed' else (df.cuff_map<65).to_numpy()
        for label,mask in masks.items():
            denom=eligible&mask;numer=denom&event
            dc=np.bincount(codes,weights=denom,minlength=n);ec=np.bincount(codes,weights=numer,minlength=n)
            specs.append((endpoint,label,int(numer.sum()),int(denom.sum()),int((dc>0).sum()),int(mask.sum()),int(np.unique(codes[mask]).size)))
            arr.extend([ec,dc])
    matrix=np.array(arr).T
    rng=np.random.default_rng(P['seed']);boots=[]
    for b in range(0,P['bootstrap_replicates'],50):
        weights=rng.multinomial(n,np.full(n,1/n),size=min(50,P['bootstrap_replicates']-b))
        boots.append(weights@matrix)
    counts=np.concatenate(boots);ne=counts[:,0::2];de=counts[:,1::2]
    br=np.divide(ne,de,out=np.full_like(ne,np.nan),where=de>0)
    # Continuity correction is restricted to boundary transforms for meta-analysis.
    boundary=(ne==0)|(ne==de)
    meta_p=np.where(boundary,(ne+.5)/(de+1),br)
    bl=logit(meta_p);rows=[]
    for k,(ep,lab,ev,den,npats,allpairs,allpats) in enumerate(specs):
        finite=np.isfinite(br[:,k]);lo,hi=np.quantile(br[finite,k],[.025,.975]) if finite.any() else (np.nan,np.nan)
        ml=bl[finite,k];ms=float(np.std(ml,ddof=1)) if len(ml)>1 else np.nan
        r=ev/den if den else np.nan
        y=float(logit((ev+.5)/(den+1) if ev in [0,den] else r)) if den else np.nan
        degenerate=bool(den and ev in [0,den] and lo==hi)
        rows.append(dict(database=db,scenario=scenario,endpoint=ep,stratum=lab,numerator=ev,denominator=den,denominator_patients=npats,all_pairs=allpairs,all_patients=allpats,estimate=r,ci_lower=lo,ci_upper=hi,bootstrap_valid=int(finite.sum()),precision='boundary_not_informative' if degenerate else 'insufficient' if hi-lo>.30 else 'adequate' if den else 'not_estimable',ci_status='boundary_degenerate' if degenerate else 'ok' if den else 'no_denominator',logit_estimate=y,logit_se=ms,ci_method='patient_cluster_percentile_2000',boundary_events=bool(den and ev in [0,den])))
    comparisons=[]
    for ep in ['missed','false_positive']:
        for name in ['bmi40','af']:
            ids=[i for i,s in enumerate(specs) if s[0]==ep and s[1] in [name+'_yes',name+'_no']]
            if len(ids)!=2: continue
            i=next(i for i in ids if specs[i][1]==name+'_yes');j=next(i for i in ids if specs[i][1]==name+'_no')
            diff=br[:,i]-br[:,j];good=np.isfinite(diff)
            lo,hi=np.quantile(diff[good],[.025,.975]) if good.any() else (np.nan,np.nan)
            comparisons.append(dict(database=db,endpoint=ep,subgroup=name,exposed_rate=rows[i]['estimate'],control_rate=rows[j]['estimate'],risk_difference=rows[i]['estimate']-rows[j]['estimate'],ci_lower=lo,ci_upper=hi,bootstrap_valid=int(good.sum()),method='unadjusted_patient_cluster_difference'))
    return pd.DataFrame(rows),pd.DataFrame(comparisons)

def loa(z):
    v=z.assign(delta=z.cuff_map-z.art_map).groupby('patient').delta.agg(['count','mean','var'])
    n=v['count'].to_numpy(float);m=v['mean'].to_numpy(float);N=n.sum();K=len(n)
    if K<2 or N<=K: return dict(n_pairs=int(N),n_patients=K,status='insufficient_repeated_clusters')
    ss=float(np.sum((n-1)*v['var'].fillna(0).to_numpy()))
    def fit(lam):
        a=n/(1+n*lam);mu=float(np.sum(a*m)/np.sum(a));S=ss+np.sum(a*(m-mu)**2);s2=S/(N-1)
        obj=(N-1)*np.log(max(s2,1e-300))+np.log1p(n*lam).sum()+np.log(a.sum())
        return obj,mu,s2,a.sum()
    opt=minimize_scalar(lambda t:fit(np.exp(t))[0],bounds=(-20,20),method='bounded')
    lam=np.exp(opt.x) if opt.fun<fit(0)[0] else 0
    _,mu,s2,A=fit(lam);tau=s2*lam;sd=np.sqrt(tau+s2);se=np.sqrt(s2/A)
    return dict(n_pairs=int(N),n_patients=K,bias=mu,bias_ci_lower=mu-1.96*se,bias_ci_upper=mu+1.96*se,patient_variance=tau,residual_variance=s2,loa_lower=mu-1.96*sd,loa_upper=mu+1.96*sd,status='ok_boundary_zero' if lam==0 else 'ok',estimation='random_intercept_REML')

def main():
    allmain=[];secondary=[];sens=[];subs=[];comp=[];ba=[];dens=[];model=[]
    for db in DBS:
        d=pd.read_parquet(RES/f'pairs_{db}.parquet')
        main,co=summarize(d,db,'primary_5min',True);allmain.append(main)
        save(main[main.endpoint.eq('missed')&~main.stratum.str.startswith(('bmi40','af'))],f'phase1_primary_{db}.csv')
        secondary.append(main[main.endpoint.eq('false_positive')&~main.stratum.str.startswith(('bmi40','af'))]);subs.append(main[main.stratum.str.startswith(('bmi40','af'))]);comp.append(co)
        for name,file in [('2min','window2'),('same_minute','same_minute')]:
            w=pd.read_parquet(RES/f'{file}_{db}.parquet')
            sd=d.drop(columns=['art_map','art_values','nearest_diff']).merge(w,on='cid',validate='one_to_one')
            a,_=summarize(sd,db,name);sens.append(a)
        first=d[d.art_map<65].sort_values(['patient','t','cid']).drop_duplicates('patient')
        a,_=summarize(first,db,'first_low_pair');sens.append(a[a.endpoint.eq('missed')])
        sd=d[~d.quality_any].copy();a,_=summarize(sd,db,'exclude_flagged');sens.append(a)
        if db=='eicu':
            a,_=summarize(d[d.nee.notna()].copy(),db,'eicu_known_dose');sens.append(a)
        d['map_bin']=pd.cut(d.art_map,[-np.inf,55,65,75,85,np.inf],right=False,labels=['<55','55-64','65-74','75-84','>=85']).astype(str)
        groups=[('overall','all',d)]
        groups += [('arterial_map',str(k),v) for k,v in d.groupby('map_bin',observed=True)]
        groups += [('dose',str(k),d[d.dose_stratum.eq(k)]) for k in [*LABELS,'unknown']]
        groups += [('arterial_map_x_dose',str(k[0])+' | '+str(k[1]),v) for k,v in d.groupby(['map_bin','dose_stratum'],observed=True)]
        for kind,label,z in groups: ba.append(dict(database=db,stratification=kind,stratum=label,**loa(z)))
        # All displayed BA points represent aggregated cells, never individual observations.
        mean=(d.cuff_map+d.art_map)/2;delta=d.cuff_map-d.art_map
        hist,xe,ye=np.histogram2d(mean,delta,bins=[np.arange(20,202,2),np.arange(-180,182,2)])
        for i,j in zip(*np.where(hist>0)):
            dens.append(dict(database=db,mean_map=(xe[i]+xe[i+1])/2,difference=(ye[j]+ye[j+1])/2,n_pairs=int(hist[i,j])))
        if db in P['model_databases']:
            m=d[d.art_map<65][['patient','cid','database','nee','art_map','age','sex','bmi40','af_model','cuff_map']].copy();m['missed']=(m.cuff_map>=65).astype(int);model.append(m)
        print(json.dumps({'database':db,'rates_and_bland_altman':'complete'}),flush=True)
    main=pd.concat(allmain,ignore_index=True)
    save(main[~main.stratum.str.startswith(('bmi40','af'))],'phase1_all_rates.csv')
    save(pd.concat(secondary,ignore_index=True),'phase1_secondary_false_positive.csv')
    save(pd.concat(subs,ignore_index=True),'phase1_secondary_subgroups.csv')
    save(pd.concat(comp,ignore_index=True),'phase1_secondary_subgroup_comparisons.csv')
    ss=pd.concat(sens,ignore_index=True)
    for name,z in ss.groupby('scenario'):save(z,f'phase1_sensitivity_{name}.csv')
    save(pd.DataFrame(ba),'phase1_secondary_bland_altman.csv');save(pd.DataFrame(dens),'phase1_bland_altman_density.csv')
    pd.concat(model,ignore_index=True).to_csv(RES/'model_input.csv',index=False)
    print('RATES_COMPLETE',flush=True)

if __name__=="__main__": main()
