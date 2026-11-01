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

"""Direct standardization with independent database bootstrap child streams."""
from common import *

def calculate():
    bins=['60-64','55-59','<55'];children=np.random.SeedSequence(P['seed']).spawn(2)
    datasets={};point={};boot={};population={};popboot={}
    for pos,db in enumerate(MAIN):
        d=read_pairs(db);datasets[db]=d
        pats=pd.Index(sorted(d.patient.unique()));code=pats.get_indexer(d.patient);N=len(pats)
        art=d.art_map.to_numpy();normal=d.cuff_map.ge(65).to_numpy()
        bm=[(art>=60)&(art<65),(art>=55)&(art<60),art<55]
        vectors=[np.bincount(code,weights=m,minlength=N) for m in bm]
        mappings={}
        scenarios=[('primary_5min',d,LABELS)]
        if db=='eicu':
            r=pd.read_csv(INP/'eicu_recovered_pairs.csv',usecols=['cid','recovered_nee']);rec=d.merge(r,on='cid',validate='one_to_one')
            use=rec.nee_raw.isna()&rec.recovered_nee.notna()&~rec.dose_weight_flag&rec.recovered_nee.between(0,5)
            rec.loc[use,'nee']=rec.loc[use,'recovered_nee'];rec['dose_stratum']=dose_strata(rec.nee)
            assert rec.cid.equals(d.cid)
            scenarios.append(('recovered_dose',rec,['all','>1.0']))
        for scenario,z,levels in scenarios:
            for lab in levels:
                mask=np.ones(len(z),bool) if lab=='all' else z.dose_stratum.eq(lab).to_numpy()
                for cat,cell in zip(bins,bm):
                    den=mask&cell;num=den&normal
                    mappings[(scenario,lab,cat)]=(len(vectors),len(vectors)+1)
                    vectors.extend([np.bincount(code,weights=num,minlength=N),np.bincount(code,weights=den,minlength=N)])
        matrix=np.array(vectors).T;rng=np.random.default_rng(children[pos]);rep=[]
        for b in range(0,2000,50):
            w=rng.multinomial(N,np.full(N,1/N),size=min(50,2000-b));rep.append(w@matrix)
        rep=np.concatenate(rep);pt=matrix.sum(axis=0)
        population[db]=pt[:3];popboot[db]=rep[:,:3]
        point[db]={key:pt[list(ids)] for key,ids in mappings.items()}
        boot[db]={key:rep[:,list(ids)] for key,ids in mappings.items()}
    counts=sum(population.values());weight=counts/counts.sum();pop=sum(popboot.values());bw=pop/pop.sum(axis=1)[:,None]
    save(pd.DataFrame(dict(arterial_bin=bins,n_pairs=counts.astype(int),weight=weight,standard_population='main_databases_all_low_arterial_pairs_including_unknown_dose')),'phase2_C_standard_weights.csv')
    sev=pd.read_csv(OUT/'phase2_C_arterial_severity.csv');rows=[];interpret=[]
    for db in MAIN:
        scenarios=[('primary_5min',LABELS)]
        if db=='eicu':scenarios.append(('recovered_dose',['all','>1.0']))
        for scenario,levels in scenarios:
            estimates=[];raw=[]
            for lab in levels:
                pp=np.array([point[db][scenario,lab,cat] for cat in bins]);num=pp[:,0];den=pp[:,1]
                rate=np.column_stack([np.divide(boot[db][scenario,lab,cat][:,0],boot[db][scenario,lab,cat][:,1],out=np.full(2000,np.nan),where=boot[db][scenario,lab,cat][:,1]>0) for cat in bins])
                bb=np.sum(bw*rate,axis=1);good=np.isfinite(bb);lower,upper=np.quantile(bb[good],[.025,.975]) if good.any() else (np.nan,np.nan)
                estimate=float(np.sum(weight*num/den)) if np.all(den>0) else np.nan;crude=num.sum()/den.sum() if den.sum() else np.nan
                n_pat=sev[sev.database.eq(db)&sev.scenario.eq(scenario)&sev.stratum.eq(lab)].iloc[0].patients
                rows.append(dict(database=db,scenario=scenario,stratum=lab,estimate=estimate,ci_lower=lower,ci_upper=upper,raw_missed_rate=crude,low_arterial_pairs=int(den.sum()),denominator_patients=int(n_pat),minimum_cell_pairs=int(den.min()),small_cell=bool((den<30).any()),wide_ci=bool(upper-lower>.30),bootstrap_valid=int(good.sum()),precision='insufficient' if (den<30).any() or upper-lower>.30 else 'adequate' if np.isfinite(estimate) else 'not_estimable',method='direct_standardization_independent_database_patient_cluster_percentile_2000'))
                estimates.append(estimate);raw.append(crude)
            if scenario=='primary_5min':
                span=max(estimates)-min(estimates);rawspan=max(raw)-min(raw);rule=span<=.05 and span<rawspan
                interpret.append(dict(database=db,raw_span=rawspan,standardized_span=span,within_5_percentage_points=bool(span<=.05),shrunk=bool(span<rawspan),rule_satisfied=bool(rule),interpretation='最高层的较低漏判率主要反映该层低血压更深' if rule else '仅报告标化后差距，不作进一步解释'))
    save(pd.DataFrame(rows),'phase2_C_standardized.csv');save(pd.DataFrame(interpret),'phase2_C_interpretation.csv')
    save(pd.DataFrame([dict(database=db,root_seed=P['seed'],child_spawn_key=str(children[i].spawn_key),replicates=2000,independent_database_child_stream=True) for i,db in enumerate(MAIN)]),'phase2_C_bootstrap_streams.csv')
    np.savez_compressed(RES/'C_bootstrap_standard_weights.npz',weights=bw)
    print('C_INDEPENDENT_BOOTSTRAP_COMPLETE',flush=True)

if __name__=='__main__':calculate()
