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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import PercentFormatter
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'legend.frameon':False,'svg.fonttype':'none','pdf.fonttype':42})
COL={'mimic':'#285C7A','eicu':'#B45A3C','pooled':'#222222'}
names={'mimic':'MIMIC-IV','eicu':'eICU','amsterdam':'AmsterdamUMCdb','sicdb':'SICdb'}
export=[]
def finish(fig,name,caption,source):
    for ext in ['png','svg','pdf']:
        path=OUT/'figures'/f'{name}.{ext}';fig.savefig(path,dpi=300,facecolor='white',bbox_inches='tight')
    (OUT/'figures'/f'{name}_caption.txt').write_text(caption,encoding='utf-8')
    export.append(dict(figure=name,source_csv=source,formats='png;svg;pdf',backend='matplotlib'))
    plt.close(fig)

rates=pd.read_csv(OUT/'phase1_all_rates.csv');pool=pd.read_csv(OUT/'phase1_pooled.csv')
fig=plt.figure(figsize=(9,7.3));gs=fig.add_gridspec(3,1,height_ratios=[5.5,.85,1.05],hspace=.6)
ax=fig.add_subplot(gs[0]);x=np.arange(5)
for db,dx,mk in [('mimic',-.18,'o'),('eicu',0,'s'),('pooled',.18,'D')]:
    z=pool[pool.endpoint.eq('missed')] if db=='pooled' else rates[rates.database.eq(db)&rates.endpoint.eq('missed')]
    z=z.set_index('stratum').loc[LABELS];y=z.estimate.to_numpy()*100
    ax.errorbar(x+dx,y,yerr=np.array([z.estimate-z.ci_lower,z.ci_upper-z.estimate])*100,fmt=mk,color=COL[db],label='Random-effects pooled' if db=='pooled' else names[db],capsize=3,markersize=5,linewidth=1.2)
ax.set(xticks=x,xticklabels=['≤0.1','>0.1–0.25','>0.25–0.5','>0.5–1.0','>1.0'],ylim=(0,70),xlim=(-.55,4.55),ylabel='Cuff missed hypotension (%)',xlabel='Norepinephrine-equivalent dose (µg/kg/min)')
ax.grid(axis='y',alpha=.18);ax.set_title('Cuff MAP ≥65 among pairs with arterial MAP <65',loc='left',pad=14,fontweight='bold');ax.legend(loc='lower left',ncol=3,fontsize=8)
ax.text(0,1.02,'Mapping review pending · estimates and 95% confidence intervals',transform=ax.transAxes,fontsize=8,color='#555555')
a=fig.add_subplot(gs[1]);a.set(xlim=(0,3),ylim=(0,1),yticks=[],xticks=[0,.5,1,2,3]);a.axvspan(.5,1,color='#A7C4B5',alpha=.7);a.spines['left'].set_visible(False);a.text(.75,.58,'Common high-dose range',ha='center',fontsize=8);a.set_xlabel('NE-eq reference scale (µg/kg/min)',fontsize=8)
b=fig.add_subplot(gs[2]);b.set(xlim=(0,3),ylim=(0,1),yticks=[],xticks=[0,1,2,2.5,3]);b.spines['left'].set_visible(False);b.axvline(2.5,color='#8E5572',ls='--',lw=1.2);b.text(2.43,.65,'EVERDAC exclusion threshold: >2.5',ha='right',fontsize=8)
b.set_xlabel('EVERDAC scale: norepinephrine tartrate + epinephrine (µg/kg/min)',fontsize=8)
fig.text(.125,.005,'The two reference scales are not interchangeable. The >1.0 NE-eq stratum is open ended.\nPooled estimates include supplementary databases; SICdb has no observations in the >1.0 stratum.',fontsize=8,color='#444444')
finish(fig,'phase1_primary_dose','Main databases and random-effects pooling of all eligible databases. Cuff missed hypotension is cuff MAP >=65 among arterial MAP <65. Database CIs use 2,000 patient-cluster bootstrap replicates. Pooled CIs use logit REML with modified Hartung-Knapp adjustment. The protocol dose ruler is distinct from NE-eq: EVERDAC used norepinephrine tartrate plus epinephrine; 2 micrograms norepinephrine tartrate equals 1 microgram base. Mapping review pending.','phase1_primary_mimic.csv;phase1_primary_eicu.csv;phase1_pooled.csv')

cur=pd.read_csv(OUT/'phase1_secondary_model_curve.csv');pts=pd.read_csv(OUT/'phase1_secondary_model_predictions.csv')
fig,ax=plt.subplots(figsize=(8.5,5.1),layout='constrained');z=cur[cur.nee<=2.5]
ax.axvspan(.5,1,color='#A7C4B5',alpha=.25,label='Common high-dose range')
ax.fill_between(z.nee,z.ci_lower*100,z.ci_upper*100,color='#285C7A',alpha=.18,label='95% confidence interval')
ax.plot(z.nee,z.estimate*100,color='#285C7A',lw=2,label='Adjusted missed hypotension')
ax.scatter(pts.nee,pts.estimate*100,color='#285C7A',s=22,zorder=3)
ax.set(xlim=(0,2.5),ylim=(0,65),xlabel='Norepinephrine-equivalent dose (µg/kg/min)',ylabel='Adjusted cuff missed hypotension (%)')
ax.set_title('Dose association after covariate adjustment',loc='left',fontweight='bold');ax.grid(axis='y',alpha=.16);ax.legend(loc='lower left',fontsize=8)
fig.supxlabel('MIMIC-IV + eICU · patient random intercept · marginal standardization\n95% delta-method CI; random-effect variance treated as plug-in. Mapping review pending.',fontsize=8)
finish(fig,'phase1_dose_curve','Association model in complete-case low-arterial-MAP pairs from MIMIC-IV and eICU. Three-knot restricted cubic spline, adjusted for arterial MAP, age, sex, BMI >40, AF ascertainment, and database; patient random intercept. Marginal standardization integrates the fitted normal random intercept. Delta-method CIs use fixed-effect covariance with random-effect variance treated as plug-in. This figure displays doses 0-2.5; the full evaluated curve to the observed maximum is in the source CSV. It does not estimate the effect or safety of monitoring strategies.','phase1_secondary_model_curve.csv;phase1_secondary_model_predictions.csv')

ba=pd.read_csv(OUT/'phase1_secondary_bland_altman.csv');den=pd.read_csv(OUT/'phase1_bland_altman_density.csv')
if 'n_patients' not in den.columns:
    raise RuntimeError('Regenerate S6 with patient-aware suppression before plotting')
den = den.loc[(den.n_pairs >= 10) & (den.n_patients >= 10)].copy()
fig,axs=plt.subplots(2,2,figsize=(9,7.2),sharex=True,sharey=True,layout='constrained')
for ax,db in zip(axs.flat,DBS):
    z=den[den.database.eq(db)&den.n_pairs.ge(5)];xb=np.arange(20,202,2);yb=np.arange(-180,182,2);h=np.full((90,180),np.nan)
    for row in z.itertuples():
        i=int(round((row.mean_map-21)/2));j=int(round((row.difference+179)/2))
        if 0<=i<90 and 0<=j<180:h[i,j]=row.n_pairs
    img=ax.pcolormesh(xb,yb,h.T,cmap='cividis',norm=LogNorm(vmin=5,vmax=max(6,z.n_pairs.max() if len(z) else 6)),rasterized=True)
    r=ba[ba.database.eq(db)&ba.stratification.eq('overall')].iloc[0]
    ax.axhline(r.bias,color='#B45A3C',lw=1.3);ax.axhline(r.loa_lower,color='#B45A3C',ls='--',lw=.9);ax.axhline(r.loa_upper,color='#B45A3C',ls='--',lw=.9)
    ax.set_title(names[db],loc='left',fontsize=10);ax.set(xlim=(30,140),ylim=(-70,70))
    ax.text(.03,.96,f'Bias {r.bias:.1f}; LoA {r.loa_lower:.1f} to {r.loa_upper:.1f} mmHg',transform=ax.transAxes,va='top',fontsize=7.5,bbox={'facecolor':'white','edgecolor':'none','alpha':.85})
    fig.colorbar(img,ax=ax,shrink=.75,pad=.02,label='Pairs per bin')
for ax in axs[1,:]:ax.set_xlabel('Mean of cuff and arterial MAP (mmHg)')
for ax in axs[:,0]:ax.set_ylabel('Cuff − arterial MAP (mmHg)')
fig.suptitle('Repeated-measures Bland–Altman analysis',x=.02,ha='left',fontweight='bold',fontsize=12)
fig.supxlabel('REML patient random-intercept LoA; density bins contain ≥5 pairs.\nDisplay limits crop the density only; all eligible pairs contribute to estimation. Mapping review pending.',fontsize=8)
finish(fig,'phase1_bland_altman','Supplementary repeated-measures Bland-Altman display. Bias and 95% population limits of agreement derive from a REML patient random-intercept model, bias +/-1.96*sqrt(patient variance + residual variance). The density display uses 2-mmHg bins and shows only bins containing at least 10 pairs from at least 10 distinct patients. Axis cropping does not exclude pairs from estimation. Limits of agreement are not confidence intervals. MAP- and dose-stratified results appear in the accompanying CSV.','phase1_secondary_bland_altman.csv;phase1_bland_altman_density.csv')
save(pd.DataFrame(export),'phase1_figure_sources.csv')
print('FIGURES_COMPLETE',flush=True)
