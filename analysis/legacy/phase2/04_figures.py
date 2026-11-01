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
from matplotlib.lines import Line2D
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'legend.frameon':False,'svg.fonttype':'none','pdf.fonttype':42})
COL={'mimic':'#285C7A','eicu':'#B45A3C','pooled':'#222222'}
exports=[]
def finish(fig,name,caption,source):
 for ext in ['png','svg','pdf']:fig.savefig(OUT/'figures'/f'{name}.{ext}',dpi=300,bbox_inches='tight',facecolor='white')
 (OUT/'figures'/f'{name}_caption.txt').write_text(caption,encoding='utf-8');plt.close(fig)
 exports.append(dict(figure=name,source_csv=source,backend='Python/Matplotlib',formats='png;svg;pdf',mapping_status='pending'))

deep=pd.read_csv(OUT/'phase2_A_deep_miss.csv');pool=pd.read_csv(OUT/'phase2_A_deep_pooled.csv')
fig=plt.figure(figsize=(11,6.8));gs=fig.add_gridspec(3,3,height_ratios=[5.0,.65,1.0],hspace=.8,wspace=.14)
tick=['≤0.1','>0.1–\n0.25','>0.25–\n0.5','>0.5–\n1.0','>1.0']
for j,db in enumerate(['mimic','eicu','pooled']):
 ax=fig.add_subplot(gs[0,j]);base=pd.read_csv(P1/'phase1_pooled.csv') if db=='pooled' else pd.read_csv(P1/f'phase1_primary_{db}.csv')
 base=base[base.endpoint.eq('missed')] if db=='pooled' else base
 dd=pool[pool.metric.eq('deep_miss_art_lt55')] if db=='pooled' else deep[deep.database.eq(db)&deep.metric.eq('deep_miss_art_lt55')]
 for q,shift,style,marker,lab in [(base,-.055,'-','o','Arterial MAP <65'),(dd,.055,'--','^','Arterial MAP <55')]:
  z=q.set_index('stratum').loc[LABELS];y=z.estimate.to_numpy()*100
  ax.errorbar(np.arange(5)+shift,y,yerr=np.array([z.estimate-z.ci_lower,z.ci_upper-z.estimate])*100,color=COL[db],ls=style,marker=marker,markersize=4.5,linewidth=1.2,capsize=2,label=lab)
 ax.set(xticks=np.arange(5),xticklabels=tick,ylim=(0,70),xlim=(-.4,4.4));ax.grid(axis='y',alpha=.18)
 ax.set_title('Random-effects pooled' if db=='pooled' else NAMES[db],color=COL[db],fontweight='bold',fontsize=10)
 if j==0:ax.set_ylabel('Cuff MAP ≥65 among low arterial MAP pairs (%)')
 else:ax.tick_params(labelleft=False)
 ax.set_xlabel('NE-eq (µg/kg/min)',fontsize=8)
fig.legend(handles=[Line2D([0],[0],color='#333333',ls='-',marker='o',label='Missed hypotension: arterial MAP <65'),Line2D([0],[0],color='#333333',ls='--',marker='^',label='Deep missed hypotension: arterial MAP <55')],loc='upper center',bbox_to_anchor=(.5,.97),ncol=2,fontsize=8)
a=fig.add_subplot(gs[1,:]);a.set(xlim=(0,3),ylim=(0,1),yticks=[],xticks=[0,.5,1,2,3]);a.spines['left'].set_visible(False);a.axvspan(.5,1,color='#A7C4B5',alpha=.6);a.text(.75,.5,'Common high-dose range',ha='center',va='center',fontsize=8);a.set_xlabel('NE-eq reference scale (µg/kg/min)',fontsize=8)
b=fig.add_subplot(gs[2,:]);b.set(xlim=(0,3),ylim=(0,1),yticks=[],xticks=[0,1,2,2.5,3]);b.spines['left'].set_visible(False);b.axvline(2.5,color='#8E5572',ls='--',lw=1.2);b.text(2.43,.56,'EVERDAC reference: 2.5 tartrate ≈1.25 base',ha='right',fontsize=8);b.set_xlabel('Pure-norepinephrine conversion reference: tartrate (µg/kg/min)',fontsize=8)
fig.text(.125,.005,'Mapping review pending · 95% patient-cluster bootstrap CIs; logit REML / modified Hartung–Knapp pooling.\nEVERDAC summed norepinephrine tartrate + epinephrine; the base conversion is approximate for norepinephrine alone, not a generic NE-eq cutoff.',fontsize=7.5,color='#444444')
base_den=[]
for db in MAIN:
 q=deep[deep.database.eq(db)&deep.metric.eq('deep_miss_art_lt55')&deep.stratum.isin(LABELS)]
 base_den.append(NAMES[db]+': '+', '.join(f'{r.stratum} n={r.denominator}' for r in q.itertuples()))
finish(fig,'phase2_primary_dose_updated','Cuff MAP >=65 conditional on arterial MAP <65 (solid) or <55 (dashed). Each curve has a different arterial-MAP denominator. The solid estimates are unchanged phase1 results. Database CIs: 2,000 patient-cluster percentile bootstrap replicates, seed20261027. Pooling uses known-dose pairs, logit REML and modified Hartung-Knapp, with eligible supplementary databases. Deep-denominator counts: '+'; '.join(base_den)+'. EVERDAC used norepinephrine tartrate plus epinephrine >2.5; approximately1.25 base is a pure-norepinephrine conversion, not a multi-drug NE-eq cutoff. Mapping review pending.','phase1_primary_mimic.csv;phase1_primary_eicu.csv;phase1_pooled.csv;phase2_A_deep_miss.csv;phase2_A_deep_pooled.csv')

B=pd.read_csv(OUT/'phase2_B_cuff_perspective.csv');fig,axes=plt.subplots(1,2,figsize=(10,5.2),sharey=True,layout='constrained')
groups=['dosegroup:<=0.25','dosegroup:>0.25-1.0','dosegroup:>1.0'];labels=['NE-eq ≤0.25','NE-eq >0.25–1.0','NE-eq >1.0'];colors=['#285C7A','#B45A3C','#222222'];markers=['o','s','^'];bins=['65-69','70-74','75-79','>=80'];sourcecounts=[]
for ax,db in zip(axes,MAIN):
 for g,label,color,marker in zip(groups,labels,colors,markers):
  q=B[B.database.eq(db)&B.stratum.eq(g)&B.metric.eq('art_lt65_given_cuff_bin')].set_index('category').loc[bins]
  ax.errorbar(np.arange(4),q.estimate*100,yerr=np.array([q.estimate-q.ci_lower,q.ci_upper-q.estimate])*100,color=color,marker=marker,linewidth=1.4,capsize=2,markersize=4.5,label=label)
  pr=B[B.database.eq(db)&B.stratum.eq(g)&B.metric.eq('arterial_lt65_prevalence')].iloc[0]
  sourcecounts.append(NAMES[db]+' '+label+': low arterial MAP prevalence '+f'{pr.estimate*100:.1f}%; cuff-bin denominators '+','.join(str(x) for x in q.denominator))
 ax.set(xticks=np.arange(4),xticklabels=['65–69','70–74','75–79','≥80'],xlabel='Observed cuff MAP (mmHg)',ylim=(0,70));ax.set_title(NAMES[db],fontweight='bold',fontsize=10);ax.grid(axis='y',alpha=.18);ax.legend(loc='upper right',fontsize=8)
 if db=='mimic':ax.set_ylabel('Paired arterial MAP <65 (%)')
 # Dose-group prevalence is included directly with each panel for interpretability.
 foot=[]
 for g,label in zip(groups,labels):
  r=B[B.database.eq(db)&B.stratum.eq(g)&B.metric.eq('arterial_lt65_prevalence')].iloc[0];foot.append(label+f': {r.estimate*100:.1f}%')
 ax.text(.03,.03,'Arterial MAP <65 in each dose group:\n'+'\n'.join(foot),transform=ax.transAxes,fontsize=7.5,va='bottom',bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
fig.suptitle('Probability of low arterial MAP when cuff MAP appears normal',fontsize=12,fontweight='bold')
fig.supxlabel('Mapping review pending · original known-dose pairs · patient-cluster 95% CIs\nProbabilities depend on the dose-group prevalence of low arterial MAP; no cross-database pooling.',fontsize=8)
finish(fig,'phase2_cuff_perspective','Probability that paired arterial MAP <65 within each normal cuff-MAP bin, by three original known-dose groups. These are observation-level conditional probabilities, not probabilities of an unobserved clinical outcome. Database panels are separate. 95% CIs use 2,000 patient-cluster bootstrap replicates, seed20261027. Dose-group arterial<65 prevalence and bin denominators: '+'; '.join(sourcecounts)+'. Mapping review pending.','phase2_B_cuff_perspective.csv')
save(pd.DataFrame(exports),'phase2_figure_sources.csv');print('PHASE2_FIGURES_COMPLETE')
