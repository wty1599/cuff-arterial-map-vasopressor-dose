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

from pathlib import Path
import sys,json
R=Path('@CUFF_WORKSPACE@');A=R/'phase6'
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Arial','font.size':9,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
d=pd.read_csv(A/'phase6_window_rates.csv');strata=['<=0.1','>0.1-0.25','>0.25-0.5','>0.5-1.0','>1.0']
labels=['≤0.1','>0.1–0.25','>0.25–0.5','>0.5–1.0','>1.0']
fig,axes=plt.subplots(3,2,figsize=(7.4,8.5),sharex=True,sharey=True)
metrics=[('missed','Cuff ≥65 given arterial <65 (%)'),('cuff_conditioned','Arterial <65 given cuff ≥65 (%)'),('cuff_65_69','Arterial <65 given cuff 65–69 (%)')]
for j,db in enumerate(['mimic','eicu']):
 for i,(metric,ylab) in enumerate(metrics):
  ax=axes[i,j]
  for k,(w,color,marker,lab) in enumerate([('5min','#0072B2','o','±5 min'),('2min','#D55E00','s','±2 min'),('same_minute','#009E73','^','Same recorded minute')]):
   q=d[(d.database==db)&(d.frame=='common')&(d.metric==metric)].set_index('stratum').loc[strata]
   q=q[q.window==w] if False else d[(d.database==db)&(d.frame=='common')&(d.metric==metric)&(d.window==w)].set_index('stratum').loc[strata]
   y=q.estimate.to_numpy()*100;lo=q.ci_lower.to_numpy()*100;hi=q.ci_upper.to_numpy()*100
   ax.errorbar(np.arange(5)+(k-1)*.15,y,yerr=[y-lo,hi-y],color=color,fmt=marker,ms=4,capsize=2,elinewidth=.8,label=lab,zorder=3)
  ax.set_ylim(0,100);ax.set_yticks([0,25,50,75,100]);ax.grid(axis='y',color='#e1e5e8',lw=.6,zorder=0)
  ax.set_xticks(range(5),labels,rotation=20,ha='right');ax.set_xlim(-.45,4.45)
  ax.text(-.13,1.03,chr(65+i*2+j),transform=ax.transAxes,fontweight='bold',fontsize=12)
  if j==0:ax.set_ylabel(ylab)
  if i==0:ax.set_title('MIMIC-IV' if db=='mimic' else 'eICU',fontweight='bold')
  if i==2:ax.set_xlabel('NE-eq (µg/kg/min)')
handles,labs=axes[0,0].get_legend_handles_labels();fig.legend(handles,labs,loc='upper center',bbox_to_anchor=(.52,.99),ncol=3,frameon=False)
fig.subplots_adjust(top=.91,bottom=.09,left=.105,right=.98,hspace=.22,wspace=.17)
for ext in ['png','svg','pdf']:fig.savefig(A/'figures'/f'phase6_common_windows.{ext}',dpi=600,facecolor='white')
caption='Supplementary Figure 5. Threshold probabilities at common cuff times across pairing windows. Identical cuff readings and dose assignments are used within each database. Points and intervals show pair-weighted proportions and 95% percentile confidence intervals from 2,000 joint patient-cluster bootstrap samples. Horizontal offsets separate windows and do not denote dose differences. Panels A–B show missed hypotension, C–D arterial hypotension given cuff MAP ≥65 mmHg, and E–F arterial hypotension given cuff MAP 65 to <70 mmHg. Only known-dose pairs enter the five dose strata. Common-time eICU values for ±2 minutes and the same recorded minute are identical; its arterial source values remain five-minute medians. Wide highest-dose intervals reflect limited numbers of readings and patients. The figure describes selected dual-monitoring observations, not probabilities validated in cuff-only care.'
(A/'figures/phase6_common_windows_caption.txt').write_text(caption,encoding='utf-8')
print('FIGURE_COMPLETE')
