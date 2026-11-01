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
import sys,csv,json
R=Path('@CUFF_WORKSPACE@');A=R/'phase5'

import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.family':'Arial','font.size':9,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
d=pd.read_csv(A/'phase5_S5_counterfactual.csv');labs=['<=0.1','>0.1-0.25','>0.25-0.5','>0.5-1.0','>1.0'];disp=['≤0.1','>0.1–0.25','>0.25–0.5','>0.5–1.0','>1.0'];trace=[]
fig,axs=plt.subplots(2,2,figsize=(7.4,5.5),sharex=True)
for row,metric in enumerate(['missed_hypotension','art_lt65_given_cuff_ge65']):
 for col,db in enumerate(['mimic','eicu']):
  ax=axs[row,col];s=d[(d.database==db)&(d.analysis=='primary')&(d.metric==metric)].set_index('stratum').loc[labs]
  for j,(key,label,color,marker,off) in enumerate([('observed','Observed','#0072B2','o',-.08),('expected','Expected under a common relation','#D55E00','s',.08)]):
   y=s[key].to_numpy()*100;lo=s[key+'_ci_lower'].to_numpy()*100;hi=s[key+'_ci_upper'].to_numpy()*100
   ax.errorbar(np.arange(5)+off,y,yerr=[y-lo,hi-y],fmt=marker+'-',color=color,markersize=4,linewidth=1,capsize=2,label=label)
   for k,st in enumerate(labs):trace.append(dict(panel=chr(65+row*2+col),database=db,metric=metric,stratum=st,series=key,estimate_percent=y[k],lower_percent=lo[k],upper_percent=hi[k],source='phase5/phase5_S5_counterfactual.csv'))
  ax.set_ylim(0,65 if row==0 else 55);ax.set_xticks(range(5),disp);ax.tick_params(axis='x',labelsize=8);ax.set_xlim(-.4,4.4)
  ax.text(-.15,1.05,chr(65+row*2+col),transform=ax.transAxes,fontweight='bold',fontsize=11)
  if row==0:ax.set_title('MIMIC-IV' if db=='mimic' else 'eICU',pad=10)
  if col==0:ax.set_ylabel('Cuff MAP ≥65 among\narterial MAP <65 (%)' if row==0 else 'Arterial MAP <65 among\ncuff MAP ≥65 (%)')
  if row==1:ax.set_xlabel('NE-eq (µg/kg/min)')
handles,leg=axs[0,0].get_legend_handles_labels();fig.legend(handles,leg,loc='lower center',ncol=2,frameon=False,bbox_to_anchor=(.52,.0))
fig.subplots_adjust(left=.12,right=.98,top=.92,bottom=.16,wspace=.26,hspace=.25)
base=A/'figures/phase5_S5_counterfactual'
for ext in ['png','svg','pdf']:fig.savefig(base.with_suffix('.'+ext),dpi=600)
pd.DataFrame(trace).to_csv(A/'figures/phase5_figure_source_trace.csv',index=False,encoding='utf-8-sig')
caption='Supplementary Figure 4. Observed and expected threshold probabilities under a common measurement relation. Panels A and B show missed hypotension; C and D show arterial hypotension given cuff MAP ≥65 mmHg. Only original known-dose pairs are included. Expected probabilities hold the within-database spline relation between arterial MAP and cuff MAP ≥65 constant across dose strata and retain each stratum’s arterial MAP distribution. Points and bars are estimates and 95% percentile intervals from 2,000 patient-cluster bootstrap samples, refitting the relation in every sample (seed 20261027). Observed intervals here come from the same joint resampling and do not replace intervals from earlier analyses. An expected risk exceeding the observed risk is retained.'
base.with_name(base.name+'_caption.txt').write_text(caption,encoding='utf-8')
print('Figure and 40-point CSV source trace saved.')
