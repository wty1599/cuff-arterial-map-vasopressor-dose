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
import sys, json, hashlib
ROOT=Path('@CUFF_WORKSPACE@')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=ROOT/'phase4/figures';OUT.mkdir(exist_ok=True)
source=ROOT/'phase4/phase4_S4_bias_by_range.csv';d=pd.read_csv(source)
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
 'font.size':7.5,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,
 'legend.frameon':False,'svg.fonttype':'none','pdf.fonttype':42})
fig,axes=plt.subplots(1,2,figsize=(7.2,3.35),sharey=True)
fig.subplots_adjust(left=.105,right=.98,top=.86,bottom=.29,wspace=.20)
groups=[('arterial_MAP','Grouped by arterial MAP','#0072B2','o',-.06),('mean_of_cuff_and_arterial_MAP','Grouped by the mean of both MAPs','#D55E00','s',.06)]
labels=['<55','55–64','65–74','75–84','85–94','≥95'];trace=[]
for i,(ax,db,title) in enumerate(zip(axes,['mimic','eicu'],['MIMIC-IV','eICU'])):
 for key,label,col,mk,offset in groups:
  z=d[d.database.eq(db)&d.axis.eq(key)].copy();assert len(z)==6
  y=z.estimate.to_numpy();lo=z.ci_lower.to_numpy();hi=z.ci_upper.to_numpy()
  ax.errorbar(np.arange(6)+offset,y,yerr=np.vstack([y-lo,hi-y]),fmt=mk+'-',color=col,
              markersize=3.5,linewidth=1.1,capsize=2,elinewidth=.8,label=label)
  for idx,r in z.iterrows():
   trace.append({'panel':title,'source_data_row':int(idx)+1,'database':db,'axis':key,
     'pressure_range':r.pressure_range,'estimate':r.estimate,'ci_lower':r.ci_lower,'ci_upper':r.ci_upper,
     'n_pairs':int(r.n_pairs),'source_csv':'phase4/phase4_S4_bias_by_range.csv'})
 ax.axhline(0,color='#777777',linewidth=.7,linestyle='--',zorder=0)
 ax.set(xticks=np.arange(6),xticklabels=labels,xlim=(-.35,5.35),ylim=(-23,17),xlabel='MAP range (mmHg)')
 ax.set_title(title,fontsize=9,pad=9)
 ax.text(-.04,1.09,'AB'[i],transform=ax.transAxes,fontweight='bold',fontsize=10)
 ax.tick_params(axis='both',length=3,width=.65)
axes[0].set_ylabel('Mean cuff minus arterial MAP (mmHg)')
fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.52,.055),ncol=2,fontsize=7.5)
for ext in ['png','svg','pdf']:fig.savefig(OUT/f'phase4_S4_bias_by_range.{ext}',dpi=600,facecolor='white')
plt.close(fig)
caption=('Supplementary Figure 3. Cuff-minus-arterial MAP differences by pressure range and grouping measurement. '
 'Points show arithmetic pair-weighted mean differences and bars show 95% percentile confidence intervals from 2,000 patient-cluster bootstrap samples (seed 20261027). '
 'Blue circles group pairs by arterial reference MAP; orange squares group the same pairs by the mean of cuff and arterial MAP. '
 'Sample counts are withheld in the repository aggregate release. '
 'Intervals use lower-inclusive, upper-exclusive boundaries, with no rounding before classification. '
 'Bin counts are in Supplementary Table 15. These arithmetic means differ from the mixed-model estimates in Supplementary Table 6. '
 'Grouping by the measurement mean is symmetric and does not establish physiological causation or remove all mathematical coupling.')
(OUT/'phase4_S4_bias_by_range_caption.txt').write_text(caption,encoding='utf-8')
pd.DataFrame(trace).to_csv(OUT/'phase4_figure_source_trace.csv',index=False)
(OUT/'figure_contract.md').write_text('# Figure contract\n\nConclusion: Positive low-pressure differences appear with arterial-reference grouping and disappear with symmetric mean grouping.\n\nTwo-panel quantitative comparison; existing Python/Matplotlib workflow; 7.2 × 3.35 inches; editable SVG text, PDF, 600-dpi PNG; unchanged original pairs. No patient-level points displayed. Source: phase4_S4_bias_by_range.csv.\n',encoding='utf-8')
print('Figure saved; 24 estimates and confidence intervals linked to source rows.')
