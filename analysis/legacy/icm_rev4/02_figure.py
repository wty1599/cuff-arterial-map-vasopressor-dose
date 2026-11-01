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
import sys,csv,json,hashlib
R=Path('@CUFF_WORKSPACE@');O=R/'phase3/icm_revision4';F=O/'figures'
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
contract=dict(conclusion='Among selected paired observations, a cuff reading above target carries greater residual arterial-hypotension risk at high doses; a common measurement relation reproduces much of this pattern.',panels={'A':'MIMIC observed and expected conditional risks plus background hypotension','B':'eICU observed and expected risks plus background hypotension','C':'MIMIC cuff-range risks in three dose groups','D':'eICU cuff-range risks in three dose groups'},archetype='quantitative grid',backend='Python Matplotlib, continuing the existing phase5 figure workflow',size_inches=[7.4,6.5],exports=['editable-text SVG','600 dpi PNG','vector PDF'],statistics='Patient-cluster 95% intervals from existing validated CSVs; no model refit; upper observed and expected intervals share phase5 bootstrap; prevalence points are descriptive.')
(F/'figure_contract.json').write_text(json.dumps(contract,indent=2),encoding='utf-8')
matplotlib.rcParams.update({'font.family':'Arial','font.size':9,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.8,'legend.frameon':False})
s5=pd.read_csv(R/'phase5/phase5_S5_counterfactual.csv');s3path=R/'phase3/revision4/data/S3_cuff_conditioned_probabilities.csv';s3=pd.read_csv(s3path)
strata=['<=0.1','>0.1-0.25','>0.25-0.5','>0.5-1.0','>1.0'];labels=['≤0.1','>0.1–0.25','>0.25–0.5','>0.5–1.0','>1.0'];bins=['65-69','70-74','75-79','>=80'];trace=[]
fig,axs=plt.subplots(2,2,figsize=(7.4,6.5))
for col,(db,dbname) in enumerate([('mimic','MIMIC-IV'),('eicu','eICU')]):
 ax=axs[0,col];q=s5[(s5.database==db)&(s5.analysis=='primary')&(s5.metric=='art_lt65_given_cuff_ge65')].set_index('stratum').loc[strata]
 for key,label,color,marker,offset,face in [('observed','Observed','#0072B2','o',-.08,'#0072B2'),('expected','Common-relation expected','#D55E00','s',.08,'white')]:
  y=q[key].to_numpy()*100;lo=q[key+'_ci_lower'].to_numpy()*100;hi=q[key+'_ci_upper'].to_numpy()*100
  ax.errorbar(np.arange(5)+offset,y,yerr=[y-lo,hi-y],fmt=marker+'-',color=color,mfc=face,mec=color,markersize=4,linewidth=1,capsize=2,label=label)
  for st in strata:
   r=q.loc[st];idx=int(s5[(s5.database==db)&(s5.analysis=='primary')&(s5.metric=='art_lt65_given_cuff_ge65')&(s5.stratum==st)].index[0])+1
   trace.append(dict(panel=chr(65+col),database=db,series=key,stratum=st,category='',estimate=float(r[key]),ci_lower=float(r[key+'_ci_lower']),ci_upper=float(r[key+'_ci_upper']),source_csv='phase5/phase5_S5_counterfactual.csv',source_data_row=idx,source_columns=f'{key};{key}_ci_lower;{key}_ci_upper'))
 prev=q.arterial_hypotension_prevalence.to_numpy()*100
 ax.plot(np.arange(5),prev,'D:',color='#8C8C8C',markersize=3.5,linewidth=.9,label='Arterial hypotension prevalence')
 for st,p in zip(strata,prev):
  idx=int(s5[(s5.database==db)&(s5.analysis=='primary')&(s5.metric=='art_lt65_given_cuff_ge65')&(s5.stratum==st)].index[0])+1
  trace.append(dict(panel=chr(65+col),database=db,series='prevalence',stratum=st,category='',estimate=p/100,ci_lower='',ci_upper='',source_csv='phase5/phase5_S5_counterfactual.csv',source_data_row=idx,source_columns='arterial_hypotension_prevalence'))
 ax.set_xticks(range(5),labels,fontsize=8);ax.set_xlabel('NE-eq (µg/kg/min)');ax.set_title(dbname,pad=12);ax.set_ylim(0,65);ax.set_xlim(-.4,4.4)
 ax=axs[1,col]
 for st,lab,color,m in [('<=0.25','≤0.25','#0072B2','o'),('>0.25-1.0','>0.25–1.0','#D55E00','s'),('>1.0','>1.0','#009E73','^')]:
  d=s3[(s3.database==dbname)&(s3.scenario=='Primary five-minute pairing')&(s3.stratum=='dosegroup:'+st)&(s3.category.isin(bins))].set_index('category').loc[bins]
  y=d.estimate.to_numpy()*100;lo=d.ci_lower.to_numpy()*100;hi=d.ci_upper.to_numpy()*100
  ax.errorbar(np.arange(4),y,yerr=[y-lo,hi-y],fmt=m+'-',color=color,markersize=4,linewidth=1,capsize=2,label=lab)
  for cat in bins:
   r=d.loc[cat];idx=int(s3[s3.row_id==r.row_id].index[0])+1
   trace.append(dict(panel=chr(67+col),database=db,series=st,stratum=st,category=cat,estimate=float(r.estimate),ci_lower=float(r.ci_lower),ci_upper=float(r.ci_upper),source_csv=s3path.relative_to(R).as_posix(),source_data_row=idx,source_columns='estimate;ci_lower;ci_upper'))
 ax.set_xticks(range(4),['65–69','70–74','75–79','≥80']);ax.set_xlabel('Cuff MAP (mmHg)');ax.set_ylim(0,75);ax.set_xlim(-.2,3.2)
for row in range(2):
 for col in range(2):
  ax=axs[row,col];ax.text(-.13,1.04,chr(65+row*2+col),transform=ax.transAxes,fontweight='bold',fontsize=12)
  ax.yaxis.grid(True,color='#E8E8E8',linewidth=.5);ax.set_axisbelow(True)
  if col==0:ax.set_ylabel('Arterial MAP <65 mmHg (%)')
handles,labels0=axs[0,0].get_legend_handles_labels();bylabel=dict(zip(labels0,handles));order=['Observed','Common-relation expected','Arterial hypotension prevalence'];fig.legend([bylabel[x] for x in order],order,loc='upper center',bbox_to_anchor=(.54,.53),ncol=3,fontsize=7.5,columnspacing=1.3)
handles,labels0=axs[1,0].get_legend_handles_labels();fig.legend(handles,labels0,loc='lower center',bbox_to_anchor=(.54,.015),ncol=3,title='NE-eq (µg/kg/min)',title_fontsize=8,fontsize=8)
fig.subplots_adjust(left=.105,right=.98,top=.925,bottom=.155,wspace=.23,hspace=.66)
for ext in ['png','svg','pdf']:fig.savefig(F/f'ICM_rev4_Figure1.{ext}',dpi=600)
pd.DataFrame(trace).to_csv(F/'figure_number_trace.csv',index=False,encoding='utf-8-sig')
caption='Figure 1. Residual arterial hypotension across dose and cuff readings. Upper panels: observed (filled) and common-relation expected (open) probabilities given cuff MAP ≥65 mmHg; grey diamonds show arterial hypotension prevalence. Lower panels: probabilities by cuff interval and dose group. Bars indicate 95% patient-cluster bootstrap intervals. Sample counts are withheld in the repository aggregate release. NE-eq, norepinephrine-equivalent dose; MAP, mean arterial pressure.'
(F/'Figure1_caption.txt').write_text(caption,encoding='utf-8')
print(json.dumps(dict(points=len(trace),caption_words=len(caption.split()),figure=str(F/'ICM_rev4_Figure1.png'))))
