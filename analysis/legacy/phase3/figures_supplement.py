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
ROOT=Path('@CUFF_WORKSPACE@')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
OUT=ROOT/'phase3/supplement/figures';OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
 'font.size':7.5,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,
 'legend.frameon':False,'svg.fonttype':'none','pdf.fonttype':42})
NAMES={'mimic':'MIMIC-IV','eicu':'eICU','amsterdam':'AmsterdamUMCdb','sicdb':'SICdb'}
sources=[];checks=[]
def ck(n,ok,d=''):
 checks.append(dict(check=n,passed=bool(ok),detail=str(d)))
 if not ok:raise AssertionError(n+': '+str(d))
def load(name):
 d=pd.read_csv(ROOT/'phase1'/name);d['_row']=np.arange(len(d))+2;d['_source']='phase1/'+name;return d
def tr(name,location,r,keys,fields):
 def cv(v):return v.item() if isinstance(v,(np.integer,np.floating)) else v
 sources.append(dict(artifact=name,location=location,source_csv=r['_source'],source_csv_row=int(r['_row']),
   source_key=json.dumps({k:cv(r[k]) for k in keys},sort_keys=True),
   source_fields=';'.join(fields),source_values=json.dumps({k:cv(r[k]) for k in fields},sort_keys=True)))
def save(fig,name,caption):
 for ext in ['png','svg','pdf']:fig.savefig(OUT/f'{name}.{ext}',dpi=400,facecolor='white')
 plt.close(fig);(OUT/f'{name}_caption.txt').write_text(caption+'\n',encoding='utf-8')

ba=load('phase1_secondary_bland_altman.csv');density=load('phase1_bland_altman_density.csv')
if 'n_patients' not in density.columns:
    raise RuntimeError('Regenerate S6 with patient-aware suppression before plotting')
density = density.loc[(density.n_pairs >= 10) & (density.n_patients >= 10)].copy()
fig,axes=plt.subplots(2,2,figsize=(7.2,6.0),sharex=True,sharey=True)
fig.subplots_adjust(left=.09,right=.96,top=.91,bottom=.13,wspace=.30,hspace=.34)
for i,(ax,db) in enumerate(zip(axes.flat,NAMES)):
 z=density[density.database.eq(db)&density.n_pairs.ge(5)]
 xbins=np.arange(20,202,2);ybins=np.arange(-180,182,2);grid=np.full((90,180),np.nan)
 for _,r in z.iterrows():
  xi=int(round((r.mean_map-21)/2));yi=int(round((r.difference+179)/2))
  if 0<=xi<90 and 0<=yi<180:
   grid[xi,yi]=r.n_pairs
   if 29<=r.mean_map<=141 and -71<=r.difference<=71:
    tr('FigureS1',db+': density',r,['database','mean_map','difference'],['n_pairs'])
 ck('density_count_eligible:'+db,np.nanmin(grid)>=5)
 image=ax.pcolormesh(xbins,ybins,grid.T,cmap='cividis',norm=LogNorm(vmin=5,vmax=max(6,z.n_pairs.max())),rasterized=True)
 q=ba[ba.database.eq(db)&ba.stratification.eq('overall')]
 ck('ba_one_overall:'+db,len(q)==1);r=q.iloc[0]
 ck('ba_loa_order:'+db,r.loa_lower<r.bias<r.loa_upper)
 ck('ba_loa_formula:'+db,np.isclose(r.loa_upper-r.bias,1.96*np.sqrt(r.patient_variance+r.residual_variance)))
 tr('FigureS1',db+': bias and LoA',r,['database','stratification','stratum'],['n_pairs','n_patients','bias','loa_lower','loa_upper'])
 ax.axhline(r.bias,color='#B45A3C',lw=1.15)
 for bound in [r.loa_lower,r.loa_upper]:ax.axhline(bound,color='#B45A3C',ls='--',lw=.85)
 ax.set(xlim=(30,140),ylim=(-70,70))
 ax.set_title(NAMES[db],fontsize=8,fontweight='bold',loc='left',pad=9)
 ax.text(-.18,1.05,chr(65+i),transform=ax.transAxes,fontweight='bold',fontsize=9,va='bottom')
 ax.text(.03,.96,f'Bias {r.bias:.1f}; LoA {r.loa_lower:.1f} to {r.loa_upper:.1f} mmHg',
         transform=ax.transAxes,va='top',fontsize=5.9,bbox=dict(facecolor='white',edgecolor='none',alpha=.88))
 cb=fig.colorbar(image,ax=ax,shrink=.76,pad=.025);cb.ax.tick_params(labelsize=6);cb.set_label('Pairs per bin',fontsize=6.2)
 for lab in ax.get_xticklabels()+ax.get_yticklabels():lab.set_fontsize(6.7)
for ax in axes[:,0]:ax.set_ylabel('Cuff − arterial reference MAP (mmHg)',fontsize=7)
for ax in axes[1,:]:ax.set_xlabel('Mean of cuff and arterial reference MAP (mmHg)',fontsize=6.8)
fig.text(.09,.97,'Repeated-measures Bland–Altman analysis',fontsize=9,fontweight='bold')
fig.text(.09,.035,'Solid line: bias. Dashed lines: 95% limits of agreement. Density cells contain ≥5 pairs.',fontsize=6.5)
save(fig,'FigureS1_bland_altman','Supplementary Figure 1. Repeated-measures Bland–Altman analysis. The horizontal solid line is the estimated cuff-minus-arterial-reference MAP bias; dashed lines are 95% population limits of agreement from a REML patient random-intercept model, calculated as bias ±1.96 × the square root of the sum of patient and residual variances. Limits of agreement describe the distribution of differences and are not confidence intervals. Colour indicates counts in 2-mmHg density bins; bins with fewer than 10 pairs or fewer than 10 distinct patients are not displayed. Displayed axes crop the density to improve readability, while all eligible pairs contribute to the estimates. The colour scale is logarithmic and database-specific.')

curve=load('phase1_secondary_model_curve.csv');points=load('phase1_secondary_model_predictions.csv')
z=curve[curve.nee.le(2.5)]
ck('curve_domain',not z.outside_observed_range.any())
ck('curve_bounds',((z.ci_lower>=0)&(z.ci_lower<=z.estimate)&(z.estimate<=z.ci_upper)&(z.ci_upper<=1)).all())
for _,r in z.iterrows():tr('FigureS2','dose curve',r,['nee'],['estimate','ci_lower','ci_upper'])
for _,r in points.iterrows():tr('FigureS2','reported dose point',r,['nee'],['estimate','ci_lower','ci_upper'])
fig,ax=plt.subplots(figsize=(7.2,4.2));fig.subplots_adjust(left=.105,right=.98,top=.87,bottom=.19)
ax.fill_between(z.nee,z.ci_lower*100,z.ci_upper*100,color='#285C7A',alpha=.18,label='95% confidence interval')
ax.plot(z.nee,z.estimate*100,color='#285C7A',lw=1.7,label='Adjusted missed hypotension rate')
ax.scatter(points.nee,points.estimate*100,color='#285C7A',s=15,zorder=4)
ax.set(xlim=(0,2.5),ylim=(0,65),xlabel='NE-eq (µg/kg/min)',ylabel='Adjusted missed hypotension rate (%)')
ck('curve_unclipped',z.ci_upper.max()*100<65)
ax.grid(axis='y',color='#DDDDDD',lw=.55);ax.legend(loc='lower left',fontsize=7.2)
ax.set_title('Dose association after adjustment',fontweight='bold',fontsize=9,loc='left',pad=12)
fig.text(.105,.055,'MIMIC-IV and eICU · patient random intercept · marginal standardization',fontsize=6.8)
save(fig,'FigureS2_dose_association','Supplementary Figure 2. Association between norepinephrine-equivalent dose and missed hypotension after covariate adjustment. The model uses complete-case pairs with arterial reference MAP <65 mmHg from MIMIC-IV and eICU, with a three-knot restricted cubic spline for dose, adjustment for arterial reference MAP, age, sex, BMI >40 kg/m², recorded atrial fibrillation and database, and a patient random intercept. The curve is marginalized over the observed covariate distribution and the fitted normal random-intercept distribution. The shaded band is a 95% delta-method confidence interval using the fixed-effect covariance with random-intercept variance treated as a plug-in estimate. The display spans 0–2.5 µg/kg/min; the full evaluated curve is available in the accompanying data. The curve describes measurement association and does not estimate the effect of a monitoring strategy.')

pd.DataFrame(sources).to_csv(OUT.parent/'supplement_figure_sources.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(checks).to_csv(OUT.parent/'supplement_figure_validation.csv',index=False,encoding='utf-8-sig')
print(json.dumps({'figures':2,'sources':len(sources),'checks':len(checks),'passed':all(x['passed'] for x in checks)}))
