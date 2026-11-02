"""Regenerate the v5 S6 display layer using credentialed local pair files.

No clinical estimates are fitted. The supplied reference CSV must contain the
existing full-data bias and limits. Entire bins below either disclosure limit
are omitted before coordinates or counts are exported. Run only under approved
source-data access; the repository supplies no patient-level input files.
"""
from pathlib import Path
import argparse,importlib.util,json,sys,csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
ROOT=Path(__file__).resolve().parents[2]
DB={'mimic':'MIMIC-IV','eicu':'eICU','amsterdam':'AmsterdamUMCdb','sicdb':'SICdb'}

def coverage(s,w,lo,hi):
 visible=s.loc[(s.mean_map+w/2>30)&(s.mean_map-w/2<130)]
 intervals=sorted((max(lo,-60,float(y)-w/2),min(hi,60,float(y)+w/2)) for y in visible.difference.unique())
 end=-np.inf;length=0.
 for a,b in intervals:
  if b>a:length+=max(0.,b-max(a,end));end=max(end,b)
 return length/(hi-lo) if hi>lo else 0.

def gate(paths):
 spec=importlib.util.spec_from_file_location('v5_gate',ROOT/'tools/small_cell_gate.py')
 mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
 hits,_,_=mod.scan_paths(paths,include_restricted=True)
 if hits:raise RuntimeError('Disclosure gate failed; no figure exported')

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--pairs-dir',type=Path,required=True)
 ap.add_argument('--reference-csv',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True)
 a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
 density=[];summary=[]
 for key,db in DB.items():
  g=pd.read_parquet(a.pairs_dir/f'pairs_{key}.parquet',columns=['patient','cuff_map','art_map'])
  if g.patient.nunique()<10:raise RuntimeError('Database cohort is below the disclosure threshold')
  mean=(g.cuff_map.to_numpy(float)+g.art_map.to_numpy(float))/2
  delta=g.cuff_map.to_numpy(float)-g.art_map.to_numpy(float)
  lo,hi=np.quantile(delta[np.isfinite(delta)],[.025,.975]);chosen=None
  for w in ([2] if key in ['mimic','eicu'] else [5,10]):
   ex=np.arange(20,200+w,w);ey=np.arange(-180,180+w,w)
   ix=np.searchsorted(ex,mean,side='right')-1;iy=np.searchsorted(ey,delta,side='right')-1
   ix[mean==ex[-1]]=len(ex)-2;iy[delta==ey[-1]]=len(ey)-2
   ok=np.isfinite(mean)&np.isfinite(delta)&(ix>=0)&(iy>=0)&(ix<len(ex)-1)&(iy<len(ey)-1)
   bins=pd.DataFrame({'ix':ix[ok],'iy':iy[ok],'patient':g.patient.to_numpy()[ok]})
   bins=bins.groupby(['ix','iy']).agg(n_pairs=('patient','size'),n_patients=('patient','nunique')).reset_index()
   # Omit the complete row before creating any exportable coordinate table.
   keep=bins.loc[(bins.n_pairs>=10)&(bins.n_patients>=10)].copy()
   keep['database']=db;keep['mean_map']=ex[keep.ix.to_numpy()]+w/2;keep['difference']=ey[keep.iy.to_numpy()]+w/2;keep['grid_width_mmHg']=w
   keep=keep[['database','mean_map','difference','grid_width_mmHg','n_pairs','n_patients']]
   cov=coverage(keep,w,lo,hi)
   if key in ['mimic','eicu'] or cov>=.8:chosen=keep;break
  mode='density' if chosen is not None else 'bias and limits only'
  if chosen is not None:density.append(chosen)
  summary.append({'database':db,'grid_width_mmHg':w,'coverage_fraction':cov,'central95_lower_mmHg':lo,'central95_upper_mmHg':hi,'display_mode':mode})
 safe=pd.concat(density,ignore_index=True);safe.to_csv(out/'S6_bland_altman_density.csv',index=False)
 pd.DataFrame(summary).to_csv(out/'S6_coverage.csv',index=False)
 refs=pd.read_csv(a.reference_csv)
 required=['database','bias','loa_lower','loa_upper']
 if any(k not in refs.columns for k in required):raise ValueError('Missing reference columns')
 # Deliberately export only the already-aggregated reference lines used here.
 with a.reference_csv.open(encoding='utf-8-sig',newline='') as f:original_refs=list(csv.DictReader(f))
 with (out/'S6_reference_lines.csv').open('w',encoding='utf-8',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=required);writer.writeheader();writer.writerows({k:r[k] for k in required} for r in original_refs)
 gate([out])
 plt.rcParams.update({'font.family':'Arial','font.size':7.5,'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.7})
 fig,axs=plt.subplots(2,2,figsize=(170/25.4,145/25.4),layout='constrained')
 for ax,rec in zip(axs.flat,summary):
  db=rec['database'];w=rec['grid_width_mmHg'];s=safe.loc[safe.database.eq(db)]
  if rec['display_mode']=='density':
   ex=np.arange(20,200+w,w);ey=np.arange(-180,180+w,w);h=np.ma.masked_all((len(ey)-1,len(ex)-1))
   for r in s.itertuples(index=False):h[int(round((r.difference+180-w/2)/w)),int(round((r.mean_map-20-w/2)/w))]=r.n_pairs
   im=ax.pcolormesh(ex,ey,h,norm=LogNorm(vmin=10,vmax=max(10,float(h.max()))),cmap='viridis',rasterized=True,shading='flat')
   fig.colorbar(im,ax=ax,label='Pairs per displayed bin',fraction=.045,pad=.025)
  else:ax.text(.5,.52,'density not shown: insufficient cell sizes',transform=ax.transAxes,ha='center',va='center',fontsize=7,wrap=True)
  ref=refs.loc[refs.database.eq(db)].iloc[0]
  ax.axhline(ref.bias,color='#303030',linewidth=1)
  for y in [ref.loa_lower,ref.loa_upper]:ax.axhline(y,color='#606060',linestyle='--',linewidth=.8)
  ax.set(xlim=(30,130),ylim=(-60,60),title=db,xlabel='Mean of cuff and arterial MAP (mmHg)',ylabel='Cuff minus arterial MAP (mmHg)')
  ax.spines[['top','right']].set_visible(False)
 for ext in ['png','svg','pdf']:fig.savefig(out/f'FigureS1_bland_altman.{ext}',dpi=600)
 plt.close(fig);gate([out]);print(json.dumps({'status':'PASS','output':str(out),'display_rows':len(safe),'estimates_refitted':False}))
if __name__=='__main__':main()
