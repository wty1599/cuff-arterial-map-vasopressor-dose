#!/usr/bin/env python3
"""Rebuild the S6 disclosure layer without refitting Bland--Altman estimates.

The original 2 by 2 mmHg bins are retained. A bin is exported and displayed
only when it contains >=10 pairs AND >=10 distinct patients. All fields of
other bins are omitted, rather than releasing their coordinates or counts.
Full-data bias and limits of agreement are reused byte-for-number from the
existing overall summary. The unsuppressed comparison is restricted-only.
No rows, IDs, timestamps or small-cell values are written to stdout.
"""
from pathlib import Path
import argparse,ast,csv,hashlib,json,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

DATABASES={'mimic':'MIMIC-IV','eicu':'eICU','amsterdam':'AmsterdamUMCdb','sicdb':'SICdb'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def plot(data,refs,edges,title):
 fig,axes=plt.subplots(2,2,figsize=(7.2,5.8),constrained_layout=True)
 plt.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42,'font.family':'Arial','font.size':8})
 for ax,(key,db) in zip(axes.flat,DATABASES.items()):
  g=data.loc[data.database.eq(db)];h=np.zeros((len(edges[0])-1,len(edges[1])-1))
  for r in g.itertuples(index=False):
   ix=int(round((r.mean_map-edges[0][0]-1)/2));iy=int(round((r.difference-edges[1][0]-1)/2));h[ix,iy]=r.n_pairs
  # Original published display already omitted <5-pair bins. The remediated
  # input has a stronger pair AND patient threshold before reaching this step.
  minimum=10 if 'n_patients' in data.columns else 5
  h=np.ma.masked_less(h,minimum)
  vmax=max(10,float(h.max()))
  im=ax.pcolormesh(edges[0],edges[1],h.T,norm=LogNorm(vmin=minimum,vmax=vmax),cmap='viridis',rasterized=True,shading='flat')
  rr=refs.loc[refs.database.eq(db)].iloc[0]
  ax.axhline(rr.bias,color='#303030',linewidth=1)
  ax.axhline(rr.loa_lower,color='#606060',linestyle='--',linewidth=.8)
  ax.axhline(rr.loa_upper,color='#606060',linestyle='--',linewidth=.8)
  ax.set(xlim=(30,130),ylim=(-60,60),title=db,xlabel='Mean of cuff and arterial MAP (mmHg)',ylabel='Cuff minus arterial MAP (mmHg)')
  ax.spines[['top','right']].set_visible(False)
  fig.colorbar(im,ax=ax,label='Pairs per displayed bin',fraction=.045,pad=.025)
 fig.suptitle(title,fontsize=11)
 return fig
def main():
 a=argparse.ArgumentParser();a.add_argument('--project',required=True);a.add_argument('--original',required=True);a.add_argument('--output',required=True);a.add_argument('--restricted-qa',required=True)
 x=a.parse_args();root=Path(x.project).resolve();out=Path(x.output).resolve();qa=Path(x.restricted_qa).resolve()
 if not any(v in {'restricted','_restricted'} for v in qa.parts):raise RuntimeError('comparison_requires_restricted_directory')
 out.mkdir(parents=True,exist_ok=True);qa.mkdir(parents=True,exist_ok=True)
 original=Path(x.original).resolve();old=pd.read_csv(original)
 # Literal bin edges are read from the original producer; no patient content
 # is evaluated. This also prevents silent changes in grid origin or width.
 src=root/'sql/phase1/02_rates.py';tree=ast.parse(src.read_text(encoding='utf-8-sig'))
 hist=[n for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='histogram2d'][0]
 bs=next(k.value for k in hist.keywords if k.arg=='bins');calls=bs.elts
 edges=[np.arange(*(ast.literal_eval(v) for v in c.args)) for c in calls]
 density=[];summary=[];sources=[]
 for key,db in DATABASES.items():
  f=root/'phase1/_restricted'/('pairs_'+key+'.parquet')
  g=pd.read_parquet(f,columns=['patient','cuff_map','art_map'])
  mean=(g.cuff_map.to_numpy(float)+g.art_map.to_numpy(float))/2
  delta=g.cuff_map.to_numpy(float)-g.art_map.to_numpy(float)
  ix=np.searchsorted(edges[0],mean,side='right')-1;iy=np.searchsorted(edges[1],delta,side='right')-1
  ix[mean==edges[0][-1]]=len(edges[0])-2;iy[delta==edges[1][-1]]=len(edges[1])-2
  ok=np.isfinite(mean)&np.isfinite(delta)&(ix>=0)&(iy>=0)&(ix<len(edges[0])-1)&(iy<len(edges[1])-1)
  cells=pd.DataFrame({'ix':ix[ok],'iy':iy[ok],'patient':g.patient.to_numpy()[ok]})
  full=cells.groupby(['ix','iy']).agg(n_pairs=('patient','size'),n_patients=('patient','nunique')).reset_index()
  full['database']=db;full['mean_map']=edges[0][full.ix.to_numpy()]+1;full['difference']=edges[1][full.iy.to_numpy()]+1
  fields=['database','mean_map','difference','n_pairs','n_patients'];full=full[fields]
  o=old.loc[old.database.eq(db),['database','mean_map','difference','n_pairs']].sort_values(['mean_map','difference']).reset_index(drop=True)
  n=full[['database','mean_map','difference','n_pairs']].sort_values(['mean_map','difference']).reset_index(drop=True)
  if not np.array_equal(o.to_numpy(),n.to_numpy()):
   (qa/'density_reproduction_failure.json').write_text(json.dumps({'database':db,'status':'FAIL','reason':'existing_density_not_reproduced','source':str(f)},ensure_ascii=False),encoding='utf-8')
   raise RuntimeError('existing_density_not_reproduced_see_restricted_evidence')
  keep=(full.n_pairs>=10)&(full.n_patients>=10)
  safe=full.loc[keep].copy();density.append(safe)
  # Raw bin membership remains restricted and is not printed or delivered.
  # Only pass/fail and source hashes are retained for suppressed bins; a
  # second unsuppressed bin-level extract is not kept after regeneration.
  summary.append({'database':db,'original_bins':len(full),'retained_bins':int(keep.sum()),'suppressed_bins':int((~keep).sum()),'original_pairs_in_grid':int(full.n_pairs.sum()),'retained_pairs_in_grid':int(safe.n_pairs.sum()),'display_threshold':'at least 10 pairs and at least 10 distinct patients'})
  sources.append({'database':db,'source':str(f),'sha256':sha(f),'existing_density_reproduced':True})
 safe=pd.concat(density,ignore_index=True)
 safe.to_csv(out/'S6_bland_altman_density.csv',index=False)
 ba=pd.read_csv(root/'phase1/phase1_secondary_bland_altman.csv');ba['database']=ba['database'].replace(DATABASES);ref=[]
 for db in DATABASES.values():
  g=ba.loc[ba.database.eq(db)&ba.stratification.eq('overall')]
  # Match the original published plot exactly: it selects the first overall
  # row for each database. The summary also contains other overall variants;
  # none is refitted or substituted during this disclosure-only change.
  if g.empty:raise RuntimeError('missing_original_overall_reference_summary')
  ref.append(g.iloc[0])
 refs=pd.DataFrame(ref).drop(columns=['row_id'],errors='ignore')
 refs.to_csv(out/'S6_reference_lines.csv',index=False)
 pd.DataFrame(summary).to_csv(out/'S6_suppression_summary.csv',index=False)
 plt.rcParams.update({'font.family':'Arial','font.size':8,'svg.fonttype':'none','pdf.fonttype':42})
 fig=plot(safe,refs,edges,'Repeated-measures Bland–Altman analysis')
 for ext in ['png','svg','pdf']:fig.savefig(out/('FigureS1_bland_altman.'+ext),dpi=600,bbox_inches='tight')
 plt.close(fig)
 before=plot(old,refs,edges,'Before remediation · restricted review only')
 before.savefig(qa/'S6_before.png',dpi=180,bbox_inches='tight');plt.close(before)
 # Build the requested side-by-side review using rasterized images. The
 # unsuppressed-side comparison never leaves this restricted directory.
 from PIL import Image
 bef=Image.open(qa/'S6_before.png').convert('RGB');aft=Image.open(out/'FigureS1_bland_altman.png').convert('RGB')
 aft=aft.resize(bef.size)
 comp=Image.new('RGB',(bef.width*2,bef.height),'white');comp.paste(bef,(0,0));comp.paste(aft,(bef.width,0));comp.save(qa/'S6_before_after_restricted.png')
 receipt={'status':'COMPLETE','no_estimates_refitted':True,'reference_source':str(root/'phase1/phase1_secondary_bland_altman.csv'),'reference_source_sha256':sha(root/'phase1/phase1_secondary_bland_altman.csv'),'original_density_sha256':sha(original),'new_density_sha256':sha(out/'S6_bland_altman_density.csv'),'suppression_rule':'n_pairs >=10 AND distinct n_patients >=10','source_pairs':sources,'summary':summary,'restricted_comparison':str(qa/'S6_before_after_restricted.png'),'versions':{'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__}}
 (out/'S6_regeneration_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
 print(json.dumps({'status':'COMPLETE','safe_density_rows':len(safe),'output':str(out),'comparison':'written to restricted directory','raw_rows_or_cell_values_printed':False},ensure_ascii=False))
if __name__=='__main__':main()
