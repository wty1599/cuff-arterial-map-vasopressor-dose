"""Rebuild Brief Report rate tables and Fig. 1 from disclosure-controlled CSVs.

No access to raw data, patient records, or a model fit is required.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import html
import json
import platform

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'outputs/aggregate'
STRATA=['<=0.1','>0.1-0.25','>0.25-0.5','>0.5-1.0','>1.0']
LABELS=['≤0.1','>0.1–0.25','>0.25–0.5','>0.5–1.0','>1.0']
DATABASES=['MIMIC-IV','eICU','AmsterdamUMCdb','SICdb']

def read(name):
    with (DATA/name).open(encoding='utf-8',newline='') as handle:
        return list(csv.DictReader(handle))

def select(rows,**where):
    selected=[r for r in rows if all(r.get(k)==v for k,v in where.items())]
    if len(selected)!=1:
        raise ValueError(f'Expected one aggregate row for {where}; found {len(selected)}')
    return selected[0]

def pct(row,key='estimate',ci=True):
    if row.get('display')=='Not pooled':return 'Not pooled'
    value=f'{100*float(row[key]):.1f}'
    if ci and row.get('ci_lower') and row.get('ci_upper'):
        value+=f" ({100*float(row['ci_lower']):.1f}–{100*float(row['ci_upper']):.1f})"
    if 'insufficient' in row.get('precision',''):value+='*'
    return value

def save_table(out,name,rows,caption):
    fields=list(rows[0])
    with (out/f'{name}.csv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    text=f'# {caption}\n\n'
    text+='| '+' | '.join(fields)+' |\n| '+' | '.join('---' for _ in fields)+' |\n'
    text+=''.join('| '+' | '.join(str(r[k]) for k in fields)+' |\n' for r in rows)
    text+='\nValues are percentages (95% confidence intervals). Counts are withheld. * Operational precision flag; this can also be triggered by an estimate less than one percentage point from a narrow CI bound.\n'
    if name=='Table1':text+='Pooled estimates use original known-dose observations, logit REML and modified Hartung–Knapp intervals. Pooled false-alarm estimates are not interpretable because of between-database heterogeneity (I² 98.9%).\n'
    (out/f'{name}.md').write_text(text,encoding='utf-8')
    doc='<html><meta charset="utf-8"><style>body{font-family:Arial;margin:2em}table{border-collapse:collapse;font-size:13px}td,th{border:1px solid #aaa;padding:7px}th{background:#eee}</style>'
    doc+='<h1>'+html.escape(caption)+'</h1><table><tr>'+''.join('<th>'+html.escape(k)+'</th>' for k in fields)+'</tr>'
    doc+=''.join('<tr>'+''.join('<td>'+html.escape(str(r[k]))+'</td>' for k in fields)+'</tr>' for r in rows)+'</table><p>Rate components only; counts withheld. See the aggregate disclosure policy.</p></html>'
    (out/f'{name}.html').write_text(doc,encoding='utf-8')

def tables(out):
    rates=read('S2_missed_and_false_alarm.csv')
    deep=read('S2_deep_missed_hypotension.csv')
    cuff=read('S3_cuff_conditioned_probabilities.csv')
    standard=read('S5_standardized_missed_hypotension.csv')
    pooled=read('S2_pooled_missed_and_false_alarm.csv')
    pooled_deep=read('S2_pooled_deep_missed_hypotension.csv')
    table1=[]
    for db in DATABASES:
        missed=select(rates,database=db,stratum='all',endpoint='Missed hypotension rate')
        dep=select(deep,database=db,stratum='all',metric='Deep missed hypotension rate (arterial reference MAP <55 mmHg)')
        alarm=select(rates,database=db,stratum='all',endpoint='False-alarm rate')
        # Layer prevalence in deep-miss rows is the full same-stratum prevalence.
        table1.append({'Database':db,'Arterial MAP <65 (%)':pct(missed,'arterial_hypotension_prevalence',False),
          'Missed hypotension (%)':pct(missed),'Deep missed hypotension (%)':pct(dep),'False alarm (%)':pct(alarm)})
    table1.append({'Database':'Pooled (known dose)','Arterial MAP <65 (%)':'Not pooled',
      'Missed hypotension (%)':pct(select(pooled,metric='Missed hypotension rate',stratum='known')),
      'Deep missed hypotension (%)':pct(select(pooled_deep,metric='Deep missed hypotension rate (arterial reference MAP <55 mmHg)',stratum='known')),
      'False alarm (%)':pct(select(pooled,metric='False-alarm rate',stratum='known'))})
    table2=[]
    for db in DATABASES[:2]:
        for st in STRATA:
            prevalence=select(cuff,database=db,stratum=st,category='all_pairs')
            conditional=select(cuff,database=db,stratum=st,category='all_normal_cuff')
            near=select(cuff,database=db,stratum=st,category='65-69',metric='Probability of arterial reference MAP <65 within cuff MAP category')
            table2.append({'Database':db,'NE-eq (µg/kg/min)':st,
              'Arterial MAP <65 (%)':pct(prevalence,ci=False),
              'Missed hypotension (%)':pct(select(rates,database=db,stratum=st,endpoint='Missed hypotension rate')),
              'Depth-standardized missed (%)':pct(select(standard,database=db,stratum=st,scenario='Primary five-minute pairing')),
              'False alarm (%)':pct(select(rates,database=db,stratum=st,endpoint='False-alarm rate')),
              'Arterial <65 given cuff ≥65 (%)':pct(conditional),
              'Arterial <65 given cuff 65–69 (%)':pct(near)})
    save_table(out,'Table1',table1,'Brief Report Table 1: overall threshold disagreement')
    save_table(out,'Table2',table2,'Brief Report Table 2: dose-stratified threshold probabilities')
    return {'Table1_rows':len(table1),'Table2_rows':len(table2)}

def plot(out):
    expected=read('phase5_S5_counterfactual.csv');cuff=read('S3_cuff_conditioned_probabilities.csv')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'svg.fonttype':'none','svg.hashsalt':'cuff-arterial-map',
                        'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})
    fig,axes=plt.subplots(2,2,figsize=(8.2,7.4))
    points=[]
    for col,db in enumerate(DATABASES[:2]):
        ax=axes[0,col];rr=[select(expected,database=db,stratum=st) for st in STRATA]
        for key,label,color,marker,offset in [('observed','Observed','#0072B2','o',-.08),('expected','Common-relation expected','#D55E00','s',.08)]:
            y=np.array([float(r[key])*100 for r in rr]);lo=np.array([float(r[key+'_ci_lower'])*100 for r in rr]);hi=np.array([float(r[key+'_ci_upper'])*100 for r in rr])
            ax.errorbar(np.arange(5)+offset,y,yerr=[y-lo,hi-y],fmt=marker+'-',color=color,mfc='white' if key=='expected' else color,markersize=4,linewidth=1,capsize=2,label=label)
            for st,r in zip(STRATA,rr):points.append({'panel':chr(65+col),'database':db,'series':key,'stratum':st,'category':'','value':r[key],'ci_lower':r[key+'_ci_lower'],'ci_upper':r[key+'_ci_upper'],'source_csv':r['source_csv'],'source_line':r['source_line']})
        prev=[float(r['arterial_hypotension_prevalence'])*100 for r in rr]
        ax.plot(range(5),prev,'D:',color='#888888',markersize=3.5,linewidth=.9,label='Arterial hypotension prevalence')
        ax.set_xticks(range(5),LABELS,fontsize=8);ax.set_xlabel('NE-eq (µg/kg/min)');ax.set_title(db,pad=10);ax.set_ylim(0,65)
        for st,r in zip(STRATA,rr):points.append({'panel':chr(65+col),'database':db,'series':'prevalence','stratum':st,'category':'','value':r['arterial_hypotension_prevalence'],'ci_lower':'','ci_upper':'','source_csv':r['source_csv'],'source_line':r['source_line']})
        ax=axes[1,col]
        for st,label,color,marker in [('<=0.25','≤0.25','#0072B2','o'),('>0.25-1.0','>0.25–1.0','#D55E00','s'),('>1.0','>1.0','#009E73','^')]:
            bins=['65-69','70-74','75-79','>=80']
            rr=[select(cuff,database=db,stratum='dosegroup:'+st,category=b,metric='Probability of arterial reference MAP <65 within cuff MAP category') for b in bins]
            y=np.array([float(r['estimate'])*100 for r in rr]);lo=np.array([float(r['ci_lower'])*100 for r in rr]);hi=np.array([float(r['ci_upper'])*100 for r in rr])
            ax.errorbar(range(4),y,yerr=[y-lo,hi-y],fmt=marker+'-',color=color,markersize=4,linewidth=1,capsize=2,label=label)
            for b,r in zip(bins,rr):points.append({'panel':chr(67+col),'database':db,'series':st,'stratum':st,'category':b,'value':r['estimate'],'ci_lower':r['ci_lower'],'ci_upper':r['ci_upper'],'source_csv':r['source_csv'],'source_line':r['source_line']})
        ax.set_xticks(range(4),['65–69','70–74','75–79','≥80']);ax.set_xlabel('Cuff MAP (mmHg)');ax.set_ylim(0,75)
    for row in range(2):
        for col in range(2):
            ax=axes[row,col];ax.text(-.12,1.04,chr(65+row*2+col),transform=ax.transAxes,fontweight='bold',fontsize=12)
            ax.yaxis.grid(True,color='#E8E8E8',linewidth=.5);ax.set_axisbelow(True)
            if col==0:ax.set_ylabel('Arterial MAP <65 mmHg (%)')
    handles,labels=axes[0,0].get_legend_handles_labels();d=dict(zip(labels,handles))
    order=['Observed','Common-relation expected','Arterial hypotension prevalence']
    fig.legend([d[x] for x in order],order,loc='upper center',bbox_to_anchor=(.54,.535),ncol=3,fontsize=7.5)
    fig.legend(*axes[1,0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.54,.025),ncol=3,title='NE-eq (µg/kg/min)',fontsize=8)
    fig.subplots_adjust(left=.10,right=.98,top=.94,bottom=.16,wspace=.23,hspace=.66)
    for ext in ['png','svg','pdf']:
        metadata={'Date':None} if ext=='svg' else ({'CreationDate':None,'ModDate':None} if ext=='pdf' else {})
        fig.savefig(out/f'Figure1.{ext}',dpi=600,metadata=metadata)
    plt.close(fig)
    with (out/'Figure1_source_trace.csv').open('w',encoding='utf-8',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(points[0]));w.writeheader();w.writerows(points)
    caption='Figure 1. Residual arterial hypotension across dose and cuff readings. Upper panels: observed (filled) and common-relation expected (open) probabilities given cuff MAP ≥65 mmHg; grey diamonds show arterial hypotension prevalence. Lower panels: probabilities by cuff interval and dose group. Bars show original 95% patient-cluster intervals, redrawn from rounded aggregates. Counts are withheld in this repository. NE-eq, norepinephrine-equivalent dose; MAP, mean arterial pressure.'
    (out/'Figure1_caption.txt').write_text(caption+'\n',encoding='utf-8')
    return {'Figure1_points':len(points),'Figure1_caption_words':len(caption.split())}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=ROOT/'reproduced');args=parser.parse_args()
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    receipt={'python':platform.python_version(),'matplotlib':matplotlib.__version__,'numpy':np.__version__,'input_kind':'rounded disclosure-controlled aggregates; no refitting'}
    receipt.update(tables(out));receipt.update(plot(out))
    receipt['source_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(DATA.glob('*.csv'))}
    (out/'reproduction_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in receipt.items() if k!='source_sha256'},indent=2))

if __name__=='__main__':main()
