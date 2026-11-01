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
import sys, json, hashlib, math
ROOT = Path('@CUFF_WORKSPACE@')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
P1, P2, OUT = ROOT/'phase1', ROOT/'phase2', ROOT/'phase3'
FIG, TAB = OUT/'figures', OUT/'tables'
for d in [OUT, FIG, TAB]: d.mkdir(parents=True, exist_ok=True)
DBS = ['mimic','eicu','amsterdam','sicdb']
MAIN = ['mimic','eicu']
NAMES = {'mimic':'MIMIC-IV','eicu':'eICU','amsterdam':'AmsterdamUMCdb','sicdb':'SICdb','pooled':'Random-effects pooled'}
STRATA = ['<=0.1','>0.1-0.25','>0.25-0.5','>0.5-1.0','>1.0']
TICKS = ['≤0.1','>0.1–0.25','>0.25–0.5','>0.5–1.0','>1.0']
COL = {'mimic':'#285C7A','eicu':'#B45A3C','pooled':'#222222'}
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
    'font.size':7.5,'axes.spines.top':False,'axes.spines.right':False,
    'axes.linewidth':.65,'legend.frameon':False,'svg.fonttype':'none','pdf.fonttype':42,
    'axes.labelsize':7.5,'xtick.labelsize':6.5,'ytick.labelsize':7,
    'savefig.facecolor':'white'})
SOURCE_ROWS, CHECKS, EXPORTS = [], [], []

def read(rel):
    df = pd.read_csv(ROOT/rel)
    df['_source'] = rel
    df['_row'] = np.arange(len(df))+2
    return df

def check(name, ok, detail=''):
    CHECKS.append(dict(check=name, passed=bool(ok), detail=str(detail)))
    if not ok: raise AssertionError(f'{name}: {detail}')

def clean(x):
    if isinstance(x, (np.integer,)): return int(x)
    if isinstance(x, (np.floating,)): return None if pd.isna(x) else float(x)
    if pd.isna(x): return None
    return x

def trace(artifact, location, r, keys, fields, derived=''):
    SOURCE_ROWS.append(dict(artifact=artifact, location=location, source_csv=r['_source'],
       source_csv_row=int(r['_row']), source_key=json.dumps({k:clean(r[k]) for k in keys},ensure_ascii=False,sort_keys=True),
       source_fields=';'.join(fields), source_values=json.dumps({k:clean(r[k]) for k in fields},ensure_ascii=False,sort_keys=True),
       display_transformation=derived))

def one(df, **filters):
    q = df
    for k,v in filters.items(): q = q[q[k].eq(v)]
    check('unique_source_row:'+','.join(f'{k}={v}' for k,v in filters.items()), len(q)==1, len(q))
    return q.iloc[0]

def save(fig, name):
    fig.canvas.draw()
    for ext in ('png','svg','pdf'):
        fig.savefig(FIG/f'{name}.{ext}', dpi=400, facecolor='white')
    EXPORTS.append(dict(figure=name,width_inches=fig.get_figwidth(),height_inches=fig.get_figheight(),
                       dpi=400,backend='Python/Matplotlib',editable_svg_text=True,formats='png;svg;pdf'))
    plt.close(fig)

PRIMARY = {db:read(f'phase1/phase1_primary_{db}.csv') for db in DBS}
PRIMARY_POOL = read('phase1/phase1_pooled.csv')
DEEP = read('phase2/phase2_A_deep_miss.csv')
DEEP_POOL = read('phase2/phase2_A_deep_pooled.csv')
B = read('phase2/phase2_B_cuff_perspective.csv')
BASELINE = read('phase2/phase2_table1.csv')
FLOW = read('phase2/phase2_flow_numbers.csv')

def primary_row(db, st):
    return one(PRIMARY_POOL, endpoint='missed',stratum=st) if db=='pooled' else one(PRIMARY[db],endpoint='missed',stratum=st)

def deep_row(db, st):
    return one(DEEP_POOL,metric='deep_miss_art_lt55',stratum=st) if db=='pooled' else one(DEEP,database=db,metric='deep_miss_art_lt55',stratum=st)

def ratio_check(r, context, pooled=False):
    n,d = float(r.numerator),float(r.denominator)
    check(context+':counts', 0<=n<=d, f'{n}/{d}')
    if d>0:
        check(context+':probability_bounds', 0<=r.estimate<=1 and 0<=r.ci_lower<=r.ci_upper<=1)
        if not pooled: check(context+':ratio', np.isclose(r.estimate,n/d,atol=1e-12),f'{r.estimate} vs {n/d}')

def dose_axes(ax, ymax, xlabel=True, shade=False, wrap=False):
    labels=['≤0.1','>0.1–\n0.25','>0.25–\n0.5','>0.5–\n1.0','>1.0'] if wrap else TICKS
    ax.set(xticks=np.arange(5),xticklabels=labels,ylim=(0,ymax),xlim=(-.45,4.45))
    ax.grid(axis='y',color='#DDDDDD',lw=.55,zorder=0)
    if xlabel: ax.set_xlabel('NE-eq (µg/kg/min)', labelpad=6)
    ax.tick_params(axis='x',length=3,pad=3)
    if shade:
        ax.axvspan(2.5,4.45,color='#EDEDED',zorder=-1)
        annotation='EVERDAC\nunder-represented\nrange' if wrap else 'EVERDAC\nunder-represented range'
        ax.text(3.47,ymax-2.4,annotation,ha='center',va='top',fontsize=5.7,color='#555555',linespacing=1.0)

def draw_rates(ax,db,artifact,panel):
    for label,getter,shift,ls,marker in [('arterial reference MAP <65 mmHg',primary_row,-.055,'-','o'),
                                      ('arterial reference MAP <55 mmHg',deep_row,.055,'--','^')]:
        rows=[getter(db,st) for st in STRATA]
        for st,r in zip(STRATA,rows):
            ratio_check(r,f'{artifact}:{panel}:{label}:{st}',db=='pooled')
            keys=['endpoint','stratum'] if 'endpoint' in r.index else ['metric','stratum']
            if db!='pooled':keys=['database',*keys]
            fields=['numerator','denominator','estimate','ci_lower','ci_upper']+(['databases'] if db=='pooled' else [])
            trace(artifact,panel+': '+label+'; '+st,r,keys,fields, 'proportions ×100; 95% CI')
        q=pd.DataFrame(rows)
        check(f'{artifact}:{panel}:{label}:unclipped',q.ci_upper.max()*100<70)
        ax.errorbar(np.arange(5)+shift,q.estimate*100,
                    yerr=np.array([q.estimate-q.ci_lower,q.ci_upper-q.estimate])*100,
                    color=COL[db],ls=ls,marker=marker,markersize=3.5,linewidth=1.1,capsize=2,zorder=4)
    dose_axes(ax,70,shade=True,wrap=artifact=='Figure2')
    if db=='pooled':
        check(artifact+':pooled_high_dose_three',int(primary_row(db,'>1.0').databases)==3 and int(deep_row(db,'>1.0').databases)==3)
        ax.text(.04,.045,'>1.0: 3 databases',transform=ax.transAxes,fontsize=6.2,color='#333333')

def draw_b(ax,db,artifact,panel,prevalence_y=-.28,show_prevalence=True):
    rows=[one(B,database=db,metric='art_lt65_given_cuff_ge65',stratum=st) for st in STRATA]
    q=pd.DataFrame(rows)
    for st,r in zip(STRATA,rows):
        ratio_check(r,f'{artifact}:{panel}:B:{st}')
        check(f'{artifact}:{panel}:{st}:layer_prevalence',np.isclose(r.layer_low_arterial_prevalence,r.layer_low_arterial_pairs/r.layer_pairs))
        trace(artifact,panel+': '+db+'; '+st,r,['database','metric','stratum','category'],
              ['numerator','denominator','estimate','ci_lower','ci_upper','layer_low_arterial_pairs','layer_pairs','layer_low_arterial_prevalence'], 'probabilities and layer prevalence ×100')
    ax.errorbar(np.arange(5),q.estimate*100,
                yerr=np.array([q.estimate-q.ci_lower,q.ci_upper-q.estimate])*100,
                color=COL[db],ls='-',marker='o' if db=='mimic' else 's',markersize=3.8,
                linewidth=1.2,capsize=2,label=NAMES[db],zorder=4)
    check(f'{artifact}:{panel}:B:unclipped',q.ci_upper.max()*100<55)
    if show_prevalence:
        for x,r in enumerate(rows):
            ax.text(x,prevalence_y,f'{r.layer_low_arterial_prevalence*100:.1f}%',transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=6.8,color=COL[db],clip_on=False)
    return q

def panel_title(ax,letter,title,subtitle=None):
    ax.set_title(title,fontweight='bold',fontsize=8,pad=8,loc='center')
    ax.text(-.13,1.06,letter,transform=ax.transAxes,fontweight='bold',fontsize=10,ha='left',va='bottom')
    if subtitle: ax.text(.5,1.015,subtitle,transform=ax.transAxes,ha='center',va='bottom',fontsize=6.4,color='#555555')

def figure2():
    fig=plt.figure(figsize=(7.2,6.3))
    gs=fig.add_gridspec(2,6,left=.09,right=.985,top=.88,bottom=.17,hspace=.70,wspace=.95)
    for j,db in enumerate(['mimic','eicu','pooled']):
        ax=fig.add_subplot(gs[0,2*j:2*j+2]);draw_rates(ax,db,'Figure2',chr(65+j))
        panel_title(ax,chr(65+j),NAMES[db] if db!='pooled' else 'Pooled')
        if j==0:ax.set_ylabel('Cuff MAP ≥65 mmHg (%)')
        else:ax.tick_params(labelleft=False)
    for j,db in enumerate(MAIN):
        ax=fig.add_subplot(gs[1,3*j:3*j+3]);draw_b(ax,db,'Figure2',chr(68+j),prevalence_y=-.34)
        dose_axes(ax,55);panel_title(ax,chr(68+j),NAMES[db])
        ax.text(.5,1.00,'Among cuff MAP ≥65 mmHg',transform=ax.transAxes,ha='center',va='bottom',fontsize=6.7,color='#555555')
        if j==0:ax.set_ylabel('Arterial reference MAP <65 mmHg (%)')
        else:ax.tick_params(labelleft=False)
        ax.text(.5,-.26,'Arterial reference MAP <65 prevalence in each dose stratum',transform=ax.transAxes,ha='center',va='top',fontsize=6.0,color='#444444',clip_on=False)
    handles=[Line2D([0],[0],color='#333333',ls='-',marker='o',markersize=3.6,label='Arterial reference MAP <65 mmHg'),
             Line2D([0],[0],color='#333333',ls='--',marker='^',markersize=3.6,label='Arterial reference MAP <55 mmHg')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.53,.978),ncol=2,fontsize=7,handlelength=2.5,columnspacing=2)
    fig.text(.09,.025,'Points: pair-level estimates; bars: 95% CIs. Dose-stratified analyses use known NE-eq.',fontsize=6.5,color='#444444')
    save(fig,'Figure2_dose_and_cuff_probability')

def figure_icm():
    fig=plt.figure(figsize=(7.2,4.3))
    gs=fig.add_gridspec(1,2,left=.085,right=.985,top=.80,bottom=.29,wspace=.36)
    a=fig.add_subplot(gs[0,0]);draw_rates(a,'pooled','ICM_Figure1','A');panel_title(a,'A','Missed hypotension: pooled')
    a.set_ylabel('Cuff MAP ≥65 mmHg (%)')
    b=fig.add_subplot(gs[0,1]);dose_axes(b,55)
    for j,db in enumerate(MAIN):
        draw_b(b,db,'ICM_Figure1','B',show_prevalence=False)
    panel_title(b,'B','Arterial hypotension given cuff MAP ≥65')
    b.set_ylabel('Arterial reference MAP <65 mmHg (%)')
    b.legend(loc='upper left',fontsize=6.7)
    b.text(.5,-.24,'Arterial reference MAP <65 prevalence (%)',transform=b.transAxes,ha='center',va='top',fontsize=6.3)
    for j,db in enumerate(MAIN):
        y=-.34-.10*j
        b.text(-.15,y,NAMES[db],transform=b.transAxes,ha='right',va='top',fontsize=6,color=COL[db],clip_on=False)
        for x,st in enumerate(STRATA):
            r=one(B,database=db,metric='art_lt65_given_cuff_ge65',stratum=st)
            b.text(x,y,f'{r.layer_low_arterial_prevalence*100:.1f}',transform=b.get_xaxis_transform(),ha='center',va='top',fontsize=6.3,color=COL[db],clip_on=False)
    fig.legend(handles=[Line2D([0],[0],color='#333333',ls='-',marker='o',markersize=3.4,label='Arterial reference MAP <65 mmHg'),
                        Line2D([0],[0],color='#333333',ls='--',marker='^',markersize=3.4,label='Arterial reference MAP <55 mmHg')],
               loc='upper left',bbox_to_anchor=(.07,.98),fontsize=6.7,ncol=2,columnspacing=2)
    fig.text(.085,.015,'Points: pair-level estimates; bars: 95% CIs. Known NE-eq pairs; no pooling in panel B.',fontsize=6.3,color='#444444')
    save(fig,'ICM_Figure1_two_panel')

def figure3():
    fig,axes=plt.subplots(1,2,figsize=(7.2,4.7),sharey=True)
    fig.subplots_adjust(left=.09,right=.985,top=.80,bottom=.34,wspace=.20)
    groups=['dosegroup:<=0.25','dosegroup:>0.25-1.0','dosegroup:>1.0']
    labels=['NE-eq ≤0.25','NE-eq >0.25–1.0','NE-eq >1.0']
    colors=['#285C7A','#B45A3C','#222222'];markers=['o','s','^'];linestyles=['-','--',':']
    bins=['65-69','70-74','75-79','>=80']
    for j,(ax,db) in enumerate(zip(axes,MAIN)):
        for g,label,color,marker,ls in zip(groups,labels,colors,markers,linestyles):
            rows=[one(B,database=db,metric='art_lt65_given_cuff_bin',stratum=g,category=bn) for bn in bins]
            q=pd.DataFrame(rows)
            for bn,r in zip(bins,rows):
                ratio_check(r,f'Figure3:{db}:{g}:{bn}')
                trace('Figure3',NAMES[db]+': '+g+'; '+bn,r,['database','metric','stratum','category'],['numerator','denominator','estimate','ci_lower','ci_upper'],'proportions ×100; 95% CI')
            check('Figure3:'+db+':'+g+':unclipped',q.ci_upper.max()*100<70)
            ax.errorbar(np.arange(4),q.estimate*100,yerr=np.array([q.estimate-q.ci_lower,q.ci_upper-q.estimate])*100,
                        color=color,marker=marker,ls=ls,linewidth=1.2,capsize=2,markersize=3.7,label=label)
        ax.set(xticks=np.arange(4),xticklabels=['65–69','70–74','75–79','≥80'],xlabel='Cuff MAP (mmHg)',ylim=(0,70),xlim=(-.2,3.2))
        panel_title(ax,chr(65+j),NAMES[db]);ax.grid(axis='y',color='#DDDDDD',lw=.55)
        if j==0:ax.set_ylabel('Arterial reference MAP <65 mmHg (%)')
        ax.text(.5,-.25,'Arterial reference MAP <65 prevalence by dose group:',transform=ax.transAxes,fontsize=6.1,ha='center',va='top',clip_on=False)
        for n,(g,label,color) in enumerate(zip(groups,labels,colors)):
            r=one(B,database=db,metric='arterial_lt65_prevalence',stratum=g)
            trace('Figure3',NAMES[db]+': prevalence; '+g,r,['database','metric','stratum','category'],['numerator','denominator','estimate'],'proportion ×100')
            ax.text(.5,-.34-.085*n,label+f': {r.estimate*100:.1f}%',transform=ax.transAxes,fontsize=6.2,color=color,ha='center',va='top',clip_on=False)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.53,.975),ncol=3,fontsize=7,handlelength=2.5)
    fig.text(.09,.018,'Known-dose pairs. Points: conditional probabilities; bars: patient-cluster 95% CIs.',fontsize=6.3,color='#444444')
    save(fig,'Figure3_cuff_reading_bins')

def figure1():
    fig=plt.figure(figsize=(7.2,7.6));ax=fig.add_axes([.01,.02,.98,.96]);ax.axis('off');ax.set(xlim=(0,1),ylim=(0,1))
    cols=[.012,.262,.512,.762];w=.224
    ys=[.89,.735,.58,.425,.27,.115];h=.105
    labels=['Source stays / admissions','Eligible adult first ICU stay','Vasopressor treatment','Cuff and arterial reference MAP\nrecorded during treatment','At least one valid ±5-min pair','Dose ascertainment among pairs']
    stages=['all_source_stays','adult_first_icu','vasopressor_patients','both_bp_during_therapy','paired5_patients']
    flow_display=[]
    for x,db in zip(cols,DBS):
        ax.text(x+w/2,.975,NAMES[db],ha='center',va='top',fontweight='bold',fontsize=8.2,color=COL.get(db,'#333333'))
        prev=None
        for k,(y,label) in enumerate(zip(ys,labels)):
            face='#ECF1F4' if k==4 else '#F7F7F7'
            ax.add_patch(FancyBboxPatch((x,y-h/2),w,h,boxstyle='round,pad=0.005,rounding_size=0.006',edgecolor='#A5ADB2',facecolor=face,linewidth=.65))
            if k<5:
                r=one(FLOW,database=db,stage=stages[k]);v=int(r.value)
                trace('Figure1',NAMES[db]+': '+stages[k],r,['database','stage'],['value','unit'],'integer count')
                unit='stays' if k==0 and db in ['mimic','eicu'] else ('admissions' if k==0 else 'patients')
                value=f'{v:,} {unit}'
                if k==4:
                    q=one(FLOW,database=db,stage='paired5_readings')
                    trace('Figure1',NAMES[db]+': paired5_readings',q,['database','stage'],['value','unit'],'integer count')
                    value+=f'\n{int(q.value):,} pairs'
                if prev is not None:check('flow_nonincreasing:'+db+':'+stages[k],v<=prev)
                prev=v
                flow_display.append(dict(database=db,display_stage=label,value=v,unit=unit,source_key=f'database={db};stage={stages[k]}'))
            else:
                known=0
                for st in STRATA:
                    q=one(FLOW,database=db,stage='phase1_dose_'+st);known+=int(q.value)
                    trace('Figure1',NAMES[db]+': known-dose sum',q,['database','stage'],['value','unit'],'sum five mutually exclusive dose strata')
                u=one(FLOW,database=db,stage='phase1_dose_unknown');total=one(FLOW,database=db,stage='paired5_readings')
                check('flow_dose_balance:'+db,known+int(u.value)==int(total.value))
                trace('Figure1',NAMES[db]+': dose unknown',u,['database','stage'],['value','unit'],'integer count')
                value=f'Known: {known:,} pairs\nUnknown: {int(u.value):,} pairs'
                flow_display.append(dict(database=db,display_stage='Known-dose pairs',value=known,unit='pairs',source_key='sum of five dose strata'))
                flow_display.append(dict(database=db,display_stage='Unknown-dose pairs',value=int(u.value),unit='pairs',source_key=f'database={db};stage=phase1_dose_unknown'))
            ax.text(x+w/2,y+.02,label,ha='center',va='center',fontsize=6.25,linespacing=1.15,color='#333333')
            ax.text(x+w/2,y-.022,value,ha='center',va='center',fontsize=7.0,fontweight='bold',linespacing=1.25)
            if k<5:ax.annotate('',xy=(x+w/2,ys[k+1]+h/2+.008),xytext=(x+w/2,y-h/2-.008),arrowprops={'arrowstyle':'-|>','lw':.7,'color':'#7A858B','mutation_scale':7})
    ax.text(.012,.015,'Patients contribute one index stay and may contribute repeated pairs. Care-setting definitions differ by database.',fontsize=6.4,va='bottom',color='#444444')
    pd.DataFrame(flow_display).to_csv(TAB/'Figure1_flow_display.csv',index=False,encoding='utf-8-sig')
    save(fig,'Figure1_cohort_flow')

def fmt_num(x,digits=1):
    if pd.isna(x):return 'Not available'
    s=f'{x:.{digits}f}'
    return s if digits==3 else s.rstrip('0').rstrip('.')

def markdown_table(df):
    escape=lambda x:str(x).replace('|','\\|').replace('\n','<br>')
    return '| '+' | '.join(map(escape,df.columns))+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+ '\n'.join('| '+' | '.join(map(escape,row))+' |' for row in df.itertuples(index=False,name=None))

def missing_display(fraction):
    return '<0.1%' if 0<float(fraction)<.001 else f'{fraction*100:.1f}%'

def table1():
    specs=[('Patients','patients',None),('Age, years','age','median_IQR'),('Female sex','sex','F'),
           ('Weight, kg','weight','median_IQR'),('BMI, kg/m²','bmi','median_IQR'),('BMI >40 kg/m²','bmi_gt40','yes'),
           ('Recorded atrial fibrillation','atrial_fibrillation','yes'),
           ('Severity score: SOFA / APACHE IVa','severity_score','median_IQR'),
           ('Norepinephrine','vasopressor_use','norepinephrine'),('Epinephrine','vasopressor_use','epinephrine'),
           ('Phenylephrine','vasopressor_use','phenylephrine'),('Dopamine','vasopressor_use','dopamine'),
           ('Vasopressin','vasopressor_use','vasopressin'),('Angiotensin II','vasopressor_use','angiotensin_ii'),
           ('Maximum observed NE-eq, µg/kg/min','max_nee_observed','median_IQR'),
           ('Treatment time with known NE-eq, %','known_duration_fraction','median_IQR'),
           ('Pairs per patient','pairs_per_patient','median_IQR'),('In-hospital death','hospital_death','yes')]
    display=[];long=[]
    for label,var,cat in specs:
        row={'Characteristic':label}
        for db in DBS:
            if var=='patients':
                r=one(BASELINE,database=db,variable='age',category='median_IQR');value=f'{int(r.patients):,}'
                fields=['patients']
            else:
                r=one(BASELINE,database=db,variable=var,category=cat)
                if int(r.n_observed)==0:
                    value='Not available; 100% missing'
                elif cat=='median_IQR':
                    mult=100 if var=='known_duration_fraction' else 1
                    dg=3 if var=='max_nee_observed' else 1
                    value=f'{fmt_num(r["median"]*mult,dg)} [{fmt_num(r.q1*mult,dg)}–{fmt_num(r.q3*mult,dg)}]; {missing_display(r.missing_fraction)} missing'
                else:
                    value=f'{int(r.numerator):,}/{int(r.denominator):,} ({r.percent:.1f}%)'
                    if r.n_missing>0:value+=f'; {missing_display(r.missing_fraction)} missing'
                fields=['patients','n_observed','n_missing','missing_fraction','numerator','denominator','percent','median','q1','q3']
            row[NAMES[db]]=value
            trace('Table1',label+': '+db,r,['database','variable','category'],fields,'median [IQR]; binary n/available denominator (%); explicit missing percentage; duration fraction ×100')
            d={k:clean(v) for k,v in r.items() if not str(k).startswith('_')};d.update(display_label=label,display_value=value);long.append(d)
        display.append(row)
    df=pd.DataFrame(display);df.to_csv(TAB/'Table1_baseline_display.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(long).to_csv(TAB/'Table1_baseline_source_values.csv',index=False,encoding='utf-8-sig')
    foot='''Values are median [interquartile range] or n/number with an available value (%). Missing percentages are shown for continuous variables and incomplete categorical variables. Drug-use categories overlap and represent documented treatment during the index stay. NE-eq, norepinephrine-equivalent dose; BMI, body mass index; SOFA, Sequential Organ Failure Assessment; APACHE, Acute Physiology and Chronic Health Evaluation.\n\nAmsterdamUMCdb age, weight, and BMI use closed-interval midpoints; open intervals are unavailable. SICdb retains source age and anthropometric coarsening. BMI >40 is not classified in AmsterdamUMCdb. Recorded atrial fibrillation uses rhythm or diagnosis records in MIMIC-IV, diagnosis coverage in eICU, and rhythm records in AmsterdamUMCdb; it is unavailable in SICdb.\n\nSOFA is the available first-day score in MIMIC-IV, and APACHE IVa is the available score in eICU. Scores are not interchangeable. In-hospital mortality is descriptive; AmsterdamUMCdb did not provide a suitable hospital-discharge outcome for this table.\n\nMaximum observed NE-eq is the maximum summed dose within ascertainable recorded or inferred treatment intervals, not a maximum restricted to blood-pressure pairs. Intervals with an unknown component remain unknown. The accompanying known-dose treatment-time fraction is calculated per patient; unobserved intervals can contain higher doses. eICU treatment intervals are inferred from positive-dose records, with continuation limited to the next same-drug record or 60 minutes.\n\nSource-specific admission and care-unit categories, original age groups, treatment duration, and dose-observation completeness are provided in the supplementary baseline detail table. No common admission-type scale is imposed.\n'''
    (TAB/'Table1_baseline.md').write_text('# Table 1. Characteristics of patients with paired cuff and arterial MAP measurements\n\n'+markdown_table(df)+'\n\n'+foot,encoding='utf-8')
    (TAB/'Table1_notes.txt').write_text(foot,encoding='utf-8')
    extra=BASELINE[BASELINE.variable.isin(['admission_source','admission_type','unit_type','age_group','treatment_hours','dose_observation_status'])].copy()
    extra.drop(columns=['_source','_row']).to_csv(TAB/'TableS_baseline_detail_source_values.csv',index=False,encoding='utf-8-sig')
    misses=BASELINE[['database','variable','patients','n_observed','n_missing','missing_fraction']].drop_duplicates()
    check('Table1_missingness_unique',not misses.duplicated(['database','variable']).any())
    for r in misses.itertuples():
        check('Table1_observed_plus_missing:'+r.database+':'+r.variable,r.n_observed+r.n_missing==r.patients)
        check('Table1_missing_fraction:'+r.database+':'+r.variable,np.isclose(r.missing_fraction,r.n_missing/r.patients))
    misses.to_csv(TAB/'TableS_baseline_missingness.csv',index=False,encoding='utf-8-sig')
    extra_names={'admission_source':'Admission source','admission_type':'Admission type',
                 'unit_type':'Care unit','age_group':'Source age group, years','treatment_hours':'Treatment duration, hours',
                 'dose_observation_status':'Dose ascertainment across treatment episodes'}
    cats={'complete_within_defined_episodes':'Complete within recorded or inferred treatment episodes',
          'no_known_segment':'No ascertainable-dose interval','partial_observation':'Partially ascertainable',
          'planned':'Planned','unplanned':'Unplanned','median_IQR':'Median [interquartile range]',
          'Andere locatie zelfde ziekenhuis, transport per ambulance':'Other site of same hospital (ambulance transfer)',
          'Anders':'Other','CCU/IC ander ziekenhuis':'Coronary/intensive care unit, another hospital',
          'CCU/IC zelfde ziekenhuis':'Coronary/intensive care unit, same hospital',
          'Eerste Hulp afdeling ander ziekenhuis':'Emergency department, another hospital',
          'Eerste Hulp afdeling zelfde ziekenhuis':'Emergency department, same hospital','Huis':'Home',
          'Operatiekamer vanaf verpleegafdeling zelfde ziekenhuis':'Operating theatre, from ward of same hospital',
          'Recovery zelfde ziekenhuis (alleen bij niet geplande IC-opname)':'Recovery unit, same hospital (unplanned ICU admission)',
          'Special/Medium care ander ziekenhuis':'Special/intermediate care, another hospital',
          'Special/Medium care zelfde ziekenhuis':'Special/intermediate care, same hospital',
          'Verpleegafdeling ander ziekenhuis':'Ward, another hospital','Verpleegafdeling zelfde ziekenhuis':'Ward, same hospital',
          '1. Medizin':'Medical department 1','2. Medizin':'Medical department 2','3. Medizin':'Medical department 3',
          'Allgemeinchirurgie':'General surgery','Dermatologie':'Dermatology','Externe Intensivstation':'External intensive care unit',
          'Externes Krankenhaus':'External hospital','Gefäßchirurgie':'Vascular surgery','Gynäkologie':'Gynaecology',
          'HNO':'Otorhinolaryngology','Kieferchirurgie':'Maxillofacial surgery','Neurologie':'Neurology',
          'Pneumologie':'Pulmonology','Unfallchirurgie':'Trauma surgery','Urologie':'Urology','Zentrale Notaufnahme':'Central emergency department'}
    extras=[]
    for _,r in extra.iterrows():
        char=extra_names[r.variable]
        if r.variable=='unit_type' and r.database=='sicdb':char='Last recorded unit (source label)'
        if r.variable=='admission_type':char={'mimic':'Hospital admission type','eicu':'ICU stay type','amsterdam':'Admission urgency','sicdb':'Surgical admission type'}[r.database]
        value=(f'{fmt_num(r["median"])} [{fmt_num(r.q1)}–{fmt_num(r.q3)}]' if r.category=='median_IQR'
               else f'{int(r.numerator):,}/{int(r.denominator):,} ({r.percent:.1f}%)')
        extras.append({'Database':NAMES[r.database],'Characteristic':char,'Category':cats.get(r.category,r.category),
                       'Value':value,'Missing':missing_display(r.missing_fraction)})
        trace('TableS_baseline_detail',r.database+': '+r.variable+'; '+r.category,r,['database','variable','category'],
              ['numerator','denominator','percent','median','q1','q3','n_missing','missing_fraction'],'source categories translated for display; no aggregation or harmonization')
    ext_df=pd.DataFrame(extras)
    ext_df.to_csv(TAB/'TableS_baseline_detail_display.csv',index=False,encoding='utf-8-sig')
    (TAB/'TableS_baseline_detail.md').write_text('# Supplementary baseline detail\n\n'+markdown_table(ext_df)+
       '\n\nValues retain the source-specific admission categories. SICdb unit labels describe the last recorded case unit and are not first-unit admission labels. Dose ascertainment refers to the recorded or inferred treatment episodes; missing intervals were not imputed. Source labels in Dutch and German are translated for presentation, with the original values retained in the accompanying data.\n',encoding='utf-8')

def rate_display(r):
    n,d=int(r.numerator),int(r.denominator)
    if d==0:return f'{n}/{d}; not estimable'
    suffix=''
    if r.get('precision')=='insufficient' or d<30:suffix+='†'
    if bool(r.get('boundary_events',False)) or bool(r.get('boundary',False)):suffix+='‡'
    return f'{n:,}/{d:,}; {r.estimate*100:.1f} ({r.ci_lower*100:.1f}–{r.ci_upper*100:.1f}){suffix}'

def table2():
    display=[];long=[]
    for db in ['mimic','eicu','pooled']:
        for st in ['all',*STRATA]:
            actual='known' if db=='pooled' and st=='all' else st
            a=primary_row(db,actual);d=deep_row(db,actual)
            for metric,r in [('missed hypotension rate',a),('deep missed hypotension rate',d)]:
                ratio_check(r,f'Table2:{db}:{actual}:{metric}',db=='pooled')
                keys=['stratum','endpoint'] if metric.startswith('missed') else ['stratum','metric']
                if db!='pooled':keys=['database',*keys]
                fields=['numerator','denominator','estimate','ci_lower','ci_upper']
                if db=='pooled':fields+=['I2','databases']
                trace('Table2',db+': '+actual+'; '+metric,r,keys,fields,'n/N; percentage (95% CI); pooled rates are meta-analytic, not summed n/N')
                long.append(dict(database=db,stratum=actual,metric=metric,numerator=int(r.numerator),denominator=int(r.denominator),
                  estimate=clean(r.estimate),ci_lower=clean(r.ci_lower),ci_upper=clean(r.ci_upper),
                  I2=clean(r.get('I2',np.nan)),databases=clean(r.get('databases',np.nan)),
                  precision=r.get('precision',''),boundary=clean(r.get('boundary_events',r.get('boundary',False))),
                  source_csv=r['_source'],source_csv_row=int(r['_row'])))
            st_label='Known-dose overall' if db=='pooled' and st=='all' else ('Overall' if st=='all' else TICKS[STRATA.index(st)])
            display.append({'Database':NAMES[db],'NE-eq, µg/kg/min':st_label,
              'Missed hypotension: n/N; % (95% CI)':rate_display(a),'Deep missed hypotension: n/N; % (95% CI)':rate_display(d),
              'I², missed / deep (%)':f'{a.I2:.1f} / {d.I2:.1f}' if db=='pooled' else '—'})
    df=pd.DataFrame(display);df.to_csv(TAB/'Table2_missed_and_deep_display.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(long).to_csv(TAB/'Table2_missed_and_deep_source_values.csv',index=False,encoding='utf-8-sig')
    foot='''The missed hypotension rate is the proportion of arterial reference MAP <65 mmHg pairs with cuff MAP ≥65 mmHg. The deep missed hypotension rate uses arterial reference MAP <55 mmHg as its denominator. Each n/N is a count of paired readings, not patients. Database-specific 95% confidence intervals use 2,000 patient-cluster bootstrap replicates.\n\nDatabase-specific overall rows include pairs with unknown NE-eq. Five dose strata and the pooled overall row include known-dose pairs. The pooled estimate uses a random-effects model on the logit scale, REML heterogeneity estimation, and modified Hartung–Knapp intervals; pooled n/N shows the summed constituent counts and is not the calculation of the pooled estimate. Four databases contribute to each pooled row except NE-eq >1.0 µg/kg/min, which includes MIMIC-IV, eICU, and AmsterdamUMCdb.\n\n† Insufficient precision: fewer than 30 denominator pairs or a 95% confidence interval wider than 30 percentage points. ‡ Boundary estimate: zero or all events in the relevant denominator; a degenerate bootstrap interval does not imply certainty. A row with no denominator is not estimable. NE-eq, norepinephrine-equivalent dose.\n'''
    foot+='\nDatabase-specific results for AmsterdamUMCdb and SICdb are provided in the Supplementary Material and its electronic data.\n'
    (TAB/'Table2_missed_and_deep.md').write_text('# Table 2. Missed and deep missed hypotension by norepinephrine-equivalent dose\n\n'+markdown_table(df)+'\n\n'+foot,encoding='utf-8')
    (TAB/'Table2_notes.txt').write_text(foot,encoding='utf-8')

def captions():
    dose_note_path=FIG/'everdac_caption_note.txt'
    dose_note=(dose_note_path.read_text(encoding='utf-8').strip() if dose_note_path.exists() else
      'EVERDAC used a norepinephrine-tartrate-plus-epinephrine dose threshold of 2.5 µg/kg/min. For norepinephrine alone, 2.5 µg/kg/min of tartrate is approximately 1.25 µg/kg/min of base; this is not a universal conversion of a combined vasopressor NE-eq dose. The shaded categories identify the dose range under-represented in EVERDAC, not a trial exclusion boundary.')
    texts={
    'Figure1_cohort_flow':('Figure 1. Cohort construction across four critical-care databases. Each database begins with its available source stays or admissions and applies its existing adult, first-stay, treatment, blood-pressure availability, and pairing definitions. Patients with at least one valid pair form the analytical cohort; patients can contribute repeated pairs. The terminal dose counts refer to pairs. MIMIC-IV uses source ICU-stay membership, eICU uses unit stays, AmsterdamUMCdb uses pure intensive-care admissions, and SICdb uses strictly ascertainable ICU cases. Counts from different database definitions are presented separately; branches are not combined as interchangeable exclusions.'),
    'Figure2_dose_and_cuff_probability':('Figure 2. Cuff and arterial hypotension discordance across norepinephrine-equivalent dose strata. Panels A–C show the proportion of cuff MAP ≥65 mmHg when arterial reference MAP is <65 mmHg (solid circles; missed hypotension rate) or <55 mmHg (dashed triangles; deep missed hypotension rate). Panels D–E show the probability of arterial hypotension given a cuff MAP ≥65 mmHg. The percentages beneath panels D–E are the arterial reference MAP <65 mmHg prevalence among all paired readings in each dose stratum, not among normal cuff readings alone. All dose strata use original known-dose pairs. Points are pair-level estimates and bars are 95% confidence intervals from 2,000 patient-cluster bootstrap replicates in each database. Panel C uses random-effects pooling on the logit scale with REML and modified Hartung–Knapp intervals; the >1.0 µg/kg/min stratum contains three databases (MIMIC-IV, eICU, and AmsterdamUMCdb), and other strata contain four. The conditional probabilities in panels D–E are not pooled. '+dose_note+' NE-eq, norepinephrine-equivalent dose; MAP, mean arterial pressure.'),
    'Figure3_cuff_reading_bins':('Figure 3. Probability of arterial hypotension across cuff MAP readings. Each point is the proportion of pairs with arterial reference MAP <65 mmHg within the indicated cuff MAP interval and known-dose group. Panels show MIMIC-IV and eICU separately, with no cross-database pooling. Bars are 95% confidence intervals from 2,000 patient-cluster bootstrap replicates. The lower annotation reports the arterial reference MAP <65 mmHg prevalence among all pairs in each dose group. Labels 65–69, 70–74, and 75–79 denote [65,70), [70,75), and [75,80) mmHg, respectively; values were not rounded before classification. NE-eq is expressed in µg/kg/min.'),
    'ICM_Figure1_two_panel':('Figure 1. Dose-stratified cuff and arterial MAP discordance. Panel A shows pooled proportions of cuff MAP ≥65 mmHg when arterial reference MAP is <65 mmHg (solid circles) or <55 mmHg (dashed triangles). Panel B shows the probability of arterial hypotension given a cuff MAP ≥65 mmHg in MIMIC-IV and eICU; the rows below report arterial reference MAP <65 mmHg prevalence in each dose stratum. All strata use original known-dose pairs. Bars are 95% confidence intervals: patient-cluster bootstrap in panel B, and logit random-effects pooling with REML and modified Hartung–Knapp intervals in panel A. Four databases contribute to panel A except the >1.0 µg/kg/min stratum, which includes three. Panel B is not pooled. '+dose_note+' NE-eq, norepinephrine-equivalent dose; MAP, mean arterial pressure.')}
    for name,txt in texts.items():(FIG/f'{name}_caption.txt').write_text(txt+'\n',encoding='utf-8')
    (FIG/'captions.md').write_text('\n\n'.join(texts.values())+'\n',encoding='utf-8')
    params=read('phase2/phase2_parameters.csv')
    for artifact in ['Figure2','Figure3','ICM_Figure1','Table2']:
        r=one(params,parameter='bootstrap_replicates')
        trace(artifact,'caption or footnote: bootstrap replicates',r,['parameter'],['value'],'count unchanged')
    dose=read('phase3/sources/everdac_dose_review.csv')
    for artifact in ['Figure2','ICM_Figure1']:
        for item in ['eligibility_and_rescue_threshold','norepinephrine_tartrate_to_base_ratio','pure_norepinephrine_base_threshold']:
            r=one(dose,item=item)
            trace(artifact,'caption: EVERDAC dose convention',r,['item'],['value','unit','verification_status','source_location'],'verified external number or pure-norepinephrine conversion')

def finish():
    trace_df=pd.DataFrame(SOURCE_ROWS)
    trace_df[trace_df.artifact.str.startswith(('Figure','ICM'))].to_csv(OUT/'figure_sources.csv',index=False,encoding='utf-8-sig')
    trace_df[trace_df.artifact.str.startswith('Table')].to_csv(OUT/'table_sources.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(CHECKS).to_csv(OUT/'figure_table_validation.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(EXPORTS).to_csv(OUT/'figure_exports.csv',index=False,encoding='utf-8-sig')
    hashes=[]
    for p in sorted(set(trace_df.source_csv)):
        hashes.append(dict(source_csv=p,sha256=hashlib.sha256((ROOT/p).read_bytes()).hexdigest()))
    pd.DataFrame(hashes).to_csv(OUT/'figure_table_input_manifest.csv',index=False,encoding='utf-8-sig')
    print(json.dumps({'figures':len(EXPORTS),'trace_rows':len(SOURCE_ROWS),'checks':len(CHECKS),'all_passed':all(x['passed'] for x in CHECKS)},ensure_ascii=False))

if __name__=='__main__':
    figure1();figure2();figure3();figure_icm();table1();table2();captions();finish()
