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

"""Resolve concepts from full dictionaries; never embed item IDs."""
from pathlib import Path
import csv,gzip,re,json
ROOT=Path('@CUFF_WORKSPACE@'); OUT=ROOT/'phase0'
def read(p):
    op=gzip.open if str(p).endswith('.gz') else open
    with op(p,'rt',encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(p,rows,fields=None):
    fields=fields or list(dict.fromkeys(k for x in rows for k in x))
    with open(p,'w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def drug(s):
    s=s.lower()
    if re.search(r'norepine|noradren|levophed',s):return 'norepinephrine'
    if re.search(r'epineph|epinefr|adrenalin',s):return 'epinephrine'
    if re.search(r'phenyleph|fenylefr|neo.?syneph',s):return 'phenylephrine'
    if re.search(r'dopamin|inotropin',s):return 'dopamine'
    if re.search(r'vasopress|pitressin',s):return 'vasopressin'
    if re.search(r'angiotensin|giapreza',s):return 'angiotensin_ii'
    return ''
def bp_class(s):
    s=s.lower()
    if re.search(r'alarm|streef|target|iapb|iabp|pulmonary|\bpap\b|apache|\ba_abp|\bra_abp|\bmca_abp',s):return 'other','',''
    if re.search(r'non.?invasive|non invasive|niet invasieve|bloodpressureni|\bnibp\b|\bnbp\b',s):mod='cuff'
    elif re.search(r'arterial blood pressure|art bp |\babp (systol|diastol|gemiddeld)|bloodpressurearterial|(?:invasive|mean|systolic|diastolic).*arterial pressure|picco ap',s):mod='arterial'
    else:return 'other','',''
    if re.search(r'systol|picco aps',s):component='sbp'
    elif re.search(r'diastol|picco apd',s):component='dbp'
    elif re.search(r'mean|gemiddeld|map|picco apm',s):component='map'
    else:return 'other','',''
    fam=('picco' if 'picco' in s else ('art2' if re.search(r'\bii\b',s) else 'art1')) if mod=='arterial' else 'cuff'
    return mod,component,fam
def main():
    allm=read(OUT/'mimic_dictionary_all.csv')
    b=[];v=[];cov=[]
    for x in allm:
        s=x['label'];mod,comp,fam=bp_class(s)
        if re.search(r'arterial|\bart\b|\babp|\bnibp|\bnbp|non.?invasive|non invasive|blood pressure',s,re.I):
            b.append({'itemid':x['itemid'],'label':s,'abbreviation':x['abbreviation'],'category':x['category'],'unitname':x['unitname'],'linksto':x['linksto'],'classification':mod,'component':comp,'family':fam,'selected':int(bool(comp) and x['linksto']=='chartevents'),'rows':'[待填：提取后计数]','human_review':'pending'})
        d=drug(s)
        if d:v.append({'itemid':x['itemid'],'label':s,'drug':d,'unitname':x['unitname'],'linksto':x['linksto'],'selected':int(x['linksto']=='inputevents'),'rows':'[待填：提取后计数]','human_review':'pending'})
        if s in ['Heart Rhythm','Admission Weight (Kg)','Admission Weight (lbs.)','Daily Weight','Height','Height (cm)']:
            cov.append({'itemid':x['itemid'],'label':s,'unitname':x['unitname']})
    write(OUT/'map_bp_mimic.csv',b);write(OUT/'map_vaso_mimic.csv',v);write(OUT/'map_covariates_mimic.csv',cov)
    a=read(OUT/'amsterdam_dictionary_all.csv');u={}
    for x in a:
        if x['itemid'] not in u:u[x['itemid']]=dict(x)
        elif x['ordercategory'] and x['ordercategory'] not in u[x['itemid']]['ordercategory']:u[x['itemid']]['ordercategory']+=' | '+x['ordercategory']
    b=[];v=[]
    for x in u.values():
        s=x['item'];mod,comp,fam=bp_class(s+' '+x['item_en'])
        if re.search(r'abp|art|bloeddruk|nibp|niet invasi|gemiddeld|blood pressure|picco ap',s+' '+x['item_en'],re.I):
            b.append({'itemid':x['itemid'],'label':s,'label_en':x['item_en'],'category':x['category'],'unitname':x['unit'],'classification':mod,'component':comp,'family':fam,'selected':int(bool(comp)),'rows':'[待填：提取后计数]','human_review':'pending'})
        d=drug(s+' '+x['item_en'])
        if d:v.append({'itemid':x['itemid'],'label':s,'drug':d,'ordercategory':x['ordercategory'],'selected':int(x['category']=='Medicatie' and 'Spuitpompen' in x['ordercategory']),'rows':'[待填：提取后计数]','human_review':'pending'})
    write(OUT/'map_bp_amsterdam.csv',b);write(OUT/'map_vaso_amsterdam.csv',v)
    sic=read(Path('@SICDB_ROOT@/d_references.csv.gz'))
    b=[];v=[]
    for x in sic:
        s=x['ReferenceValue'];mod,comp,fam=bp_class(s)
        if re.search(r'bloodpressure|arterial|noninvasive|\bnibp\b|\babp\b',s+' '+x['ReferenceDescription'],re.I):
            b.append({'itemid':x['ReferenceGlobalID'],'label':s,'category':x['ReferenceName'],'unitname':x['ReferenceUnit'],'classification':mod,'component':comp,'family':fam,'selected':int(bool(comp) and x['ReferenceName']=='SignalFloat'),'rows':'[待填：提取后计数]','human_review':'pending'})
        d=drug(s)
        if d:v.append({'itemid':x['ReferenceGlobalID'],'label':s,'drug':d,'unitname':x['ReferenceUnit'],'category':x['ReferenceName'],'selected':int(x['ReferenceName']=='Drug' and not re.search(r'lidoc|bupiv|articain',s,re.I)),'rows':'[待填：提取后计数]','human_review':'pending'})
    write(OUT/'map_bp_sicdb.csv',b);write(OUT/'map_vaso_sicdb.csv',v)
    e=read(OUT/'eicu_drug_dictionary_all.csv');v=[]
    for x in e:
        d=drug(x['drugname'])
        if d:v.append({'itemid':x['drugname'],'label':x['drugname'],'drug':d,'unitname':' | '.join(re.findall(r'\(([^()]*)\)',x['drugname'])),'selected':int(not re.search(r'lidoc|bupiv|volume\s*\(ml\)',x['drugname'],re.I)),'rows':x['rows'],'human_review':'pending'})
    write(OUT/'map_vaso_eicu.csv',v)
    b=[]
    for table,mod,prefix in [('vitalperiodic','arterial','systemic'),('vitalaperiodic','cuff','noninvasive')]:
        for suffix,comp in [('mean','map'),('systolic','sbp'),('diastolic','dbp')]:
            b.append({'itemid':table+'.'+prefix+suffix,'label':prefix+suffix,'classification':mod,'component':comp,'unitname':'mmHg','selected':1,'rows':'[待填：提取后计数]','human_review':'pending'})
    write(OUT/'map_bp_eicu.csv',b)
    print(json.dumps({k:sum(x.get('selected')=='1' for x in read(OUT/k)) for k in ['map_bp_mimic.csv','map_bp_amsterdam.csv','map_bp_sicdb.csv','map_vaso_eicu.csv']},ensure_ascii=False))
if __name__=='__main__':main()
