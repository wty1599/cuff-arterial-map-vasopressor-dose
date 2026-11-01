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
import sys,csv,json,re,math,datetime,traceback

import duckdb,pandas as pd,numpy as np
ROOT=Path('@CUFF_WORKSPACE@');OUT=ROOT/'phase0';R=OUT/'_restricted'
SIC=Path('@SICDB_ROOT@')
P=json.loads((ROOT/'sql/phase0_parameters.json').read_text(encoding='utf-8'))
def connect(db):
    c=duckdb.connect(str(R/f'{db}.duckdb'));c.execute("SET memory_limit='2GB';SET threads=2;SET preserve_insertion_order=false;")
    c.execute(f"SET temp_directory='{R.as_posix()}/{db}_tmp';")
    return c
def csvq(p,**args):
    args={'all_varchar':True,**args};opts=','.join(k+'='+('true' if v is True else 'false' if v is False else "'"+str(v).replace("'","''")+"'") for k,v in args.items())
    return f"read_csv_auto('{Path(p).as_posix()}',{opts})"
def save(c,table,df):
    df=df.copy()
    for col in df:
        if isinstance(df[col].dtype,pd.StringDtype):df[col]=df[col].astype(object).where(df[col].notna(),None)
    c.register('_frame',df);c.execute(f'CREATE OR REPLACE TABLE {table} AS SELECT * FROM _frame');c.unregister('_frame')
def read(p,**kw):return pd.read_csv(p,keep_default_na=False,**kw)
def number(x):
    try:v=float(x);return v if np.isfinite(v) else np.nan
    except (ValueError,TypeError):return np.nan
def num(s):return pd.to_numeric(s,errors='coerce')
def epoch(s):
    d=pd.to_datetime(s,errors='coerce');return d.astype('datetime64[s]').astype('int64').astype(float).where(d.notna(),np.nan)
def midpoint(x):
    nums=re.findall(r'\d+(?:\.\d+)?',str(x))
    return (float(nums[0])+float(nums[1]))/2 if len(nums)==2 else np.nan
def bounds(x):
    nums=re.findall(r'\d+(?:\.\d+)?',str(x))
    if len(nums)==2:return float(nums[0]),float(nums[1])
    if len(nums)==1:return (float(nums[0]),np.inf) if '+' in str(x) else (0,float(nums[0]))
    return np.nan,np.nan
def normalize_unit(s):
    return str(s).lower().replace('μ','mc').replace('µ','mc').replace('micrograms','mcg').replace('nanograms','ng').replace('hour','hr').replace('uur','hr').replace('minutes','min').replace('minute','min').replace('units','u').replace('unit','u').replace('ie','u').replace(' ','')
def convert_rate(rate,unit,weight,drug):
    r=number(rate);w=number(weight);u=normalize_unit(unit)
    if not np.isfinite(r) or r<=0:return np.nan,'nonpositive_or_missing_rate'
    if re.search(r'ml|unknown',u) or not u:return np.nan,'unit_or_concentration_unknown'
    perkg='/kg/' in u;perhour=u.endswith('/hr')
    if not (u.endswith('/min') or perhour):return np.nan,'time_unit_unknown'
    if drug=='vasopressin':
        if not u.startswith('u/'):return np.nan,'vasopressin_unit_unknown'
        if perkg:
            if not np.isfinite(w) or w<=0:return np.nan,'weight_missing'
            r*=w
        return r/(60 if perhour else 1)*P['vasopressin_units_min_factor'],'known'
    mass={'mcg':1,'ng':.001,'mg':1000,'g':1000000}.get(u.split('/')[0])
    if mass is None:return np.nan,'mass_unit_unknown'
    r=r*mass/(60 if perhour else 1)
    if not perkg:
        if not np.isfinite(w) or w<=0:return np.nan,'weight_missing'
        r/=w
    factor=P['angiotensin_ii_ng_kg_min_factor']*1000 if drug=='angiotensin_ii' else P['nee_factors_mcg_kg_min'][drug]
    return r*factor,'known'
def run_guard(db,fn):
    try:fn(db)
    except Exception as e:
        (R/f'error_{db}.txt').write_text(traceback.format_exc(),encoding='utf-8')
        print(json.dumps({'db':db,'status':'error','error_class':type(e).__name__,'detail_file':str(R/f'error_{db}.txt')}));raise SystemExit(1)
def mapping(db,kind):return read(OUT/f'map_{kind}_{db}.csv',dtype=str)
def exported(c,name,q):
    d=c.execute(q).df();d.to_csv(OUT/name,index=False);return d
