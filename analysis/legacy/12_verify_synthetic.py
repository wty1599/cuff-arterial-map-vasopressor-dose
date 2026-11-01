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

from importlib import import_module
g=import_module('08_core');p=import_module('10_pair_and_counts')
tests=[(.1,'mcg/kg/min',70,'norepinephrine',.1),(.04,'units/min',70,'vasopressin',.1),(2.4,'units/hour',70,'vasopressin',.1),(10,'ng/kg/min',70,'angiotensin_ii',.1),(.01,'mcg/kg/min',70,'angiotensin_ii',.1),(1,'mg/hr',50,'norepinephrine',1/3),(5,'mcg/kg/min',70,'phenylephrine',.5),(10,'mcg/kg/min',70,'dopamine',.1),(0.000007,'g/min',70,'norepinephrine',.1)]
for rate,unit,weight,drug,expected in tests:
    actual,status=g.convert_rate(rate,unit,weight,drug);assert abs(actual-expected)<1e-10 and status=='known'
assert g.convert_rate(1,'ml/hr',70,'norepinephrine')[1]!='known'
assert g.convert_rate(1,'mcg/min',float('nan'),'norepinephrine')[1]=='weight_missing'
c=g.duckdb.connect();c.execute('CREATE TABLE drugs(sid VARCHAR,drug VARCHAR,rawid VARCHAR,start DOUBLE,"end" DOUBLE,end30 DOUBLE,dose DOUBLE)')
c.executemany('INSERT INTO drugs VALUES(?,?,?,?,?,?,?)',[('s','A','a',0,60,60,.2),('s','A','b',120,180,180,.2),('s','B','c',70,110,110,.1)])
p.episodes(c,'e60',True,60,'end');p.episodes(c,'e0',True,0,'end')
assert c.execute("SELECT count(*) FROM e60 WHERE drug='A'").fetchone()[0]==1
assert c.execute("SELECT count(*) FROM e0 WHERE drug='A'").fetchone()[0]==2
expected=c.execute('SELECT count(DISTINCT drug) FROM e60 WHERE start<=90 AND stop>90').fetchone()[0]
active=c.execute('SELECT count(DISTINCT drug) FROM drugs WHERE start<=90 AND "end">90').fetchone()[0]
assert expected==2 and active==1 # The total NE-eq at 90 s must remain unknown.
assert c.execute('SELECT median(v) FROM (VALUES(50.0),(70.0)) x(v)').fetchone()[0]==60
g.pd.DataFrame([{'check':'units_and_component_gap_boundaries','passed':True,'fixed_expectations':len(tests)+6}]).to_csv(g.OUT/'synthetic_checks.csv',index=False)
print('Synthetic scientific boundary checks passed.')
