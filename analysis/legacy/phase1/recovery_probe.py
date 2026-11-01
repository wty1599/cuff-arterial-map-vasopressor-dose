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

from common import *
from pathlib import Path
import json
c=con('eicu')
c.execute(f"CREATE TEMP VIEW raw AS SELECT * FROM read_csv_auto('{(SRC/'eicu_drugs.csv').as_posix()}',all_varchar=true)")
q="""WITH eligible AS (SELECT DISTINCT d.sid,d.drug,d.start,d.rawid,d.dose,d.dose_reason FROM drugs d JOIN pairs p ON d.sid=p.sid AND p.t>=d.start AND p.t<d.end WHERE p.is_primary AND p.nee IS NULL)
SELECT r.drugname,e.drug,e.dose_reason,count(*) segments,
count(*) FILTER(WHERE try_cast(r.drugrate AS DOUBLE)>0) drugrate_positive,
count(*) FILTER(WHERE try_cast(r.infusionrate AS DOUBLE)>0) infusionrate_positive,
count(*) FILTER(WHERE try_cast(r.drugamount AS DOUBLE)>0 AND try_cast(r.volumeoffluid AS DOUBLE)>0) amount_volume_positive,
count(*) FILTER(WHERE regexp_matches(lower(r.drugamount),'mg|mcg|unit')) amount_unit_text,
count(*) FILTER(WHERE try_cast(r.patientweight AS DOUBLE) BETWEEN 30 AND 300) infusion_weight_valid
FROM eligible e JOIN raw r ON e.rawid=r.infusiondrugid GROUP BY 1,2,3 ORDER BY segments DESC"""
d=c.execute(q).df();save(d,'phase1_eicu_recovery_probe.csv')
print(d.head(25).to_string(index=False))
c.execute(f"COPY (SELECT DISTINCT sid FROM pairs WHERE is_primary AND nee IS NULL) TO '{(RES/'recovery_candidate_stays.csv').as_posix()}' (HEADER,DELIMITER ',')")
print(json.dumps({'raw_columns':c.execute('DESCRIBE raw').df().column_name.tolist(),'probe_rows':len(d),'eligible_source_segments':int(d.segments.sum())}))
c.close()
med=RES/'recovery_medication_source.csv'
if med.exists():
    m=pd.read_csv(med,dtype=str,keep_default_na=False)
    print(json.dumps({'medication_columns':list(m.columns),'rows':len(m)}))
    q=m.groupby(['drugname','routeadmin','frequency'],dropna=False).size().reset_index(name='records').sort_values('records',ascending=False)
    save(q,'phase1_eicu_recovery_medication_dictionary.csv')
    print(q.head(35).to_string(index=False))
