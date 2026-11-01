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

import os as _cuff_env
from pathlib import Path
import sys,json,hashlib
ROOT=Path('@CUFF_WORKSPACE@');OUT=ROOT/'phase5';RES=OUT/'_restricted'
import pandas as pd,numpy as np
summary=[];flow=[]
for db in ['mimic','eicu']:
 d=pd.read_parquet(ROOT/f'phase1/_restricted/pairs_{db}.parquet')
 z=d.loc[d.nee.notna(),['patient','cid','art_map','cuff_map','nee']].copy()
 z.to_csv(RES/f'audit_original_pairs_{db}.csv',index=False)
 pats,patient_index=np.unique(z.patient.to_numpy(),return_inverse=True)
 patient_copies=np.random.default_rng(20261027).multinomial(len(pats),np.full(len(pats),1/len(pats)))
 pd.DataFrame({'copies':patient_copies[patient_index]}).to_csv(RES/f'audit_first_cluster_weights_{db}.csv',index=False)
 assert patient_copies.sum()==len(pats)
 all_low=d.art_map.lt(65);known_low=z.art_map.lt(65)
 summary.append(dict(database=db,all_pairs=len(d),all_patients=d.patient.nunique(),known_pairs=len(z),known_patients=z.patient.nunique(),known_pair_fraction=len(z)/len(d),known_patient_fraction=z.patient.nunique()/d.patient.nunique(),all_low_arterial_pairs=int(all_low.sum()),all_low_arterial_fraction=all_low.mean(),known_low_arterial_pairs=int(known_low.sum()),known_low_arterial_fraction=known_low.mean(),known_cuff_ge65=int(z.cuff_map.ge(65).sum()),known_joint_low_normal=int((known_low&z.cuff_map.ge(65)).sum()),pair_fields_missing=int(z.isna().sum().sum()),arterial_min=z.art_map.min(),arterial_max=z.art_map.max(),cuff_min=z.cuff_map.min(),cuff_max=z.cuff_map.max(),nee_min=z.nee.min(),nee_max=z.nee.max(),unique_pair_ids=d.cid.is_unique,original_pair_sha256=hashlib.sha256((ROOT/f'phase1/_restricted/pairs_{db}.parquet').read_bytes()).hexdigest()))
 # Explicit intervals provide a path independent of production searchsorted.
 masks=[z.nee.le(.1),z.nee.gt(.1)&z.nee.le(.25),z.nee.gt(.25)&z.nee.le(.5),z.nee.gt(.5)&z.nee.le(1),z.nee.gt(1)]
 labels=['<=0.1','>0.1-0.25','>0.25-0.5','>0.5-1.0','>1.0']
 assert np.column_stack(masks).sum(axis=1).tolist()==[1]*len(z)
 for label,mask in zip(labels,masks):
  a=z[mask];flow.append(dict(database=db,stratum=label,pairs=len(a),patients=a.patient.nunique(),low_arterial_pairs=int(a.art_map.lt(65).sum()),normal_cuff_pairs=int(a.cuff_map.ge(65).sum()),joint_low_normal=int((a.art_map.lt(65)&a.cuff_map.ge(65)).sum())))
pd.DataFrame(summary).to_csv(OUT/'phase5_independent_population.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(flow).to_csv(OUT/'phase5_independent_dose_counts.csv',index=False,encoding='utf-8-sig')
print(json.dumps({'known_pairs':{x['database']:x['known_pairs'] for x in summary},'row_level_exports':'restricted only'}))
