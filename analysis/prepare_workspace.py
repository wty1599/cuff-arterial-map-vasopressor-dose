"""Render legacy code into an external credentialed-data workspace.

Only source code, parameters and dictionary mappings are copied. No analysis runs.
The destination must be outside the repository. Paths come from environment
variables CUFF_WORKSPACE, AMSTERDAM_ROOT and SICDB_ROOT.
"""
from pathlib import Path
import argparse
import json
import os

REPO=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workspace',type=Path);args=parser.parse_args()
    configured=str(args.workspace) if args.workspace else os.environ.get('CUFF_WORKSPACE','')
    if not configured:parser.error('Set CUFF_WORKSPACE or provide --workspace')
    work=Path(configured).expanduser().resolve()
    if work==REPO or REPO in work.parents:parser.error('The private workspace must be outside this repository')
    values={'CUFF_WORKSPACE':work.as_posix(),
            'AMSTERDAM_ROOT':os.environ.get('AMSTERDAM_ROOT',''),
            'SICDB_ROOT':os.environ.get('SICDB_ROOT','')}
    if any(not x for x in values.values()):parser.error('Set AMSTERDAM_ROOT and SICDB_ROOT to credentialed source locations')
    for key in ['AMSTERDAM_ROOT','SICDB_ROOT']:values[key]=Path(values[key]).expanduser().resolve().as_posix()
    if any("'" in value or '\n' in value for value in values.values()):parser.error('Paths may not contain a single quote or newline')
    work.mkdir(parents=True,exist_ok=True)
    def render(src,dest):
        text=src.read_text(encoding='utf-8')
        for key,value in values.items():text=text.replace('@'+key+'@',value)
        if dest.exists():
            if dest.read_text(encoding='utf-8')==text:return
            raise FileExistsError('Existing changed workspace file; preserve it before preparing a new workspace: '+dest.name)
        dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(text,encoding='utf-8')
    for folder in ['sql','analysis/legacy']:
        for src in (REPO/folder).rglob('*'):
            if src.is_file():render(src,work/'sql'/src.relative_to(REPO/folder))
    for record in json.loads((REPO/'config/parameter_targets.json').read_text()):render(REPO/record['file'],work/record['workspace_target'])
    for src in (REPO/'config/mappings').glob('*.csv'):render(src,work/'phase0'/src.name)
    for phase in ['phase0','phase1','phase2','phase3','phase4','phase5','phase6','phase3/icm_revision4']:
        for directory in ['_restricted','_logs','figures','qa','analysis']:(work/phase/directory).mkdir(parents=True,exist_ok=True)
    print('Private workspace source templates prepared. No data were read and no analysis was run.')

if __name__=='__main__':main()
