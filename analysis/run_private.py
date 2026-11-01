"""Run one prepared private-workspace script after checking its output boundary."""
from pathlib import Path
import argparse
import os
import subprocess
import sys

REPO=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--workspace',type=Path)
    parser.add_argument('--script',required=True,help='Prepared script path relative to workspace/sql')
    parser.add_argument('--database',help='Local PostgreSQL database for SQL extraction')
    parser.add_argument('--rscript',default='Rscript',help='Installed Rscript executable')
    parser.add_argument('--check',action='store_true',help='Validate the boundary without executing')
    parser.add_argument('script_args',nargs=argparse.REMAINDER,help='Arguments passed unchanged to Python/R after --')
    args=parser.parse_args()
    supplied=str(args.workspace) if args.workspace else os.environ.get('CUFF_WORKSPACE','')
    if not supplied or '@' in supplied or not Path(supplied).is_absolute():parser.error('Set an absolute prepared private workspace')
    work=Path(supplied).resolve()
    if work==REPO or REPO in work.parents:parser.error('Workspace must be outside the repository checkout')
    path=(work/'sql'/args.script).resolve()
    if not path.is_relative_to(work/'sql'):parser.error('Script must be inside the prepared workspace/sql directory')
    text=path.read_text(encoding='utf-8')
    if any(t in text for t in ['@CUFF_WORKSPACE@','@AMSTERDAM_ROOT@','@SICDB_ROOT@']):parser.error('Unrendered template; prepare the workspace first')
    if path.suffix.lower()=='.sql':
        if not args.database:parser.error('Specify --database for PostgreSQL extraction')
        command=['psql','-X','-v','ON_ERROR_STOP=1','-v','CUFF_PREPARED_WORKSPACE=1','-d',args.database,'-f',str(path)]
    elif path.suffix.lower()=='.py':command=[sys.executable,str(path)]
    elif path.suffix.lower()=='.r':command=[args.rscript,str(path)]
    else:parser.error('Only Python, R or PostgreSQL extraction scripts are supported')
    forwarded=args.script_args[1:] if args.script_args[:1]==['--'] else args.script_args
    if path.suffix.lower()=='.sql' and forwarded:parser.error('SQL options belong to the launcher; argument forwarding is for Python/R scripts')
    command+=forwarded
    if args.check:print('Prepared script boundary verified; no data read or analysis executed.');return
    subprocess.run(command,cwd=work,check=True)

if __name__=='__main__':main()
