#!/usr/bin/env python3
"""Install local Git data-use hooks without modifying repository history."""
from pathlib import Path
import argparse,subprocess,json
def main():
 a=argparse.ArgumentParser();a.add_argument('--repo',default='.');args=a.parse_args();root=Path(args.repo).resolve()
 if not (root/'tools/small_cell_gate.py').is_file():raise RuntimeError('gate_not_found')
 h=root/'.githooks';h.mkdir(exist_ok=True)
 scripts={
  'pre-commit':'#!/bin/sh\nset -eu\ncd "$(git rev-parse --show-toplevel)"\npython tools/small_cell_gate.py --staged --repo . --report "$(git rev-parse --git-path small_cell_gate/pre_commit.json)"\n',
  'pre-push':'#!/bin/sh\nset -eu\ncd "$(git rev-parse --show-toplevel)"\npython tools/small_cell_gate.py --pre-push --repo . --report "$(git rev-parse --git-path small_cell_gate/pre_push.json)"\n'
 }
 for name,body in scripts.items():
  target=h/name;target.write_text(body,encoding='utf-8',newline='\n');target.chmod(0o755)
 subprocess.run(['git','-C',str(root),'config','core.hooksPath','.githooks'],check=True,capture_output=True)
 print(json.dumps({'status':'INSTALLED','repo':str(root),'hooks':sorted(scripts)}))
if __name__=='__main__':main()
