#!/usr/bin/env python3
"""Byte-for-byte synthetic stage-artifact parity: notebook cells vs CLI.
UNRUN when Rscript is not present. No real-data access or fitting.
"""
import json, os, pathlib, shutil, subprocess, tempfile, hashlib, sys
ROOT=pathlib.Path(__file__).resolve().parents[1]

def main():
 if not shutil.which('Rscript'):
  print('UNRUN: Rscript absent; notebook/CLI byte test NOT verified')
  return 2
 nb=json.loads((ROOT/'notebooks/01_explore_one_galaxy.ipynb').read_text())
 with tempfile.TemporaryDirectory() as temp:
  # Resolve macOS's /var alias before the output safety check sees parents.
  temp=pathlib.Path(temp).resolve()
  a=temp/'notebook'; b=temp/'cli'
  a.mkdir();b.mkdir()
  # Notebook invocation uses exactly the published code cells, until G2.
  cells=[c for c in nb['cells'] if c['cell_type']=='code'][:4]
  text='\n'.join(''.join(c['source']) for c in cells)
  # Replace the notebook's output folder only; no science call is modified.
  text=text.replace('ctx <- bd_context(id, cfg, run_id, mode)',
                    'ctx <- bd_context(id, cfg, run_id, mode, output_override=Sys.getenv("BD_OUTPUT_OVERRIDE"))')
  nbscript=pathlib.Path(temp)/'notebook_run.R';nbscript.write_text(text)
  env=dict(os.environ,BD_PROFILE='mac',BD_ID='W42710765002459',
           BD_RUN_ID='byte_equiv',BD_MODE='synthetic',BD_OUTPUT_OVERRIDE=str(a))
  n=subprocess.run(['Rscript',str(nbscript)],cwd=ROOT,env=env,capture_output=True,text=True)
  if n.returncode:raise RuntimeError('Notebook test invocation failed: '+n.stderr)
  # CLI calls bd_run_one() which delegates the exact same bd_stage() stages.
  # Using an output root override only for this isolated synthetic test.
  env['BD_OUTPUT_OVERRIDE']=str(b)
  c=subprocess.run(['Rscript','scripts/run_one_galaxy.R','--profile=mac',
       '--id=W42710765002459','--mode=synthetic','--run-id=byte_equiv',
       '--stages=inputs,g2'],cwd=ROOT,env=env,capture_output=True,text=True)
  if c.returncode:raise RuntimeError('CLI parity invocation failed: '+c.stderr)
  for name in ['inputs.rds','g2.rds']:
   p=a/'byte_equiv'/'W42710765002459'/name
   q=b/'byte_equiv'/'W42710765002459'/name
   if not p.exists() or not q.exists(): raise AssertionError('Missing artifacts: '+name)
   if p.read_bytes()!=q.read_bytes():raise AssertionError('NOT BYTE IDENTICAL: '+name)
  print('PASS notebook vs CLI BYTE-IDENTICAL synthetic inputs.rds + g2.rds')
 return 0
if __name__=='__main__':sys.exit(main())
