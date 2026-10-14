#!/usr/bin/env python3
"""Run optional numerical search, separately from the supplied exact proofs."""
from pathlib import Path
from datetime import datetime
import argparse
import os
import subprocess
import sys
import tempfile
from repository_layout import ROOT, stage

CASES=('351_140_73_44','375_102_45_21','405_132_63_33','441_88_35_13','476_133_60_28','540_147_66_30','550_162_75_36','575_112_45_16','703_182_81_35','1344_221_88_26')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases',nargs='+',choices=CASES,default=list(CASES),help='Full v_k_lambda_mu parameter tuples')
    parser.add_argument('--output',type=Path,help='New output directory; defaults under results/runs/')
    parser.add_argument('--check-inputs',action='store_true',help='Check syntax and inputs, without solving')
    parser.add_argument('--prepare-only',action='store_true',help='Stage a new run without solving')
    args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='srg_discovery_') as tmp:
        sandbox=stage(Path(tmp)/'repository')
        cmd=[sys.executable,str(sandbox/'discovery/run.py'),'--cases',*[c.split('_')[0] for c in args.cases]]
        if args.check_inputs:cmd+=['--check-inputs']
        else:
            output=(args.output or ROOT/'results/runs'/('discovery_'+datetime.now().strftime('%Y%m%d_%H%M%S'))).resolve()
            cmd+=['--output',str(output)]
            if args.prepare_only:cmd+=['--prepare-only']
        run=subprocess.run(cmd,cwd=sandbox,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True)
        text=(run.stdout+run.stderr).replace(str(sandbox)+os.sep,'').replace('discovery/requirements.txt','code/discovery/requirements.txt')
        print(text,end='');raise SystemExit(run.returncode)

if __name__=='__main__':main()
