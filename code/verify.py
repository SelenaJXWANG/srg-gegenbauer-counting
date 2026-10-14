#!/usr/bin/env python3
"""Check all supplied proofs and their actual TeX tables, using exact arithmetic.

Python 3.9+ standard library only. No SDP solver or internet is required.
Run from any directory: python3 -S /path/to/repository/code/verify.py
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from repository_layout import ROOT, stage, input_hashes

STAGES=[
 ('ten_parameter_certificates','data/legacy_ten/run_all.py'),
 ('ordinary_counting_405_and_all_degree_comparison','data/pure405/run_all.py'),
 ('printed_short_certificates','verification/verify_displayed_tables.py'),
 ('printed_expanded_certificates','verification/verify_expanded_certificates.py'),
 ('article_parameter_and_structural_alignment','verification/verify_article_metadata.py'),
]

def clean(value,sandbox):
    if isinstance(value,str):
        return value.replace(str(sandbox)+os.sep,'').replace(str(ROOT)+os.sep,'')
    if isinstance(value,list):return [clean(x,sandbox) for x in value]
    if isinstance(value,dict):return {clean(k,sandbox):clean(v,sandbox) for k,v in value.items()}
    return value

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reports',type=Path,default=ROOT/'results/latest_verification')
    parser.add_argument('--keep-logs',action='store_true')
    args=parser.parse_args()
    if sys.flags.optimize:raise SystemExit('Run without -O: the exact engines use assertions.')
    if sys.version_info<(3,9):raise SystemExit('Python 3.9 or later is required.')
    args.reports.mkdir(parents=True,exist_ok=True)
    before=input_hashes();started=time.monotonic()
    report={'result':'RUNNING','started_utc':datetime.now(timezone.utc).isoformat(),
            'python':sys.version.split()[0],'arithmetic':'fractions.Fraction; Python standard library',
            'fixed_inputs_sha256':before,'stages':[]}
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
    with tempfile.TemporaryDirectory(prefix='srg_exact_') as tmp:
        sandbox=stage(Path(tmp)/'repository')
        for i,(label,rel) in enumerate(STAGES,1):
            print('[%s/%s] %s'%(i,len(STAGES),label),flush=True)
            then=time.monotonic();p=sandbox/rel
            run=subprocess.run([sys.executable,'-B','-S',str(p)],cwd=p.parent,env=env,capture_output=True,text=True)
            output=clean(run.stdout+run.stderr,sandbox)
            item={'check':label,'exit_code':run.returncode,'elapsed_seconds':round(time.monotonic()-then,3),'output_sha256':sha256(output.encode()).hexdigest()}
            report['stages'].append(item)
            if args.keep_logs or run.returncode:(args.reports/(label+'.log')).write_text(output)
            if run.returncode:
                report['result']='FAIL';report['failed_check']=label
                (args.reports/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
                print(output[-16000:]);raise SystemExit('Exact verification failed: '+label)
            print('  PASS (%.3f s)'%item['elapsed_seconds'],flush=True)
        reports={
          'data/legacy_ten/last_three_local_replay.json':'structural_lemmas.json',
          'data/legacy_ten/verification_1344_reconstructed.json':'1344_221_88_26_exact_reconstruction.json',
          'data/pure405/405_sevenpoint_core_independent_audit_20261003.json':'405_132_63_33_counting_model.json',
          'data/pure405/405_kernel_all_degree_verify_20261003.json':'405_132_63_33_all_degree_comparison.json',
          'data/pure405/data/degree2/pure_exact_replay.json':'405_132_63_33_exact_certificate.json',
        }
        cases={375:'375_102_45_21',441:'441_88_35_13',575:'575_112_45_16'}
        for v,case in cases.items():
            for source,label in [('independent_psd_completion_audit_%s.json'%v,'kernel_completion'),('moments_certificate_audit_%s.json'%v,'integer_moments_and_margin'),('rebuild_%s/affine_source_audit.json'%v,'counting_equations')]:
                reports['data/legacy_ten/'+source]=case+'_'+label+'.json'
        details=args.reports/'details';details.mkdir(exist_ok=True)
        report['generated_reports']=[]
        for source,name in reports.items():
            p=sandbox/source
            if not p.is_file():raise RuntimeError('Expected verification result missing: '+source)
            (details/name).write_text(json.dumps(clean(json.loads(p.read_text()),sandbox),indent=2)+'\n')
            report['generated_reports'].append('details/'+name)
    report['fixed_inputs_unchanged']=(input_hashes()==before)
    if not report['fixed_inputs_unchanged']:raise RuntimeError('A fixed input changed during verification.')
    report['result']='PASS';report['elapsed_seconds']=round(time.monotonic()-started,3)
    report['meaning']='All ten exclusions, the ordinary-counting seven-point 405 proof, its all-degree five-point comparison, and the specified actual TeX certificate/parameter/structural tables passed exact checks.'
    (args.reports/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: five exact verification stages; fixed inputs unchanged.',flush=True)
    print('Report: '+str(args.reports/'summary.json'))

if __name__=='__main__':main()
