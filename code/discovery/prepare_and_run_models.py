#!/usr/bin/env python3
"""Optional numerical discovery in a new output directory.

The six five/six-point models are rebuilt from parameters. The four
seven-point models use retained exact graph/count/structural source data.
No certificate direction or multiplier is supplied to the fresh solve.
This wrapper does not modify the checked paper certificates.
"""
from pathlib import Path
from hashlib import sha256
from importlib import metadata
from datetime import datetime
import argparse
import ast
import json
import os
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
CASES=(351,375,405,441,476,540,550,575,703,1344)
SIX=(540,550,703,1344,351,476)
SEVEN=(405,375,441,575)

def inputs():
    pinned = 0
    for package in ("legacy_ten", "pure405"):
        manifest = ROOT/"data"/package/"INPUTS.sha256"
        for line in manifest.read_text().splitlines():
            digest, name = line.split("  ", 1)
            candidate = manifest.parent/name
            if not candidate.is_file() or sha256(candidate.read_bytes()).hexdigest() != digest:
                raise SystemExit("Pinned model input changed: data/"+package+"/"+name)
            pinned += 1
    paths=[ROOT/"data/pure405/atlas/graphs.txt"]
    paths += [ROOT/("data/pure405/atlas/map%d.bin"%n) for n in range(1,8)]
    paths += [ROOT/"data/pure405/data/degree2/core_recipes_405.json"]
    for v in (375,441,575):
        paths += [ROOT/("data/legacy_ten/rebuild_%d/%s"%(v,name)) for name in
                  ("affine_recipes_%d.json"%v,"forbidden_%d.json"%v,"integer_moments.json")]
    for p in paths:
        if not p.is_file():
            raise SystemExit("Missing discovery model input: "+str(p.relative_to(ROOT)))
    programs=sorted((HERE/"engine").glob("*.py"))
    for p in programs:
        ast.parse(p.read_text(),filename=p.name)
    return {"pinned_inputs_checked":pinned, "model_inputs":len(paths), "programs_syntax_checked":len(programs),
            "model_sha256":{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths}}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases",nargs="+",type=int,choices=CASES,default=list(CASES))
    parser.add_argument("--output",type=Path,help="New directory for staged programs, model inputs and fresh results")
    parser.add_argument("--check-inputs",action="store_true",help="Check syntax, model paths and installed dependency metadata, without solving")
    parser.add_argument("--prepare-only",action="store_true",help="Stage a new run directory without importing optimization packages")
    args=parser.parse_args()
    report=inputs();versions={}
    for package in ("cvxpy","clarabel","numpy","scipy"):
        try: versions[package]=metadata.version(package)
        except metadata.PackageNotFoundError: versions[package]=None
    report["installed_versions"]=versions
    if args.check_inputs:
        report["check"]="PASS syntax and model input paths; no optimization executed"
        print(json.dumps(report,indent=2));return
    if not args.prepare_only and any(v is None for v in versions.values()):
        raise SystemExit("Install discovery/requirements.txt in your numerical Python environment first, or use --prepare-only.")
    destination=(args.output or ROOT/"runs"/("discovery_"+datetime.now().strftime("%Y%m%d_%H%M%S"))).resolve()
    if destination.exists():
        raise SystemExit("Choose a new output directory; existing results will not be overwritten: "+str(destination))
    destination.mkdir(parents=True)
    shutil.copytree(HERE/"engine",destination/"engine",ignore=shutil.ignore_patterns("__pycache__","*.pyc"))
    shutil.copytree(ROOT/"data",destination/"article_replay",ignore=shutil.ignore_patterns("__pycache__","*.pyc"))
    (destination/"fresh").mkdir()
    report.update({"requested_cases":args.cases,"optimization_executed":False})
    commands=[]
    selected=[str(v) for v in SIX if v in args.cases]
    if selected:commands.append([sys.executable,str(destination/"engine/run_fresh.py"),*selected])
    for v in SEVEN:
        if v in args.cases:
            commands += [[sys.executable,str(destination/"engine/fresh_seven.py"),str(v)],
                         [sys.executable,"-S",str(destination/"engine/verify_fresh_seven.py"),str(v)]]
    if set(args.cases)==set(CASES):
        commands.append([sys.executable,"-S",str(destination/"engine/export_catalogue.py")])
    report["commands"]=[[Path(c[0]).name]+[str(Path(a).relative_to(destination)) if a.startswith(str(destination)) else a for a in c[1:]] for c in commands]
    (destination/"run_configuration.json").write_text(json.dumps(report,indent=2)+"\n")
    if args.prepare_only:
        print("Prepared discovery inputs at "+str(destination));return
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE="1")
    for command in commands:
        subprocess.run(command,cwd=destination,env=env,check=True)
    report["optimization_executed"]=True
    report["result"]="Completed numerical search and the associated exact checks."
    (destination/"run_configuration.json").write_text(json.dumps(report,indent=2)+"\n")
    print("Finished. Fresh results: "+str(destination/"fresh"))

if __name__=="__main__":
    main()
