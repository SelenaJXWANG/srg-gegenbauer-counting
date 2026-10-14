#!/usr/bin/env python3
"""Public command line for the article's exact proofs and optional search."""
import argparse
import subprocess
import sys
from pathlib import Path

def main():
    parser=argparse.ArgumentParser(description=__doc__,epilog='Use `verify --help` or `discover --help` for options.')
    parser.add_argument('command',choices=('verify','discover'))
    args=parser.parse_args(sys.argv[1:2])
    remaining=sys.argv[2:]
    script='verify.py' if args.command=='verify' else 'run_discovery.py'
    command=[sys.executable,'-B']
    if args.command=='verify':command+=['-S']
    raise SystemExit(subprocess.call([*command,str(Path(__file__).with_name(script)),*remaining]))

if __name__=='__main__':main()
