#!/usr/bin/env python3
"""Replay every certificate and all reconstructed source layers, offline."""
from pathlib import Path
from hashlib import sha256
import subprocess,sys,time

ROOT=Path(__file__).resolve().parent


def main():
    started=time.time();checked=0
    for line in (ROOT/'INPUTS.sha256').read_text().splitlines():
        expected,name=line.split('  ',1)
        assert sha256((ROOT/name).read_bytes()).hexdigest()==expected, 'Input changed: '+name
        checked+=1
    print('Supplement input manifest verified:',checked,'files',flush=True)
    commands=[
        ['six_exact_data/verify_all.py'],
        ['verify_1344_reconstructed.py'],
        ['verify_local_ingredients.py'],
        ['verify_affine_recipes.py','375','441','575'],
        ['verify_psd_completion.py','375'],
        ['verify_psd_completion.py','441'],
        ['verify_psd_completion.py','575'],
        ['verify_moments_and_certificate.py','375','441','575'],
    ]
    for command in commands:
        print('\nRUN',*command,flush=True)
        subprocess.run([sys.executable,'-S',str(ROOT/command[0]),*command[1:]],cwd=ROOT,check=True)
    print('\nPASS: all ten parameter cases, including every reconstructed source layer.',flush=True)
    print('Elapsed seconds:',round(time.time()-started,3),flush=True)


if __name__=='__main__':main()
