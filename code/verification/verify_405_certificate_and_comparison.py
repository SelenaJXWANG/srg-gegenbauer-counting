#!/usr/bin/env python3
"""Replay the portable pure405 proof and its five-point comparison."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent

def main():
    started = time.monotonic()
    files = 0
    for line in (HERE/'INPUTS.sha256').read_text().splitlines():
        digest, name = line.split('  ', 1)
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == digest, name
        files += 1
    print('PASS input manifest:', files, 'files', flush=True)
    scripts = [
        '405_sevenpoint_core_independent_audit_20261003.py',
        '405_kernel_certificate_verify_20261003.py',
        '405_kernel_all_degree_verify_20261003.py',
    ]
    records = []
    for script in scripts:
        start = time.monotonic()
        result = subprocess.run([sys.executable, '-S', str(HERE/script)],
                                cwd=HERE, capture_output=True, text=True)
        print(result.stdout, end='', flush=True)
        if result.stderr:
            print(result.stderr, file=sys.stderr, end='', flush=True)
        (HERE/(script+'.log')).write_text(result.stdout+result.stderr)
        if result.returncode:
            raise RuntimeError('Verification failed: '+script)
        records.append({'script': script, 'exit_code': result.returncode,
                        'elapsed_seconds': round(time.monotonic()-start, 3)})
    record = {'status': 'PASS', 'input_files': files, 'checks': records,
              'elapsed_seconds': round(time.monotonic()-started, 3),
              'meaning': 'Pure seven-point nonexistence certificate and exact all-degree five-point relaxation witness both verified.'}
    (HERE/'verification_summary.json').write_text(json.dumps(record, indent=2)+'\n')
    print('PASS: pure405 proof and five-point comparison; all source layers verified.', flush=True)

if __name__ == '__main__':
    main()
