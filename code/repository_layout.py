"""Resolve public repository paths and stage byte-preserved research engines.

The public filenames are descriptive. Audited engines retain their original
imports and relative paths; a temporary layout adapter keeps that renaming
independent of mathematical computations. See extras/file_mapping.json.
"""
from pathlib import Path
from hashlib import sha256
import json
import shutil

ROOT = Path(__file__).resolve().parents[1]

def mapping():
    return json.loads((ROOT/'extras/file_mapping.json').read_text())['entries']

def checked_entries():
    entries = mapping()
    for entry in entries:
        path = ROOT/entry['path']
        if not path.is_file() or sha256(path.read_bytes()).hexdigest() != entry['sha256']:
            raise RuntimeError('Fixed input changed or missing: '+entry['path'])
    return entries

def stage(destination):
    """Copy pinned inputs into an isolated compatibility tree, without edits."""
    destination = Path(destination)
    entries = checked_entries()
    for entry in entries:
        target = destination/entry['runtime_path']
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/entry['path'], target)
    # Old manifests included obsolete prose/reports. These manifests enumerate
    # only the byte-preserved inputs used by the retained mathematical checks.
    for rel, name in [('data/legacy_ten/six_exact_data','MANIFEST.sha256'),
                      ('data/legacy_ten','INPUTS.sha256'),
                      ('data/pure405','INPUTS.sha256')]:
        folder=destination/rel
        paths=sorted(p for p in folder.rglob('*') if p.is_file() and p.name!=name)
        (folder/name).write_text(''.join(sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(folder))+'\n' for p in paths))
    return destination

def input_hashes():
    return {e['path']:e['sha256'] for e in checked_entries()}
