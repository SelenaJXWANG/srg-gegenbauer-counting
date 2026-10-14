"""Offline exact checks of NEW seven-point certificates; standard library only.

Uses the article's independent source verifiers, with fresh certificates and
the already audited identical count/source data. No discovery module or solver
is imported. Complete verification also runs article_replay/run_all.py.
"""
from pathlib import Path
from fractions import Fraction as F
import importlib.util,json,sys,shutil
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)

def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def verify(v):
    dest=ROOT/'fresh'/str(v);record=json.loads((dest/'certificate.json').read_text())
    assert record['parameters'][0]==v and F(record['exact_margin'])>0
    if v==405:
        src=ROOT/'article_replay/pure405'
        recipes=src/'data/degree2/core_recipes_405.json'
        shutil.copy2(recipes,dest/recipes.name)
        converted=dict(record,combined_core_coefficients=record['combined_coefficients'],core_density_upper_bounds=record['bounds'])
        (dest/'pure_exact_certificate.json').write_text(json.dumps(converted,indent=2)+'\n')
        mod=load('pure_verify',src/'405_kernel_certificate_verify_20261003.py')
        mod.ROOT=src;mod.OUT=dest
        mod.main()
        result=json.loads((dest/'pure_exact_replay.json').read_text())
        assert result['status']=='PASS'
    else:
        src=ROOT/'article_replay/legacy_ten';model=src/f'rebuild_{v}'
        rec={int(i):{int(j):F(q) for j,q in row.items()} for i,row in json.loads((model/f'affine_recipes_{v}.json').read_text()).items()}
        dim=len(record['combined_coefficients']);rows=[]
        for item in record['terms']:
            full={int(i):F(q) for i,q in item['full_row'].items()};r=[F(0)]*dim
            for i,q in full.items():
                for j,p in rec[i].items():r[j]+=q*p
            scale=F(item['scale']);assert scale==max(map(abs,r)) and scale>0
            rows.append(dict(source=item['source'],weight=item['weight'],normalization=item['scale'],full_y_row=item['full_row'],reduced_normalized_row=[str(q/scale) for q in r]))
        converted=dict(v=v,dimension=dim,variable_upper_bound=1024 if v==375 else 256,rows=rows,active_row_count=len(rows),combined_coefficients=record['combined_coefficients'],exact_margin=record['exact_margin'])
        # Place only the verified model inputs beside the fresh certificate.
        local=dest/f'rebuild_{v}';local.mkdir(exist_ok=True)
        for name in [f'forbidden_{v}.json','completion_constraints.json','integer_moments.json',f'affine_recipes_{v}.json',f'exact_elimination_{v}.json']:
            shutil.copy2(model/name,local/name)
        (local/'recovered_exact_certificate.json').write_text(json.dumps(converted,indent=2)+'\n')
        psd=load('psd_verify',src/'verify_psd_completion.py');psd.DATA=dest;psd.CATALOGUE=src/'atlas_graphs.txt'
        psd.main(v)
        arithmetic=load('moment_verify',src/'verify_moments_and_certificate.py');arithmetic.DATA=dest;arithmetic.CATALOGUE=src/'atlas_graphs.txt'
        result=arithmetic.main(v)
        assert result['status']=='PASS'
    (dest/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print('FRESH_SEVEN_VERIFIED',v,flush=True)

if __name__=='__main__':
    for v in map(int,sys.argv[1:] or [405,375,441,575]):verify(v)
