"""Fresh full-cone SDP discovery, rational reconstruction and exact verification."""
from pathlib import Path
import json,sys,time,traceback
from build_model import build
from solve_sdp import run
from recover_rational import main as recover
from verify_certificate import verify
ROOT=Path(__file__).resolve().parents[1]
CASES={540:(540,147,66,30),550:(550,162,75,36),703:(703,182,81,35),1344:(1344,221,88,26),351:(351,140,73,44),476:(476,133,60,28)}
results=[]
for v in map(int,sys.argv[1:] or CASES):
    started=time.monotonic();out=ROOT/'fresh'/str(v);p=6 if v in (351,476) else 5
    try:
        print('START',v,flush=True)
        data=build(CASES[v],out,p=p)
        run(out/'data.json',out/'numerical_solution.json',variable_scale=10**(10 if p==5 else 12))
        recover(out)
        certificate=json.loads((out/'certificate.json').read_text())
        audit=verify(json.loads((out/'data.json').read_text()),certificate)
        (out/'verification.json').write_text(json.dumps(audit,indent=2)+'\n')
        result=dict(v=v,status='EXACT_VERIFIED',kernel_terms=len(certificate['terms']),count_terms=sum(t!='0' for t in certificate['z']),p=p,seconds=time.monotonic()-started)
    except Exception as e:
        traceback.print_exc();result=dict(v=v,status='FAILED',error=str(e),seconds=time.monotonic()-started)
    results.append(result);print('FINAL_CASE',json.dumps(result),flush=True)
    (ROOT/'fresh'/f'batch_{sys.argv[1] if len(sys.argv)>1 else "all"}.json').write_text(json.dumps(results,indent=2)+'\n')
if any(result['status']!='EXACT_VERIFIED' for result in results):
    raise SystemExit('One or more fresh cases failed; see batch report and log.')
