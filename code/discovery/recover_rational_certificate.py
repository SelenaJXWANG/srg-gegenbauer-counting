"""Article 5.3: round NEW SDP eigendirections, solve LP, reconstruct exactly.

The LP selects a basis only. Its coefficients are solved again over Q and
accepted only if every multiplier is nonnegative and the identity is exact.
"""
from pathlib import Path
import sys,json,math,time,argparse
from fractions import Fraction as Q
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'python'))
import numpy as np
from scipy.optimize import linprog
from build_model import rref,serial

def primitive(f,den):
    a=[int(round(float(x)*den/max(abs(t) for t in f))) for x in f]
    g=math.gcd(*a)
    if not g:return None
    a=[x//g for x in a]
    if next(x for x in a if x)!=abs(next(x for x in a if x)):
        a=[-x for x in a]
    return a

def main(folder):
    start=time.monotonic()
    data=json.loads((folder/'data.json').read_text())
    numerical=json.loads((folder/'numerical_solution.json').read_text())
    blocks={(b['H'],b['ell']):[[[Q(t) for t in row] for row in mat] for mat in b['M']] for b in data['blocks']}
    affine=[[Q(t) for t in row] for row in data['affine']]
    q=len(data['free_indices']); target=[Q(-1)]+[Q(0)]*q
    eigs=[e for e in numerical['directions'] if e['eigenvalue']>1e-6]
    candidates=[];seen=set();trials=[]
    for den in [1,2,3,4,5,6,8,10,12,16,20,30,50,100,200,500,1000,10000]:
        for eig in eigs:
            f=primitive(eig['f'],den)
            if f is None:continue
            key=(eig['H'],eig['ell'],tuple(f))
            if key in seen:continue
            seen.add(key)
            M=blocks[key[:2]]
            kap=[sum(f[a]*f[b]*mm[a][b] for a in range(len(f)) for b in range(len(f))) for mm in M]
            candidates.append(dict(H=key[0],ell=key[1],f=f,kappa=kap,rounding_denominator=den))
        columns=[c['kappa'] for c in candidates]+affine
        varscale=[Q(1)]+[Q(str(numerical['variable_scaling']['variable_scale']))]*q
        scales=[max(abs(x*s) for x,s in zip(col,varscale)) or Q(1) for col in columns]
        B=np.array([[float(x*s/z) for x,s in zip(col,varscale)] for col,z in zip(columns,scales)]).T
        res=linprog(np.ones(len(columns)),A_eq=B,b_eq=np.array([-1.]+[0.]*q),bounds=(0,None),method='highs',
                    options=dict(primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))
        trial=dict(denominator=den,candidate_count=len(candidates),status=int(res.status),message=res.message)
        trials.append(trial);print(trial,flush=True)
        if not res.success:continue
        support=[i for i,x in enumerate(res.x) if x>1e-8]
        # Re-solve the original, unscaled coefficient identity with rational Gaussian elimination.
        aug=[[columns[i][j] for i in support]+[target[j]] for j in range(q+1)]
        red,piv=rref(aug,len(support))
        if len(piv)!=len(support) or any(row[-1] for row in red[len(piv):]):
            trial['exact_failure']='basis_rank_or_consistency';continue
        weights=[Q(0)]*len(support)
        for row,i in zip(red,piv): weights[i]=row[-1]
        if min(weights)<0:
            trial['exact_failure']='negative_weight';continue
        total=[sum(w*columns[i][j] for i,w in zip(support,weights)) for j in range(q+1)]
        assert total==target
        terms=[];z=[Q(0)]*len(affine)
        for i,w in zip(support,weights):
            if not w:continue
            if i<len(candidates): terms.append(dict(candidates[i],eta=w))
            else:z[i-len(candidates)]=w
        cert=dict(parameters=data['parameters'],terms=terms,z=z,identity=total,
                  source='Fresh SDP eigendirections; no prior certificate f or weights read',
                  exact_reconstruction='Rational Gaussian elimination on LP-selected columns',
                  numerical_file='numerical_solution.json',trials=trials,seconds=time.monotonic()-start)
        (folder/'certificate.json').write_text(json.dumps(serial(cert),indent=2)+'\n')
        print('EXACT CERTIFICATE',len(terms),'kernel terms',sum(t!=0 for t in z),'count terms',flush=True)
        for t in terms:print(t['H'],t['ell'],t['f'],'eta',str(t['eta']),flush=True)
        return
    (folder/'rational_recovery_failed.json').write_text(json.dumps(trials,indent=2)+'\n')
    raise RuntimeError('No exact certificate recovered from rounded fresh directions')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--folder',type=Path,default=HERE);args=ap.parse_args();main(args.folder)
