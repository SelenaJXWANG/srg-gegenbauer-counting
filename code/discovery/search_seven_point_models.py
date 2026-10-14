"""Fresh seven-point cone SDP; no saved f, cut, dual or certificate is read.

Verified archived count eliminations and proved structural source rows are
used as model data. New cone dual eigendirections are rounded, an LP selects
rows, and every selected row is regenerated over the rationals. Bounds pay
ALL residual coefficients exactly; no floating residual is ignored.
"""
from pathlib import Path
import sys,json,time,importlib.util,itertools,math,warnings
from fractions import Fraction as F
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;BASE=ROOT.parents[1]
sys.path.insert(0,str(HERE/'python'))
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog
import cvxpy as cp
from row_sources import ExactModel

def setup(v):
    return ExactModel(v)

def numeric_blocks(eng,P):
    blocks=[]
    for hid,A in enumerate(eng.aa):
        m=len(A)
        if m>5:break
        if hid in eng.forbidden:continue
        eta=float(eng.eta);G=A+3*np.eye(m)-eta
        if eng.v==405:G=G/(3-eta)
        if np.linalg.eigvalsh(G)[0]<=1e-10:continue
        pats=np.array(list(itertools.product((0,1),repeat=m)))
        u=(pats-eta)/(3-eta) if eng.v==405 else pats-eta
        proj=u@np.linalg.inv(G)@u.T
        rho=(1 if eng.v==405 else 3-eta)-np.diag(proj)
        singles=[eng.extended(A,p) for p in pats]
        keep=list(range(len(pats))) if eng.v==405 else [i for i,h in enumerate(rho) if h>=-1e-10 and singles[i] not in eng.forbidden]
        if not keep:continue
        pats=pats[keep];rho=rho[keep];proj=proj[np.ix_(keep,keep)];singles=[singles[i] for i in keep]
        n=len(pats);D=eng.f-m
        rr=[];cc=[];values=[[],[],[]]
        for i,p in enumerate(pats):
            rr.append(i*n+i);cc.append(singles[i])
            for ell in range(3):values[ell].append(rho[i]**ell*(1-1/D if ell==2 and eng.v!=405 else 1))
            for j,q in enumerate(pats):
                for edge in (0,1):
                    target=eng.extended(A,p,q,edge)
                    if target in eng.forbidden:continue
                    t=((edge-eta)/(3-eta) if eng.v==405 else edge-eta)-proj[i,j]
                    vals=[1,t,(D*t*t-rho[i]*rho[j])/(D-1) if eng.v==405 else t*t-rho[i]*rho[j]/D]
                    rr.append(i*n+j);cc.append(target)
                    for ell in range(3):values[ell].append((eng.v-m-1)/2**(m+1)*vals[ell])
        for ell in range(3):
            K=sp.csr_matrix((values[ell],(rr,cc)),shape=(n*n,len(eng.aa)))@P
            M=np.array(K).reshape(n,n,eng.nvar).transpose(2,0,1)
            M=(M+M.transpose(0,2,1))/2
            scale=float(np.max(abs(M)))
            if scale<1e-8:continue
            blocks.append(dict(anchor=hid,degree=ell,m=m,size=n,M=M,scale=scale,profiles=pats.tolist()))
    return blocks

def run(v):
    started=time.monotonic();out=ROOT/'fresh'/str(v);out.mkdir(parents=True,exist_ok=True)
    eng=setup(v)
    P=np.zeros((len(eng.aa),eng.nvar))
    for i,row in eng.recipes.items():
        for j,q in row.items():P[i,j]=float(q)
    blocks=numeric_blocks(eng,P)
    print('BUILT',v,len(blocks),'blocks',eng.nvar,'coordinates',time.monotonic()-started,flush=True)
    counts=[i for i in range(len(P)) if i not in eng.forbidden and max(abs(P[i]))>1e-14]
    base=[-P[i] for i in counts];sources=[dict(kind='nonnegative_count',graph_id=i) for i in counts]
    for i,row in enumerate(eng.moments['ineq']):
        exact=eng.reduce({int(k):F(c) for k,c in row.items()})
        if not any(exact):continue
        base.append(np.array([float(q) for q in exact]))
        sources.append(dict(kind='integer_moment',row_index=i))
    base=np.array(base);base/=np.max(abs(base),axis=1)[:,None]
    if '--recover-only' in sys.argv:
        directions=json.loads((out/'numerical_solution.json').read_text())['directions']
    else:
        x=cp.Variable(eng.nvar-1);t=cp.Variable(nonneg=True);xx=cp.hstack([1,x])
        cons=[base@xx<=t];matrix_cons=[]
        for block in blocks:
            M=block['M']/block['scale'];n=block['size']
            mm=cp.reshape(M.reshape(eng.nvar,n*n).T@xx,(n,n),order='C')
            c=mm+t*np.eye(n)>>0;cons.append(c);matrix_cons.append(c)
        prob=cp.Problem(cp.Minimize(t),cons)
        prob.solve(solver='CLARABEL',verbose=True,tol_gap_abs=1e-11,tol_feas=1e-11,tol_gap_rel=1e-11,max_iter=400,time_limit=180)
        print('SDP',v,prob.status,prob.value,flush=True)
        assert t.value is not None and float(t.value)>1e-10
        directions=[]
        for block,c in zip(blocks,matrix_cons):
            Y=(c.dual_value+c.dual_value.T)/2;ev,U=np.linalg.eigh(Y)
            for j in range(len(ev)):
                if ev[j]>1e-9:
                    f=np.rint(U[:,j]/np.max(abs(U[:,j]))*1000000).astype(np.int64)
                    g=math.gcd(*map(int,f));f=f//g
                    directions.append(dict(kind='PSD',anchor=block['anchor'],degree=block['degree'],direction=f.tolist(),eigenvalue=float(ev[j])))
        (out/'numerical_solution.json').write_text(json.dumps(dict(v=v,status=prob.status,phase_one_margin=float(t.value),blocks=len(blocks),coordinates=eng.nvar,directions=directions,provenance='Fresh cone SDP on verified count/source data; no archived directions, cuts or weights read'),indent=2)+'\n')
    print('DIRECTIONS',v,len(directions),flush=True)
    bykey={(b['anchor'],b['degree']):b for b in blocks}
    rows=list(base);all_sources=list(sources)
    for source in directions:
        b=bykey[source['anchor'],source['degree']];w=np.array(source['direction'],dtype=float)
        row=-np.einsum('a,jab,b->j',w,b['M'],w);scale=max(abs(row))
        if scale>1e-15:rows.append(row/scale);all_sources.append(source)
    rows=np.array(rows)
    # Replace all selected floating rows by exact reconstructions before accepting.
    cache={};checks=0
    for retry in range(40):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            lp=linprog(np.r_[np.zeros(eng.nvar-1),1.],A_ub=np.c_[rows[:,1:],-np.ones(len(rows))],b_ub=-rows[:,0],
                       bounds=[(None,None)]*(eng.nvar-1)+[(0,None)],method='highs-ipm',
                       options={'time_limit':120,'threads':1,'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10})
        assert lp.success,lp.message
        weights=-lp.ineqlin.marginals
        active=[i for i,w in enumerate(weights) if w>1e-13]
        todo=[i for i in active if i not in cache]
        print('LP_ACTIVE',v,retry,len(active),'new_exact',len(todo),'margin',lp.fun,flush=True)
        for index in todo:
            full=eng.full_row(all_sources[index]);reduced=eng.reduce(full);scale=max(map(abs,reduced))
            cache[index]=(full,reduced,scale)
            rows[index]=[float(q/scale) for q in reduced] if scale else np.zeros(eng.nvar)
            checks+=1
            if checks%20==0:print('EXACT',v,checks,flush=True)
        if not todo:break
    else:raise RuntimeError('Exact active-row stabilization did not finish')
    terms=[];total=[F(0)]*eng.nvar
    for index in active:
        full,reduced,scale=cache[index]
        if not scale:continue
        source=all_sources[index];w=F(repr(float(weights[index])))
        for j,q in enumerate(reduced):total[j]+=w*q/scale
        terms.append(dict(source=source,weight=str(w),scale=str(scale),full_row={str(i):str(q) for i,q in full.items()}))
    if v==405:bounds=[2**(len(eng.aa[i])*(len(eng.aa[i])-1)//2) for i in eng.cores]
    else:bounds=[1]+[1024 if v==375 else 256]*(eng.nvar-1)
    margin=total[0]+sum((q*b for q,b in zip(total[1:],bounds[1:]) if q<0),F(0))
    report=dict(parameters=[eng.v,eng.k,eng.lam,eng.mu],status='EXACT_POSITIVE_MARGIN' if margin>0 else 'RECOVERY_FAILED',
                terms=terms,combined_coefficients=list(map(str,total)),bounds=bounds,exact_margin=str(margin),decimal_margin=float(margin),
                seconds=time.monotonic()-started,provenance='Fresh cone SDP -> new rounded f -> fresh LP -> exact regenerated rows and bounded residual')
    (out/'certificate.json').write_text(json.dumps(report,indent=2)+'\n')
    print('FINAL',v,report['status'],len(terms),float(margin),report['seconds'],flush=True)
    assert margin>0

if __name__=='__main__':run(int(sys.argv[1]))
