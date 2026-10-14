"""Build article (2.9), (2.6), (4.6) exactly, with LABELLED embedding counts.

No previously saved affine data or certificate directions are read.
Python standard library only. Graph bits use 01,02,...,12,... ordering.
"""
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import comb, isqrt
from pathlib import Path
from functools import lru_cache
import argparse, json, time

HERE = Path(__file__).resolve().parent

def adjacency(n, bits):
    a = [[0]*n for _ in range(n)]
    for (i,j), t in zip(combinations(range(n),2),bits):
        a[i][j] = a[j][i] = int(t)
    return a

def code(a, order):
    return ''.join(str(a[i][j]) for i,j in combinations(order,2))

def catalogue(p):
    graphs, lookup = [], {}
    for n in range(1,p+1):
        for b in product('01', repeat=comb(n,2)):
            bits=''.join(b)
            if (n,bits) in lookup: continue
            a=adjacency(n,bits)
            orbit=[code(a,o) for o in permutations(range(n))]
            assert bits==min(orbit)
            index=len(graphs)
            graphs.append(dict(n=n,bits=bits,aut=orbit.count(bits)))
            for t in orbit: lookup[n,t]=index
    return graphs,lookup

def rref(a, n):
    a=[[Q(x) for x in r] for r in a]
    piv=[]
    for j in range(n):
        k=len(piv)
        i=next((i for i in range(k,len(a)) if a[i][j]),None)
        if i is None: continue
        a[k],a[i]=a[i],a[k]
        c=a[k][j]; a[k]=[x/c for x in a[k]]
        for i in range(len(a)):
            if i!=k and a[i][j]:
                c=a[i][j]; a[i]=[x-c*y for x,y in zip(a[i],a[k])]
        piv.append(j)
    return a,piv

def inverse(g):
    a,p=rref([r+[Q(i==j) for j in range(len(g))] for i,r in enumerate(g)],len(g))
    assert len(p)==len(g)
    return [r[len(g):] for r in a]

def pd_pivots(g):
    a=[r[:] for r in g]; ps=[]
    for i in range(len(g)):
        t=a[i][i]; ps.append(t)
        if t<=0: return None
        for j in range(i+1,len(g)):
            for k in range(i+1,len(g)):
                a[j][k]-=a[j][i]*a[i][k]/t
    return ps

def serial(x):
    if isinstance(x,Q): return str(x)
    if isinstance(x,list): return [serial(t) for t in x]
    if isinstance(x,dict): return {k:serial(v) for k,v in x.items()}
    return x

def build(par, out, p=5, anchor_sizes=None):
    start=time.monotonic()
    v,k,lam,mu=par
    assert k*(k-lam-1)==(v-k-1)*mu
    delta=(lam-mu)**2+4*(k-mu); sq=isqrt(delta)
    assert sq*sq==delta
    r=Q(lam-mu+sq,2); s=Q(lam-mu-sq,2)
    d=Q(-k-(v-1)*s,r-s); assert d.denominator==1
    d=int(d)
    # Projector E_r=(A-sI-(k-s)J/v)/(r-s), then rescale diagonal to 1.
    diag=Q(-s-Q(k-s,v),r-s)
    alpha=(-Q(k-s,v)/(r-s))/diag
    beta=((1-Q(k-s,v))/(r-s))/diag
    anchor_sizes=anchor_sizes or ([3] if p==5 else [3,4]); degrees=[1,2]
    graphs,lookup=catalogue(p); n=len(graphs)
    assert n in (52,208)
    mats=[adjacency(g['n'],g['bits']) for g in graphs]
    A=[]; b=[]; labels=[]
    for hi,g in enumerate(graphs):
        h=g['n']; H=mats[hi]
        if h>=p: continue
        exts=[]
        for prof in product((0,1),repeat=h):
            T=[row+[prof[i]] for i,row in enumerate(H)]+[list(prof)+[0]]
            exts.append((prof,lookup[h+1,code(T,range(h+1))]))
        for B in [()]+[(i,) for i in range(h)]+list(combinations(range(h),2)):
            if not B: e=v-h
            elif len(B)==1: e=k-sum(H[B[0]])
            else:
                i,j=B
                e=(lam if H[i][j] else mu)-sum(H[i][z]*H[j][z] for z in range(h))
            row=[Q(0)]*n; row[hi]=-e
            for prof,fi in exts:
                if all(prof[i] for i in B): row[fi]+=1
            A.append(row); b.append(Q(0)); labels.append(dict(H=g,B=B))
    for t in range(1,p+1):
        A.append([Q(1,g['aut']) if g['n']==t else Q(0) for g in graphs])
        b.append(Q(comb(v,t))); labels.append(dict(total_order=t))
    red,piv=rref([row+[rhs] for row,rhs in zip(A,b)],n)
    assert not any(row[-1] for row in red[len(piv):])
    free=[i for i in range(n) if i not in piv]; q=len(free)
    affine=[[Q(0)]*(q+1) for _ in graphs]
    for j,i in enumerate(free): affine[i][j+1]=Q(1)
    for row,i in zip(red,piv): affine[i]=[row[-1]]+[-row[j] for j in free]
    for row,rhs in zip(A,b):
        assert [sum(row[i]*affine[i][j] for i in range(n) if row[i]) for j in range(q+1)]==[rhs]+[0]*q
    blocks=[]
    for hi,g in enumerate(graphs):
        if g['n'] not in anchor_sizes: continue
        m=g['n']; profs=list(product((0,1),repeat=m)); size=len(profs)
        H=mats[hi]
        G=[[Q(1) if i==j else beta if H[i][j] else alpha for j in range(m)] for i in range(m)]
        pivots=pd_pivots(G)
        if pivots is None: continue
        Gi=inverse(G)
        u=[[beta if bit else alpha for bit in pr] for pr in profs]
        proj=[[sum(u[a][i]*Gi[i][j]*u[c][j] for i in range(m) for j in range(m)) for c in range(size)] for a in range(size)]
        rho=[1-proj[i][i] for i in range(size)]
        def extension(a,c=None,e=None):
            nn=m+1+(c is not None)
            T=[[0]*nn for _ in range(nn)]
            for i in range(m):
                for j in range(m): T[i][j]=H[i][j]
                T[i][m]=T[m][i]=profs[a][i]
                if c is not None: T[i][m+1]=T[m+1][i]=profs[c][i]
            if c is not None: T[m][m+1]=T[m+1][m]=e
            return lookup[nn,code(T,range(nn))]
        singles=[extension(a) for a in range(size)]
        pairs=[[[extension(a,c,e) for e in (0,1)] for c in range(size)] for a in range(size)]
        D=d-m
        for ell in degrees:
            C=[[[Q(0)]*size for _ in range(size)] for _ in range(n)]
            for a in range(size):
                C[singles[a]][a][a]+=rho[a]**ell
                for c in range(size):
                    for e in (0,1):
                        t=(beta if e else alpha)-proj[a][c]
                        val=t if ell==1 else (D*t*t-rho[a]*rho[c])/(D-1)
                        C[pairs[a][c][e]][a][c]+=val
            M=[[[sum(affine[i][j]*C[i][a][c] for i in range(n) if C[i][a][c]) for c in range(size)] for a in range(size)] for j in range(q+1)]
            assert all(cc[a][c]==cc[c][a] for cc in C for a in range(size) for c in range(size))
            blocks.append(dict(m=m,H=g['bits'],ell=ell,G=G,G_pivots=pivots,C=C,M=M))
    data=dict(parameters=par,d=d,alpha=alpha,beta=beta,m=max(anchor_sizes),anchor_sizes=anchor_sizes,p=p,degrees=degrees,
              graphs=graphs,A=A,b=b,affine=affine,free_indices=free,blocks=blocks,
              eigenvalues=[r,s],count_convention='labelled induced embeddings (article 2.1)',
              profile_order=[''.join(map(str,t)) for t in profs],rank=len(piv),
              equation_count=len(A),generation_seconds=time.monotonic()-start)
    out.mkdir(exist_ok=True,parents=True)
    (out/'data.json').write_text(json.dumps(serial(data),indent=2)+'\n')
    print(json.dumps({k:serial(data[k]) for k in ['parameters','d','alpha','beta','rank','equation_count','free_indices','generation_seconds']}),flush=True)
    return data

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--parameters',type=int,nargs=4,default=[540,147,66,30]);ap.add_argument('--out',type=Path,default=HERE)
    args=ap.parse_args();build(args.parameters,args.out)
