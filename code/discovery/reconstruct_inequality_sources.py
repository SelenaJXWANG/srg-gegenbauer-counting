"""Portable exact source construction for the four seven-point models.

Reads only graph atlas, affine count formulas and proved structural conditions.
Does not read any archived certificate, f direction, cut or multiplier.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]

def inverse(g):
    n=len(g);a=[row[:]+[F(i==j) for j in range(n)] for i,row in enumerate(g)]
    for i in range(n):
        assert a[i][i]>0,'reference Gram not positive definite'
        pivot=a[i][i];a[i]=[q/pivot for q in a[i]]
        for j in range(n):
            if j!=i:
                t=a[j][i];a[j]=[q-t*r for q,r in zip(a[j],a[i])]
    return [row[n:] for row in a]

class ExactModel:
    def __init__(self,v):
        self.v=v
        self.k,self.lam,self.mu,self.f={375:(102,45,21,34),405:(132,63,33,30),441:(88,35,13,44),575:(112,45,16,46)}[v]
        self.eta=F(self.k+3,v)
        atlas=ROOT/'article_replay/pure405/atlas'
        self.aa=[]
        for line in (atlas/'graphs.txt').read_text().splitlines():
            i,n,mask=map(int,line.split());assert i==len(self.aa)
            a=np.zeros((n,n),dtype=np.int64)
            for j in range(n):
                for k in range(j):a[k,j]=a[j,k]=(mask>>(j*(j-1)//2+k))&1
            self.aa.append(a)
        self.maps={n:np.fromfile(atlas/f'map{n}.bin',dtype='<i2') for n in range(1,8)}
        if v==405:
            raw=json.loads((ROOT/'article_replay/pure405/data/degree2/core_recipes_405.json').read_text())
            recipes=raw['expressions'];self.cores=raw['cores'];self.forbidden=set();self.moments={'ineq':[]}
        else:
            folder=ROOT/f'article_replay/legacy_ten/rebuild_{v}'
            recipes=json.loads((folder/f'affine_recipes_{v}.json').read_text())
            self.forbidden={r['id'] for r in json.loads((folder/f'forbidden_{v}.json').read_text())}
            self.moments=json.loads((folder/'integer_moments.json').read_text())
        self.recipes={int(i):{int(j):F(q) for j,q in row.items()} for i,row in recipes.items()}
        self.nvar=1+max(j for row in self.recipes.values() for j in row)
        self.cache={}

    def identify(self,A):
        mask=sum(int(A[i,j])<<(j*(j-1)//2+i) for j in range(len(A)) for i in range(j))
        return int(self.maps[len(A)][mask])

    def extended(self,A,p,q=None,edge=0):
        n=len(A);size=n+1+(q is not None);B=np.zeros((size,size),dtype=np.int64);B[:n,:n]=A
        B[n,:n]=B[:n,n]=p
        if q is not None:B[n+1,:n]=B[:n,n+1]=q;B[n,n+1]=B[n+1,n]=edge
        return self.identify(B)

    def geometry(self,hid):
        if hid in self.cache:return self.cache[hid]
        A=self.aa[hid];m=len(A);c=3-self.eta
        G=[[F(int(A[i,j])+3*(i==j))-self.eta for j in range(m)] for i in range(m)]
        if self.v==405:G=[[q/c for q in row] for row in G]
        inv=inverse(G);pats=[];u=[];h=[];singles=[]
        for p in product((0,1),repeat=m):
            target=self.extended(A,p)
            if target in self.forbidden:continue
            vector=[(F(t)-self.eta)/(c if self.v==405 else 1) for t in p]
            ht=(1 if self.v==405 else c)-sum(vector[i]*inv[i][j]*vector[j] for i in range(m) for j in range(m))
            if self.v!=405 and ht<0:continue
            pats.append(p);u.append(vector);h.append(ht);singles.append(target)
        n=len(pats)
        proj=[[sum(u[i][a]*inv[a][b]*u[j][b] for a in range(m) for b in range(m)) for j in range(n)] for i in range(n)]
        pairs=[[[self.extended(A,p,q,e) for e in (0,1)] for q in pats] for p in pats]
        result=(m,pats,h,singles,proj,pairs);self.cache[hid]=result;return result

    def full_row(self,source):
        kind=source['kind']
        if kind=='nonnegative_count':return {source['graph_id']:F(-1)}
        if kind=='integer_moment':return {int(i):F(q) for i,q in self.moments['ineq'][source['row_index']].items()}
        assert kind=='PSD'
        m,pats,h,ones,proj,pairs=self.geometry(source['anchor']);ell=source['degree'];w=source['direction']
        assert len(w)==len(pats);D=self.f-m;c=3-self.eta;row={};factor=F(self.v-m-1,2**(m+1))
        def add(i,q):
            if q and i not in self.forbidden:row[i]=row.get(i,F(0))-q
        for i,wi in enumerate(w):
            diag=h[i]**ell*(1-F(1,D) if ell==2 and self.v!=405 else 1)
            add(ones[i],wi*wi*diag)
            for j,wj in enumerate(w):
                for e,target in enumerate(pairs[i][j]):
                    t=(F(e)-self.eta)/(c if self.v==405 else 1)-proj[i][j]
                    val=[F(1),t,(D*t*t-h[i]*h[j])/(D-1) if self.v==405 else t*t-h[i]*h[j]/D][ell]
                    add(target,wi*wj*factor*val)
        return {i:q for i,q in row.items() if q}

    def reduce(self,full):
        row=[F(0)]*self.nvar
        for i,q in full.items():
            for j,r in self.recipes[i].items():row[j]+=q*r
        return row
