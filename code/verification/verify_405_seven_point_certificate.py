"""Standard-library exact replay of the PURE seven-point 405 certificate.

Independently expands every selected degree 0,1,2 profile-kernel energy.
Atlas canonicalization and every rational elimination recipe are checked
separately by the companion pure-seven-point audit.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
from array import array
from collections import Counter
import json
import sys

if hasattr(sys,'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)
HERE = Path(__file__).resolve().parent
ROOT = HERE
OUT = HERE/'data/degree2'


def matrix(n,mask):
    return [[int(i!=j and bool(mask & (1 << (max(i,j)*(max(i,j)-1)//2+min(i,j)))))
             for j in range(n)] for i in range(n)]


def invert_positive(Q):
    n = len(Q)
    gram = [row[:] for row in Q]
    for _ in range(n):
        pivot = gram[0][0]
        assert pivot>0
        gram = [[gram[i][j]-gram[i][0]*gram[0][j]/pivot
                 for j in range(1,len(gram))] for i in range(1,len(gram))]
    rows = [list(Q[i])+[F(i==j) for j in range(n)] for i in range(n)]
    for i in range(n):
        pivot = rows[i][i]
        assert pivot>0
        rows[i] = [q/pivot for q in rows[i]]
        for j in range(n):
            if j==i:
                continue
            factor = rows[j][i]
            rows[j] = [a-factor*b for a,b in zip(rows[j],rows[i])]
    return [row[n:] for row in rows]


def main():
    graphs = []
    for line in (ROOT/'atlas/graphs.txt').read_text().splitlines():
        i,n,mask = map(int,line.split())
        assert i==len(graphs)
        graphs.append((n,mask))
    maps = {}
    for n in range(1,8):
        data = array('h')
        data.frombytes((ROOT/f'atlas/map{n}.bin').read_bytes())
        if sys.byteorder != 'little':
            data.byteswap()
        assert len(data)==2**(n*(n-1)//2)
        maps[n]=data
    recipes = json.loads((OUT/'core_recipes_405.json').read_text())
    coreids = recipes['cores']
    recipes = {int(i):{int(j):F(q) for j,q in row.items()}
               for i,row in recipes['expressions'].items()}
    record = json.loads((OUT/'pure_exact_certificate.json').read_text())
    assert record['parameters']==[405,132,63,33]
    geometry = {}
    def extension(hid,p,q=None,e=0):
        n,mask=graphs[hid]
        value = mask+sum(int(b)<<(n*(n-1)//2+i) for i,b in enumerate(p))
        if q is None:
            return maps[n+1][value]
        value += sum(int(b)<<(n*(n+1)//2+i) for i,b in enumerate(q))
        value += e<<(n*(n+1)//2+n)
        return maps[n+2][value]
    def geom(hid):
        if hid in geometry:
            return geometry[hid]
        n,mask=graphs[hid]
        assert n<=5
        A=matrix(n,mask)
        Q=[[F(9*(i==j)+3*A[i][j]-1,8) for j in range(n)] for i in range(n)]
        inv=invert_positive(Q)
        profiles=list(product((0,1),repeat=n))
        u=[[F(3*b-1,8) for b in p] for p in profiles]
        projection=[[sum((u[i][a]*inv[a][b]*u[j][b] for a in range(n) for b in range(n)),F(0))
                     for j in range(len(profiles))] for i in range(len(profiles))]
        heights=[1-projection[i][i] for i in range(len(profiles))]
        singles=[extension(hid,p) for p in profiles]
        pairs=[[[extension(hid,p,q,e) for e in (0,1)] for q in profiles] for p in profiles]
        geometry[hid]=(n,projection,heights,singles,pairs)
        return geometry[hid]
    total=[F(0)]*len(coreids)
    inventory=Counter()
    for index,term in enumerate(record['terms']):
        source=term['source']
        row={}
        if source['kind']=='nonnegative_count':
            row={source['graph_id']:F(-1)}
            inventory['nonnegative counts']+=1
        else:
            assert source['kind']=='PSD'
            hid,ell,w=source['anchor'],source['degree'],source['direction']
            assert ell in (0,1,2) and all(isinstance(v,int) for v in w)
            n,projection,h,singles,pairs=geom(hid)
            assert len(w)==len(h)
            D=30-n
            for i,wi in enumerate(w):
                row[singles[i]]=row.get(singles[i],F(0))-wi*wi*h[i]**ell
                for j,wj in enumerate(w):
                    for e,target in enumerate(pairs[i][j]):
                        z=F(3*e-1,8)-projection[i][j]
                        kernel=(F(1),z,(D*z*z-h[i]*h[j])/(D-1))[ell]
                        coefficient=-F(404-n,2**(n+1))*wi*wj*kernel
                        row[target]=row.get(target,F(0))+coefficient
            row={i:q for i,q in row.items() if q}
            inventory[f'degree {ell} PSD energies']+=1
        assert row=={int(i):F(q) for i,q in term['full_row'].items()},('source mismatch',index)
        reduced=[F(0)]*len(coreids)
        for i,q in row.items():
            for j,r in recipes[i].items():
                reduced[j]+=q*r
        scale=max(abs(q) for q in reduced)
        assert scale==F(term['scale']) and scale>0
        weight=F(term['weight'])
        assert weight>0
        for j,q in enumerate(reduced):
            total[j]+=weight*q/scale
    assert total==list(map(F,record['combined_core_coefficients']))
    bounds=[2**(graphs[i][0]*(graphs[i][0]-1)//2) for i in coreids]
    assert bounds==record['core_density_upper_bounds']
    margin=total[0]+sum((q*b for q,b in zip(total[1:],bounds[1:]) if q<0),F(0))
    assert margin==F(record['exact_margin'])
    assert margin>F(1,20000000)
    report={'status':'PASS','parameters':[405,132,63,33],
            'terms_checked':len(record['terms']),'inventory':dict(inventory),
            'exact_margin':str(margin),'certified_simple_lower_bound':'1/20000000',
            'decimal_margin':float(margin),'pre_forbidden_counts':0,
            'maximum_anchor_size':max(graphs[i][0] for i in geometry),
            'maximum_harmonic_degree':2,'maximum_count_vertices':7,
            'scope':'every energy expansion, positive coefficient, exact reduction and final residual payment replayed; combine with independent atlas/extension audit'}
    (OUT/'pure_exact_replay.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: pure seven-point 405 certificate; exact margin > 1/20000000')
    print(dict(inventory))


if __name__=='__main__':
    main()
