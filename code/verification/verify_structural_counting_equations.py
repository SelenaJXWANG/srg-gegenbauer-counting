#!/usr/bin/env python3
"""Standard-library audit of reconstructed last-three affine coordinates.

Regenerates graph isomorphism classes, checks triangular deletion recipes,
every short-vector/completion witness, all eliminated equations, the free
coordinate identity, and a modular lower bound on the exact equation rank.
Together these prove that the affine parametrization covers every solution
of the necessary equations; checking substitution alone would not suffice.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,product
from math import gcd,lcm,prod
from functools import reduce
import importlib.util,json,sys,time

ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('raw_graph_orbits',ROOT/'six_exact_data/verify_raw_v2_certificate.py')
raw=importlib.util.module_from_spec(spec);spec.loader.exec_module(raw)


def adjacency(n,mask):
    a=[[0]*n for _ in range(n)]
    for j in range(n):
        for i in range(j):a[i][j]=a[j][i]=(mask>>(j*(j-1)//2+i))&1
    return a


def encode(a):return sum(a[i][j]<<(j*(j-1)//2+i) for j in range(len(a)) for i in range(j))


def old_code(a):
    code=0
    for i,j in combinations(range(len(a)),2):code=2*code+a[i][j]
    return code


def graph_catalog():
    graphs={};canonical_to_id={}
    for line in (ROOT/'atlas_graphs.txt').read_text().splitlines():
        hid,n,mask=map(int,line.split());a=adjacency(n,mask)
        key=n,raw.canonical(n,old_code(a))
        assert key not in canonical_to_id
        graphs[hid]=a;canonical_to_id[key]=hid
    assert sorted(graphs)==list(range(1252))
    for n,count in enumerate((1,2,4,11,34,156,1044),1):
        assert sum(len(a)==n for a in graphs.values())==count
        assert {key[1] for key in canonical_to_id if key[0]==n}==set(raw.universe(n)[1])
    def identify(a):return canonical_to_id[len(a),raw.canonical(len(a),old_code(a))]
    return graphs,identify


def extend(a,alpha):return [row+[q] for row,q in zip(a,alpha)]+[list(alpha)+[0]]


def combine(rows,coefficients):
    out={}
    for row,c in zip(rows,coefficients):
        for j,q in row.items():out[j]=out.get(j,F(0))+c*q
    return {j:q for j,q in out.items() if q}


def primitive(row):
    row={int(i):F(q) for i,q in row.items() if q}
    den=lcm(*(q.denominator for q in row.values()))
    ints={i:int(q*den) for i,q in row.items()};g=gcd(*ints.values())
    if ints[min(ints)]<0:g=-g
    return {i:q//g for i,q in sorted(ints.items())}


def modrank(matrix,p):
    a=[r[:] for r in matrix];rank=0
    for j in range(len(a[0])):
        pivot=next((i for i in range(rank,len(a)) if a[i][j]),None)
        if pivot is None:continue
        a[rank],a[pivot]=a[pivot],a[rank];inv=pow(a[rank][j],-1,p)
        a[rank]=[(q*inv)%p for q in a[rank]]
        for i in range(rank+1,len(a)):
            c=a[i][j]
            if c:a[i]=[(q-c*r)%p for q,r in zip(a[i],a[rank])]
        rank+=1
    return rank


def audit(v,graphs,identify):
    started=time.time();folder=ROOT/f'rebuild_{v}'
    if not folder.exists():folder=ROOT.parent/f'support/rebuild_{v}'
    k,lam,mu,r,f={375:(102,45,21,27,34),441:(88,35,13,25,44),575:(112,45,16,32,46)}[v]
    core=json.loads((folder/f'core_recipes_{v}.json').read_text())
    ids=core['cores'];R={int(i):{int(j):F(q) for j,q in row.items()} for i,row in core['expressions'].items()}
    assert ids==[i for i,a in graphs.items() if len(a)==1 or min(map(sum,a))>=3]
    assert len(ids)==174
    for j,hid in enumerate(ids):assert R[hid]=={j:F(1)}
    final={int(i):{int(j):F(q) for j,q in row.items()} for i,row in json.loads((folder/f'affine_recipes_{v}.json').read_text()).items()}
    elimination=json.loads((folder/f'exact_elimination_{v}.json').read_text())
    T=[[F(q) for q in row] for row in elimination['T']];nfree=len(T[0])
    assert len(T)==174 and nfree==111
    for j,hid in enumerate(ids):assert final[hid]=={i:q for i,q in enumerate(T[j]) if q}
    free=[0]+elimination['free_core_ids'];assert len(free)==nfree and len(set(free))==nfree
    for j,hid in enumerate(free):assert final[hid]=={j:F(1)}
    checked=0;seen=set()
    for hid in sorted(graphs,key=lambda i:(len(graphs[i]),-sum(map(sum,graphs[i])))):
        a=graphs[hid];n=len(a)
        if hid in ids:seen.add(hid);continue
        x=min(range(n),key=lambda i:sum(a[i]));keep=[i for i in range(n) if i!=x]
        h=[[a[i][j] for j in keep] for i in keep];parent=identify(h)
        S=[i for i,j in enumerate(keep) if a[x][j]];assert len(S)<=2
        if len(S)==0:c=v-n+1
        elif len(S)==1:c=k-sum(h[S[0]])
        else:
            i,j=S;c=(lam if h[i][j] else mu)-sum(h[i][z]*h[j][z] for z in range(n-1))
        row={parent:F(2**(n-1)*c,v-n+1)}
        rest=[i for i in range(n-1) if i not in S]
        for values in product((0,1),repeat=len(rest)):
            if not any(values):continue
            alpha=[int(i in S) for i in range(n-1)]
            for i,b in zip(rest,values):alpha[i]=b
            target=identify(extend(h,alpha));row[target]=row.get(target,F(0))-1
        assert all(i in seen for i in row)
        assert combine([R[i] for i in row],row.values())==R[hid]
        assert combine([final[i] for i in row],row.values())==final[hid]
        checked+=1;seen.add(hid)
    forbidden=json.loads((folder/f'forbidden_{v}.json').read_text());equations=[]
    for record in forbidden:
        hid=record['id'];a=graphs[hid];h=record['root'];n=len(a)
        assert record['mask']==encode(a) and h is not None and len(h)==n
        z=[3*h[i]+sum(a[i][j]*h[j] for j in range(n)) for i in range(n)]
        assert sum(x*y for x,y in zip(h,z)) in (1,2)
        equations.append({hid:F(1)})
    complete=json.loads((folder/'completion_constraints.json').read_text());eq=[{int(i):F(q) for i,q in row.items()} for row in complete['eq']]
    used=set()
    for w in complete['witnesses']:
        hid=w['anchor'];a=graphs[hid];n=len(a);h=w['h'];p=w['profile']
        assert n==w['n'] and encode(a)==w['anchor_mask'] and sum(h)==1
        z=[3*h[i]+sum(a[i][j]*h[j] for j in range(n)) for i in range(n)]
        assert z==p and all(x in (0,1) for x in z) and sum(x*y for x,y in zip(h,z))==3
        full=extend(a,p)
        if w['type']=='unique_profile_completion':
            target=identify(full);assert target==w['completed_graph'] and n<=6
            row={hid:F(-1),target:F(v-n,2**n)}
        else:
            assert w['type']=='seven_card_bijection' and n==7 and encode(full)==w['completed_8_mask']
            null=h+[-1];assert null==w['null'] and sum(null)==0
            assert all(3*null[i]+sum(full[i][j]*null[j] for j in range(8))==0 for i in range(8))
            deleted=w['deleted_vertex'];keep=[i for i in range(8) if i!=deleted]
            assert keep==w['card_vertices'] and null[deleted] in (-1,1)
            card=[[full[i][j] for j in keep] for i in keep];target=identify(card)
            assert target==w['card']
            hc=[-null[i]//null[deleted] for i in keep];assert hc==w['card_h'] and sum(hc)==1
            zc=[3*hc[i]+sum(card[i][j]*hc[j] for j in range(7)) for i in range(7)]
            assert zc==[full[i][deleted] for i in keep] and sum(x*y for x,y in zip(hc,zc))==3
            row={hid:F(1),target:F(-1)}
        assert primitive(row)==primitive(eq[w['row_index']]);used.add(w['row_index'])
    assert used==set(range(len(eq)))
    equations+=eq
    for row in equations:assert not combine([final[i] for i in row],row.values())
    # Every necessary equation is reconstructed above; no unproved extra
    # equalities can be silently introduced by the saved affine recipes.
    p=1000000007;E=[]
    for equation in equations:
        row=combine([R[i] for i in equation],equation.values())
        assert all(q.denominator%p for q in row.values())
        E.append([(row.get(j,F(0)).numerator*pow(row.get(j,F(0)).denominator,-1,p))%p for j in range(174)])
    rank=modrank(E,p);assert rank==174-nfree==63
    # T has 111 identity rows and ET=0, so rank(E)<=63 over Q. The modular
    # lower bound is therefore sharp, proving image(T)=ker(E) over Q.
    largest=max(F(2**(n*(n-1)//2)*k**(n-1),prod(range(v-n+1,v))) for n in range(4,8))
    B=1024 if v==375 else 256;assert largest<B
    for hid in free[1:]:assert min(map(sum,graphs[hid]))>=3
    result={'status':'PASS','v':v,'graphs':1252,'core_coordinates':174,'free_coordinates':111,
            'triangular_deletion_rows':checked,'forbidden_witnesses':len(forbidden),
            'completion_rows':len(eq),'completion_witnesses':len(complete['witnesses']),
            'modular_rank_prime':p,'modular_rank':rank,'exact_rank':rank,
            'spanning_tree_bound':str(largest),'safe_free_coordinate_bound':B,
            'seconds':time.time()-started,'meaning':'All source equations valid under the independently proved complete-shortest-shell theorem. The saved affine map parametrizes their entire solution space, not just a subset.'}
    (folder/'affine_source_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
    return result


if __name__=='__main__':
    started=time.time();graphs,identify=graph_catalog();print('Exact graph atlas reconstructed',round(time.time()-started,3),flush=True)
    for v in map(int,sys.argv[1:] or ('375','441','575')):audit(v,graphs,identify)
