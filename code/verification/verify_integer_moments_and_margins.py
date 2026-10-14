"""Independent standard-library audit of integer moments and final arithmetic.

No discovery/recovery modules are imported. Affine coverage and PSD source
validity are certified by the companion independent audits. This audit checks
all saved moment sources, every active row's affine substitution, and the
strictly positive exact bounded-residual margin. Graph lookups are checked
against catalogue representatives by an independent isomorphism search.
"""
from array import array
from collections import Counter
from fractions import Fraction as F
from functools import reduce
from hashlib import sha256
from itertools import combinations
from math import factorial, gcd, lcm
from pathlib import Path
import json, sys, time
if hasattr(sys, 'set_int_max_str_digits'): sys.set_int_max_str_digits(0)
BASE = Path(__file__).resolve().parents[1]
PORTABLE = Path(__file__).resolve().parent
PACKAGED = (PORTABLE/'atlas_graphs.txt').exists()
DATA = PORTABLE if PACKAGED else BASE/'output/paper/support'
CATALOGUE = PORTABLE/'atlas_graphs.txt' if PACKAGED else BASE/'work/atlas/graphs.txt'


def decode(n, mask):
    a = [[0]*n for _ in range(n)]
    for j in range(n):
        for i in range(j): a[i][j] = a[j][i] = (mask >> (j*(j-1)//2+i)) & 1
    return a


def encode(a):
    return sum(a[i][j] << (j*(j-1)//2+i) for j in range(len(a)) for i in range(j))


def iso(a,b):
    if len(a) != len(b): return False
    n=len(a);da=list(map(sum,a));db=list(map(sum,b))
    if sorted(da)!=sorted(db): return False
    candidates=[[j for j in range(n) if da[i]==db[j]] for i in range(n)]
    order=sorted(range(n),key=lambda i:(len(candidates[i]),-da[i]));assigned={};used=set()
    def rec(pos):
        if pos==n:return True
        i=order[pos]
        for j in candidates[i]:
            if j in used or any(a[i][k]!=b[j][l] for k,l in assigned.items()):continue
            assigned[i]=j;used.add(j)
            if rec(pos+1):return True
            del assigned[i];used.remove(j)
        return False
    return rec(0)


def falling(t,j):
    q=1
    for a in range(j):q*=t-a
    return q


def evaluate(c,t):return sum(q*t**j for j,q in enumerate(c))


def primitive(row,eq=False):
    row={int(i):F(q) for i,q in row.items() if q}
    if not row:return {}
    den=lcm(*(q.denominator for q in row.values()));r={i:int(q*den) for i,q in row.items()}
    g=reduce(gcd,map(abs,r.values()))
    if eq and r[min(r)]<0:g=-g
    return {i:q//g for i,q in r.items()}


class Auditor:
    def __init__(self,v):
        self.v=v;self.k,self.lam,self.mu={375:(102,45,21),441:(88,35,13),575:(112,45,16)}[v]
        self.out=DATA/f'rebuild_{v}'
        self.graphs={}
        for line in CATALOGUE.read_text().splitlines():
            i,n,m=map(int,line.split());self.graphs[i]=decode(n,m)
        self.candidates={}
        for hid,a in self.graphs.items():
            self.candidates.setdefault((len(a),tuple(sorted(map(sum,a)))),[]).append(hid)
        self.lookups={}
        self.forbidden={w['id'] for w in json.loads((self.out/f'forbidden_{v}.json').read_text())}

    def identify(self,a):
        key=len(a),encode(a)
        if key not in self.lookups:
            matches=[hid for hid in self.candidates[(len(a),tuple(sorted(map(sum,a))))] if iso(a,self.graphs[hid])]
            assert len(matches)==1,('not a unique graph representative',key,matches)
            self.lookups[key]=matches[0]
        return self.lookups[key]

    def moments(self):
        data=json.loads((self.out/'integer_moments.json').read_text());byanchor={};domains={}
        for record in data['anchors']:
            hid=record['anchor'];a=self.graphs[hid];assert len(a)==3 and encode(a)==record['anchor_mask']
            assert record['override'] is None,'External domain overrides require additional audit'
            marginals={0:(self.v-3,0),7:(0,1)}
            for i in range(3):marginals[1<<i]=(self.k-sum(a[i]),0)
            for i,j in combinations(range(3),2):
                marginals[(1<<i)|(1<<j)]=((self.lam if a[i][j] else self.mu)-sum(a[i][s]*a[j][s] for s in range(3)),0)
            cells=[]
            for bits in range(8):
                c=s=0
                for sup,(mc,ms) in marginals.items():
                    if bits&sup==bits:
                        sign=(-1)**(bin(sup).count('1')-bin(bits).count('1'));c+=sign*mc;s+=sign*ms
                p=[(bits>>i)&1 for i in range(3)]
                b=[row+[p[i]] for i,row in enumerate(a)]+[p+[0]]
                target=self.identify(b)
                cell={'profile':p,'constant':c,'slope':s,'extension':target,'forced_zero':target in self.forbidden}
                cells.append(cell)
            assert cells==record['cells']
            upper=min(marginals[q][0] for q in (0,3,5,6))
            domain=[t for t in range(upper+1) if all(c['constant']+c['slope']*t>=0 and (not c['forced_zero'] or c['constant']+c['slope']*t==0) for c in cells)]
            assert domain==record['domain']==record['base_domain'];domains[hid]=domain
            maps=[{hid:F(1)}];recipes=[{'order':0,'factor':'1','targets':[hid]}]
            for j in range(1,5):
                q=F(falling(self.v-3,j),2**(3*j+j*(j-1)//2));row={};targets=[]
                edges=list(combinations(range(3,3+j),2))
                for bits in range(1<<len(edges)):
                    b=[[0]*(3+j) for _ in range(3+j)]
                    for x in range(3):
                        for y in range(3):b[x][y]=a[x][y]
                        for y in range(3,3+j):b[x][y]=b[y][x]=1
                    for bit,(x,y) in enumerate(edges):b[x][y]=b[y][x]=(bits>>bit)&1
                    target=self.identify(b);targets.append(target);row[target]=row.get(target,F(0))+q
                maps.append(row);recipes.append({'order':j,'factor':str(q),'targets':targets})
            assert maps==[{int(i):F(q) for i,q in r.items()} for r in record['moment_rows']]
            assert recipes==record['moment_recipes'];byanchor[hid]=maps
        checked={'eq':set(),'ineq':set()}
        for w in data['witnesses']:
            kind=w['kind'];idx=w['row_index'];hid=w['anchor'];domain=domains[hid]
            p=list(map(F,w['polynomial']));assert len(p)<=5
            values=[evaluate(p,t) for t in domain]
            assert all(x==0 for x in values) if kind=='eq' else all(x>=0 for x in values)
            # Independently compute the Newton coefficient using its binomial formula.
            c=[sum((-1)**(j-t)*F(factorial(j),factorial(t)*factorial(j-t))*evaluate(p,t) for t in range(j+1))/factorial(j) for j in range(len(p))]
            assert c==list(map(F,w['falling_coefficients']))
            assert all(evaluate(p,t)==sum(q*falling(t,j) for j,q in enumerate(c)) for t in range(-5,10))
            row={}
            for q,m in zip(c,byanchor[hid]):
                for i,x in m.items():row[i]=row.get(i,F(0))+(1 if kind=='eq' else -1)*q*x
            row=primitive(row,kind=='eq')
            assert row=={int(i):int(q) for i,q in data[kind][idx].items()}
            checked[kind].add(idx)
        assert all(checked[kind]==set(range(len(data[kind]))) for kind in checked)
        return data,{'anchors':len(domains),'domains':domains,'equalities':len(checked['eq']),'inequalities':len(checked['ineq']),'witnesses':len(data['witnesses'])}

    def final(self,moments):
        cert=json.loads((self.out/'recovered_exact_certificate.json').read_text())
        rec={int(i):{int(j):F(q) for j,q in r.items()} for i,r in json.loads((self.out/f'affine_recipes_{self.v}.json').read_text()).items()}
        dim=cert['dimension'];B=cert['variable_upper_bound'];assert B==(1024 if self.v==375 else 256)
        elim=json.loads((self.out/f'exact_elimination_{self.v}.json').read_text())
        ids=elim['free_core_ids'];assert len(ids)==dim-1
        bounds=[]
        for coord,hid in enumerate(ids,1):
            a=self.graphs[hid];n=len(a);assert min(map(sum,a))>=3 and 4<=n<=7
            reached={0}
            while True:
                new=reached|{j for i in reached for j in range(n) if a[i][j]}
                if new==reached:break
                reached=new
            assert len(reached)==n
            upper=F(2**(n*(n-1)//2)*self.k**(n-1),falling(self.v-1,n-1))
            assert upper<B and rec[hid]=={coord:F(1)};bounds.append(upper)
        total=[F(0)]*dim;types=Counter();active_quartics=[]
        for item in cert['rows']:
            source=item['source'];kind=source['kind'];types[kind]+=1
            weight=F(item['weight']);assert weight>=0
            full=item['full_y_row'];r=[F(0)]*dim
            if kind=='upper_bound':
                assert full is None and source['bound']==B and 0<source['coordinate']<dim
                r[0]=-B;r[source['coordinate']]=1
            else:
                full={int(i):F(q) for i,q in full.items()}
                if kind=='nonnegative_count':assert full=={source['graph_id']:F(-1)}
                elif kind=='integer_moment':
                    assert full=={int(i):F(q) for i,q in moments['ineq'][source['row_index']].items()}
                    ws=[w for w in moments['witnesses'] if w['kind']=='ineq' and w['row_index']==source['row_index']]
                    assert ws
                    w=ws[0];nodes=w['nodes'];assert len(nodes)==4 and nodes[1]==nodes[0]+1 and nodes[3]==nodes[2]+1 and w['sign']==1
                    p=[F(1)]
                    for t in nodes:
                        nxt=[F(0)]*(len(p)+1)
                        for j,q in enumerate(p):nxt[j]-=t*q;nxt[j+1]+=q
                        p=nxt
                    assert p==list(map(F,w['polynomial']))
                    active_quartics.append({'anchor':w['anchor'],'nodes':nodes,'globally_nonnegative_on_integers':True})
                else:assert kind=='PSD' # independently source-audited by the companion tool.
                for i,q in full.items():
                    for j,p in rec[i].items():r[j]+=q*p
            norm=max(map(abs,r));assert norm>0 and norm==F(item['normalization'])
            r=[q/norm for q in r];assert r==list(map(F,item['reduced_normalized_row']))
            for j,q in enumerate(r):total[j]+=weight*q
        assert total==list(map(F,cert['combined_coefficients']))
        margin=total[0]+B*sum(x for x in total[1:] if x<0)
        assert margin==F(cert['exact_margin']) and margin>0
        assert len(cert['rows'])==cert['active_row_count']
        return {'dimension':dim,'free_core_counts':len(ids),'bound':B,'maximum_spanning_tree_bound':str(max(bounds)),
                'active_types':dict(types),'active_rows':len(cert['rows']),'active_quartics':active_quartics,'exact_margin':str(margin),'margin_decimal':float(margin),
                'certificate_sha256':sha256((self.out/'recovered_exact_certificate.json').read_bytes()).hexdigest()}


def main(v):
    start=time.time();a=Auditor(v);moments,report=a.moments();final=a.final(moments)
    result={'v':v,'status':'PASS','scope':'All integer-moment sources; all active affine substitutions and final rational arithmetic. Requires companion affine-coverage and PSD source audits.',
            'moments':report,'final':final,'independently_checked_graph_lookups':len(a.lookups),'seconds':time.time()-start}
    target=DATA/f'moments_certificate_audit_{v}.json';target.write_text(json.dumps(result,indent=2))
    print(v,'PASS',final['active_rows'],'rows','margin',final['margin_decimal'],'seconds',round(result['seconds'],3),flush=True)
    return result
if __name__=='__main__':
    for v in map(int,sys.argv[1:] or [375,441,575]):main(v)
