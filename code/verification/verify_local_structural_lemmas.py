"""Standard-library exact replay of local ingredients, NOT final certificates."""
from fractions import Fraction as F
from itertools import product,permutations,combinations
from math import comb,prod,lcm
from pathlib import Path
import json
OUT=Path(__file__).resolve().parent

def inverse(m):
    n=len(m);a=[[F(x) for x in row]+[F(i==j) for j in range(n)] for i,row in enumerate(m)]
    for j in range(n):
        k=next(i for i in range(j,n) if a[i][j]);a[j],a[k]=a[k],a[j]
        p=a[j][j];a[j]=[x/p for x in a[j]]
        for i in range(n):
            if i!=j:
                p=a[i][j];a[i]=[x-p*y for x,y in zip(a[i],a[j])]
    return [row[n:] for row in a]

def quadratic(m,h):return sum(h[i]*m[i][j]*h[j] for i in range(len(h)) for j in range(len(h)))

cycle=[]
for d in (2,3,4):
    n=2*d;v=[1]*d+[-1]*d;D=[[2*int(i==j) for j in range(n)] for i in range(n)]
    for i in range(d):
        for j in {i,(i+1)%d}:D[i][d+j]=D[d+j][i]=1
    inv=inverse([[D[i][j]+v[i]*v[j] for j in range(n)] for i in range(n)])
    energies=[]
    for p in product((0,1),repeat=n):
        if sum(p[:d])!=sum(p[d:]):continue
        energies.append(quadratic(inv,p))
    denom=lcm(*(x.denominator for x in energies))
    assert denom=={2:4,3:6,4:8}[d]
    cycle.append({'length':n,'balanced_profiles':len(energies),'denominator_lcm':denom,
                  'energies':sorted(set(map(str,energies)))})

eta=F(13,63);c=28;beta=F(1,2)-eta
for T in range(29):
    energy=F(7,3) # Any rational value suffices for the polynomial identity.
    projected=energy-eta*T+eta*eta*c/2-beta*(T-c*eta)**2/(4*(1+beta*c/2))
    gamma=3-eta-projected
    assert gamma==F(37*T*T+104*T+6896,2576)-energy
    assert (37*T*T+104*T+6896)%7==(2*(T-2)**2)%7
    assert (37*T*T+104*T+6896)%23==(14*(T-16)**2)%23
possible=[T for T in range(29) if (T-2)%7==0 and (T-16)%23==0]
assert possible==[16] and F(37*16**2+104*16+6896,2576)==7 and F(16,2)>7
types=[]
for n2 in range(15):
    for n1 in range(15-n2):
        T=n1+2*n2;E=F(3,4)*n1+n2;gamma=F(37*T*T+104*T+6896,2576)-E
        if gamma<0:continue
        gap=gamma-F(11,161)-F(237,8372)*n2
        assert gap>=0
        types.append({'n1':n1,'n2':n2,'gamma':str(gamma),'gap':str(gap)})
assert len(types)==14
total=385*F(11,161)+98*F(237,8372)
assert total==F(17389,598)>28

parameters={441:(88,35,13,25,44),575:(112,45,16,32,46),375:(102,45,21,27,34)}
core_bounds={}
for v,(k,lam,mu,r,f) in parameters.items():
    assert k==3*r+mu and lam-mu==r-3
    assert (3*(v-1)-k)//(r+3)==f
    bound=max(F(2**comb(h,2)*k**(h-1),prod(range(v-h+1,v))) for h in range(4,8))
    safe=1024 if v==375 else 256
    assert bound<safe
    core_bounds[str(v)]={'exact_spanning_tree_bound':str(bound),'safe_bound':safe}
v,k,lam=375,102,45
NI3=comb(v,3)-v*k*(v-2)//2+v*comb(k,2)-v*k*lam//6
mass=F(8*6*NI3,v*(v-1)*(v-2))
assert mass==F(12160,4103)
for t in range(22):assert t+3*(21-t)+3*(60+t)+(129-t)==372

# The two seven-cards share one explicit eight-vertex completion.
a=[[0]*8 for _ in range(8)]
edges=list(combinations(range(4),2))+[(4,6),(5,7)]+[(i,j) for i in (4,5) for j in (0,1)]+[(i,j) for i in (6,7) for j in (2,3)]
for i,j in edges:a[i][j]=a[j][i]=1
cards=[]
for deleted,h in [(3,[1,1,-1,-1,-1,1,1]),(7,[-1,-1,1,1,1,1,-1])]:
    keep=[i for i in range(8) if i!=deleted]
    m=[[a[i][j]+3*int(i==j) for j in keep] for i in keep]
    assert sum(h)==1 and quadratic(m,h)==3
    z=[sum(m[i][j]*h[j] for j in range(7)) for i in range(7)]
    assert z==[a[i][deleted] for i in keep]
    adj=[[a[i][j] for j in keep] for i in keep]
    aut=sum(all(adj[i][j]==adj[p[i]][p[j]] for i in range(7) for j in range(i)) for p in permutations(range(7)))
    assert aut==4
    cards.append({'deleted_vertex':deleted,'h':h,'profile':z,'automorphisms':aut})
report={'status':'PASS local ingredients only; final Farkas arithmetic is outside this verifier',
        'cycle_energies':cycle,'441_zero_height_T_candidates':possible,
        '441_nonnegative_profiles':types,'441_total_height_lower_bound':str(total),
        'core_upper_bounds':core_bounds,'375_I3_mass':str(mass),'seven_card_completion':cards}
(OUT/'last_three_local_replay.json').write_text(json.dumps(report,indent=2))
print('PASS: cycle denominators, 441 modular and 120-profile checks, all core bounds, 375 I3 mass, two seven-card completions.')
print('Final global nonexistence certificates were NOT checked by this script.')
