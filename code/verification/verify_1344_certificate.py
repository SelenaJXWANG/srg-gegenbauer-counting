#!/usr/bin/env python3
"""Reconstruct the five-point 1344 certificate (and the alternate 540 one).

Python standard library only. The five directions and final rational values
are transcribed from srg_540_1344_proofs_zh.pdf, pp.7-8. All graph orbits,
count equations, elimination coefficients and kernel coefficients are rebuilt.
No original discovery or verification code is imported. This is a new
reconstruction, not a recovered copy of the missing original 1344 program.
"""
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from pathlib import Path
from math import factorial
import json


def matrix(n,bits):
    a=[[0]*n for _ in range(n)]
    for (i,j),bit in zip(combinations(range(n),2),bits):a[i][j]=a[j][i]=int(bit)
    return a


def code(a,order=None):
    if order is None:order=range(len(a))
    return ''.join(str(a[i][j]) for i,j in combinations(order,2))


def atlas():
    graphs=[];lookup={};autos={}
    for n,expected in enumerate((1,2,4,11,34),1):
        count=0
        for number in range(1<<(n*(n-1)//2)):
            bits=format(number,f'0{n*(n-1)//2}b') if n>1 else ''
            if (n,bits) in lookup:continue
            a=matrix(n,bits);orbit=[code(a,p) for p in permutations(range(n))]
            representative=min(orbit);assert representative==bits
            key=(n,bits);graphs.append(key);count+=1
            for x in orbit:lookup[n,x]=key
            autos[key]=orbit.count(bits)
            assert len(set(orbit))*autos[key]==factorial(n)
        assert count==expected
    return graphs,lookup,autos


GRAPHS,LOOKUP,AUT=atlas()
INDEX={key:i for i,key in enumerate(GRAPHS)}


def extend(a,alpha):
    b=[row+[x] for row,x in zip(a,alpha)]
    return b+[list(alpha)+[0]]


def count_system(params):
    v,k,lam,mu=params;rows=[];rhs=[]
    row=[Q(0)]*52;row[0]=Q(1);rows.append(row);rhs.append(Q(1))
    for key in GRAPHS:
        n,bits=key
        if n==5:continue
        a=matrix(*key)
        for subset in [()]+[(i,) for i in range(n)]+list(combinations(range(n),2)):
            if not subset:t=v-n
            elif len(subset)==1:t=k-sum(a[subset[0]])
            else:
                i,j=subset;t=(lam if a[i][j] else mu)-sum(a[i][z]*a[j][z] for z in range(n))
            row=[Q(0)]*52;row[INDEX[key]]=-AUT[key]*t
            for alpha in product((0,1),repeat=n):
                if not all(alpha[i] for i in subset):continue
                child=LOOKUP[n+1,code(extend(a,alpha))]
                row[INDEX[child]]+=AUT[child]
            rows.append(row);rhs.append(Q(0))
    assert len(rows)==160
    return rows,rhs


def solve_counts(rows,rhs):
    m=[list(map(Q,row))+[b] for row,b in zip(rows,rhs)];pivot=0;pivots=[]
    for c in range(52):
        j=next((j for j in range(pivot,len(m)) if m[j][c]),None)
        if j is None:continue
        m[pivot],m[j]=m[j],m[pivot];p=m[pivot][c]
        m[pivot]=[x/p for x in m[pivot]]
        for j in range(len(m)):
            if j!=pivot and m[j][c]:
                p=m[j][c];m[j]=[x-p*y for x,y in zip(m[j],m[pivot])]
        pivots.append(c);pivot+=1
    assert pivots==list(range(48))
    assert all(not any(row) for row in m[48:])
    affine=[[Q(0)]*5 for _ in range(52)]
    for i in range(48):affine[i]=[m[i][-1]]+[-m[i][j] for j in range(48,52)]
    for i in range(4):affine[48+i][i+1]=Q(1)
    for row,b in zip(rows,rhs):
        assert [sum(row[i]*affine[i][j] for i in range(52)) for j in range(5)]==[b,0,0,0,0]
    return affine


def inverse(a):
    n=len(a);m=[list(map(Q,row))+[Q(i==j) for j in range(n)] for i,row in enumerate(a)]
    piv=[]
    for i in range(n):
        p=m[i][i];assert p>0;piv.append(p);m[i]=[x/p for x in m[i]]
        for j in range(n):
            if j!=i:
                p=m[j][i];m[j]=[x-p*y for x,y in zip(m[j],m[i])]
    inv=[row[n:] for row in m]
    assert all(sum(a[i][k]*inv[k][j] for k in range(n))==int(i==j) for i in range(n) for j in range(n))
    return inv,piv


def kernel_row(params,dimension,a,b,term):
    bits,ell,den,integers,weight=term;H=matrix(3,bits);f=[Q(x,den) for x in integers]
    gram=[[Q(1) if i==j else b if H[i][j] else a for j in range(3)] for i in range(3)]
    inv,piv=inverse(gram);D=dimension-3;assert D>1
    profiles=list(product((0,1),repeat=3));uv=[[b if bit else a for bit in p] for p in profiles]
    project=[[sum(u[i]*inv[i][j]*v[j] for i in range(3) for j in range(3)) for v in uv] for u in uv]
    heights=[1-project[i][i] for i in range(8)]
    row=[Q(0)]*52
    # This implementation enumerates embeddings into each host graph,
    # rather than appending profiles and looking up automorphism ratios.
    for idx,(n,hostbits) in enumerate(GRAPHS):
        if n not in (4,5):continue
        host=matrix(n,hostbits);value=Q(0)
        for root in permutations(range(n),3):
            if any(host[root[i]][root[j]]!=H[i][j] for i,j in combinations(range(3),2)):continue
            outside=[x for x in range(n) if x not in root]
            for x in outside:
                for y in outside:
                    if len(set(root+(x,y)))!=n:continue
                    px=sum(host[root[i]][x]<<(2-i) for i in range(3))
                    py=sum(host[root[i]][y]<<(2-i) for i in range(3))
                    if x==y:z=heights[px]
                    else:z=(b if host[x][y] else a)-project[px][py]
                    kernel=z if ell==1 else (D*z*z-heights[px]*heights[py])/(D-1)
                    value+=f[px]*f[py]*kernel
        row[idx]=value/AUT[3,bits]
    return row,{'anchor':bits,'degree':ell,'positive_elimination_pivots':list(map(str,piv))}


CASES=[
 {'parameters':[1344,221,88,26],'dimension':56,'a':'-1/17','b':'5/17',
  'terms':[
   ['001',1,30,[30,17,17,0,24,0,0,0],11],
   ['011',2,30,[21,7,20,30,20,30,30,0],9],
   ['011',1,30,[20,-7,30,0,30,0,0,0],2],
   ['011',1,30,[21,-6,30,0,30,0,0,0],4],
   ['001',1,30,[30,17,17,0,25,0,0,0],1]],
  'expected':['-735463695439/4904219511462297600','-20935675404/171401572225','-246307918/3986083075','-45666157858/171401572225','-946035944/1105816595']},
 {'parameters':[540,147,66,30],'dimension':35,'a':'-5/49','b':'13/49',
  'terms':[
   ['001',1,10,[10,6,6,0,7,0,0,0],102],
   ['111',2,10,[3,-4,-4,10,-4,10,10,0],80],
   ['011',1,10,[10,4,9,0,9,0,0,0],61],
   ['011',1,10,[10,5,8,0,8,0,0,0],37],
   ['011',1,10,[10,4,8,0,8,0,0,0],22]],
  'expected':['-32843546233/17050267530000000','-519736215717/31197753650000','-132260706/227388875','-8705719132623/15598876825000','-1025687358819/623955073000']}
]


def main():
    reports=[]
    for case in CASES:
        params=case['parameters'];v,k,lam,mu=params;dim=case['dimension'];a=Q(case['a']);b=Q(case['b'])
        assert k*(k-lam-1)==(v-k-1)*mu
        # Exact Bose-Mesner multiplication verifies the proposed Gram.
        ci=1-a;ca=b-a;cj=a
        coeff=[ci*ci+ca*ca*(k-mu),2*ci*ca+ca*ca*(lam-mu),ca*ca*mu+2*ci*cj+2*k*ca*cj+v*cj*cj]
        assert coeff==[Q(v,dim)*x for x in (ci,ca,cj)] and ci+k*ca+v*cj==0
        rows,rhs=count_system(params);affine=solve_counts(rows,rhs)
        total=[Q(0)]*52;terms=[]
        for term in case['terms']:
            row,diagnostic=kernel_row(params,dim,a,b,term)
            expanded=[sum(row[i]*affine[i][j] for i in range(52)) for j in range(5)]
            expanded[0]/=v**4
            terms.append({**diagnostic,'weighted_direction':term,'L_expansion':list(map(str,expanded))})
            for i in range(52):total[i]+=term[-1]*row[i]
        final=[sum(total[i]*affine[i][j] for i in range(52)) for j in range(5)];final[0]/=v**4
        assert final==list(map(Q,case['expected'])) and all(x<0 for x in final)
        report={'verified':True,'parameters':params,'dimension':dim,'graph_counts_by_order':[1,2,4,11,34],
                'equations':160,'rank':48,'free_graphs':GRAPHS[48:],'normalization':'L=v^-5 times embedding-averaged ordered kernel sum; u=N_F/v^5; X=N_F/v.',
                'terms':terms,'strict_negative_final_coefficients':list(map(str,final)),
                'all_count_affine_rows':[{'n':n,'bits':bits,'coefficients':list(map(str,row))} for (n,bits),row in zip(GRAPHS,affine)],
                'weighted_kernel_raw_coefficients':list(map(str,total))}
        reports.append(report)
        print(json.dumps({k:report[k] for k in ('verified','parameters','equations','rank','strict_negative_final_coefficients')}),flush=True)
    target=Path(__file__).resolve().with_name('verification_1344_reconstructed.json')
    target.write_text(json.dumps(reports,indent=2)+'\n')


if __name__=='__main__':main()
