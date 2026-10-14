#!/usr/bin/env python3
"""All-degree FIVE-POINT relaxation feasibility for SRG(405,132,63,33).

Standard library only. Reconstructs graph counts, profile pair masses,
all low-degree PSD blocks, and the exact uniform tail certificate.
No optimizer, root-free fact, forbidden graph, or existing exclusion
certificate is used. This does not establish existence of an SRG.
"""

from fractions import Fraction as F
from itertools import combinations, permutations
from math import comb
from pathlib import Path
import importlib.util
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE
spec = importlib.util.spec_from_file_location('exact_count_reconstruction',
    ROOT/'verify_1344_reconstructed.py')
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
from projected_harmonic_kernels_20261003 import projected_kernel_table

PARAMETERS = (405,132,63,33)
FREE_COUNTS = list(map(F,[988060,134777,287429,68038]))


def values(D,z,hx,hy,last):
    out = [F(1),z]
    for degree in range(2,last+1):
        out.append(((2*degree+D-4)*z*out[-1]
                    -(degree-1)*hx*hy*out[-2])/(degree+D-3))
    return out[:last+1]


def exact_psd(matrix):
    a = [row[:] for row in matrix]
    assert all(a[i][j] == a[j][i] for i in range(len(a)) for j in range(len(a)))
    pivots = []
    for i in range(len(a)):
        p = a[i][i]
        assert p >= 0, ('negative pivot',i,p)
        if p == 0:
            assert all(a[i][j] == 0 for j in range(i+1,len(a))),('nonzero zero-pivot row',i)
            continue
        pivots.append(p)
        for j in range(i+1,len(a)):
            for k in range(i+1,len(a)):
                a[j][k] -= a[j][i]*a[i][k]/p
    return pivots


def tail(D,s,q):
    moment,result = F(1),F(0)
    for j in range(q+1):
        if j:
            moment *= F(2*j-1,D+2*j-3)
        result += comb(q,j)*s**(q-j)*(1-s)**j*moment
    return result


def encode(value):
    if isinstance(value,F):
        return str(value)
    if isinstance(value,dict):
        return {str(k):encode(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):
        return [encode(v) for v in value]
    return value


def main():
    rows,rhs = base.count_system(PARAMETERS)
    affine = base.solve_counts(rows,rhs)
    ex = [F(1)]+FREE_COUNTS
    counts = [sum((a*b for a,b in zip(row,ex)),F(0)) for row in affine]
    assert len(counts) == 52 and all(x >= 0 for x in counts)
    assert counts[0] == 1
    assert all(sum((a*b for a,b in zip(row,counts)),F(0)) == rhs[i]
               for i,row in enumerate(rows))
    record = {'status':'PASS','parameters':PARAMETERS,
        'statement':'The real five-point count/profile-kernel relaxation is feasible at EVERY harmonic degree.',
        'scope':'standard r-eigenspace embedding; all anchors of size at most three; ordinary SRG extension equations; nonnegative real graph counts',
        'not_claimed':['existence of an SRG','feasibility of integer graph counts',
                       'feasibility with counts on six or more vertices','feasibility with additional structural constraints'],
        'free_graphs':base.GRAPHS[48:],'free_counts':FREE_COUNTS,
        'all_52_counts':[{'n':n,'bits':bits,'X':x} for (n,bits),x in zip(base.GRAPHS,counts)],
        'equations_verified':len(rows),'anchor_records':[]}
    for t,bits in base.GRAPHS:
        if t > 3:
            break
        H = base.matrix(t,bits)
        table = projected_kernel_table(33,33,H,0)
        profiles = table['profiles']
        d,D,N = len(profiles),30-t,405-t
        pair_lookup = {(row['profile_x'],row['profile_y'],row['adjacency']):row
                       for row in table['records']}
        one = [F(0)]*d
        pairs = [[[F(0),F(0)] for _ in range(d)] for _ in range(d)]
        # This assembles counts directly from host embeddings and the
        # 52 witness counts, without the optimizer's affine kernel blocks.
        for host_idx,(n,hostbits) in enumerate(base.GRAPHS):
            if n not in (t+1,t+2):
                continue
            host = base.matrix(n,hostbits)
            amount = counts[host_idx]/base.AUT[t,bits]
            for root in permutations(range(n),t):
                if any(host[root[i]][root[j]] != H[i][j]
                       for i,j in combinations(range(t),2)):
                    continue
                outside = [x for x in range(n) if x not in root]
                if n == t+1:
                    x = outside[0]
                    px = sum(host[root[i]][x] << (t-1-i) for i in range(t))
                    one[px] += amount
                else:
                    for x,y in (outside,outside[::-1]):
                        px = sum(host[root[i]][x] << (t-1-i) for i in range(t))
                        py = sum(host[root[i]][y] << (t-1-i) for i in range(t))
                        pairs[px][py][host[x][y]] += amount
        assert all(x >= 0 for x in one)
        assert all(pairs[i][j][e] >= 0 and pairs[i][j][e] == pairs[j][i][e]
                   for i in range(d) for j in range(d) for e in (0,1))
        row_masses = [sum((pairs[i][j][e] for j in range(d) for e in (0,1)),F(0))
                      for i in range(d)]
        assert row_masses == [(N-1)*x for x in one]
        heights = table['heights']
        assert all(h > 0 for h in heights)
        s = max(row['z']**2/(row['hx']*row['hy']) for row in table['records']
                if row['adjacency'] in (0,1))
        assert s < 1
        q = 0
        while (N-1)*tail(D,s,q) > 1:
            q += 1
        cutoff = 2*q
        kernel_values = {(i,j,e):values(D,pair_lookup[i,j,e]['z'],heights[i],heights[j],cutoff-1)
                         for i in range(d) for j in range(d) for e in (0,1)}
        checks = []
        for ell in range(cutoff):
            matrix = [[(one[i]*heights[i]**ell if i == j else F(0))
                        +sum((pairs[i][j][e]*kernel_values[i,j,e][ell] for e in (0,1)),F(0))
                       for j in range(d)] for i in range(d)]
            pivots = exact_psd(matrix)
            checks.append({'degree':ell,'rank':len(pivots),'positive_pivots':pivots})
        record['anchor_records'].append({'anchor_size':t,'anchor_bits':bits,'D':D,'N':N,
            'anchor_positive_pivots':table['anchor_positive_pivots'],
            'profiles':profiles,'one_point_masses':one,'ordered_distinct_pair_masses':pairs,
            'row_masses':row_masses,'row_identity_verified':True,'residual_heights':heights,
            'max_squared_angle':s,'automatic_degrees_from':cutoff,
            'tail_row_bound':(N-1)*tail(D,s,q),'low_degree_exact_PSD_checks':checks})
        print('PASS anchor',t,bits,'low degrees',cutoff,'tail from',cutoff,flush=True)
    # Unanchored full Bose--Mesner spectra, for completeness. These are
    # determined by the parameters and do not depend on the primal point.
    ua,ub = values(30,-F(1,8),1,1,5),values(30,F(1,4),1,1,5)
    spectra = [[1+132*b+272*a,1+33*b-34*a,1-3*b+2*a] for a,b in zip(ua,ub)]
    assert all(x >= 0 for spectrum in spectra for x in spectrum)
    assert 404*tail(30,F(1,16),3) <= 1
    record['unanchored'] = {'degrees_0_through_5_spectra':spectra,
        'automatic_degrees_from':6,'max_squared_angle':F(1,16),
        'tail_row_bound':404*tail(30,F(1,16),3)}
    record['all_degrees_feasible'] = True
    target = Path(__file__).with_suffix('.json')
    target.write_text(json.dumps(encode(record),indent=2)+'\n')
    print('PASS: EXACT ALL-DEGREE FIVE-POINT RELAXATION WITNESS',target,flush=True)


if __name__ == '__main__':
    main()
