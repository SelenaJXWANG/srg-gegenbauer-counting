#!/usr/bin/env python3
"""Exact homogeneous projected Gegenbauer kernels and a profile exporter.

Importing this module performs no search, optimization, or file write.
Uses only Python's standard library. D is the residual ambient dimension.
The harmonic recurrence requires D >= 2; do not substitute D=0 or D=1.

Example:
  python3 projected_harmonic_kernels_20261003.py --r 1 --mu 3 \
      --anchor-bits 000 --max-degree 6 --output kg62_kernel_table.json
"""

from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import argparse
import json


def kernel(degree, D, z, hx, hy):
    """Return K_degree^D(z,hx,hy), normalized by K_l(h,h,h)=h**l.

    z=<qx,qy>, hx=||qx||^2, hy=||qy||^2. Arithmetic evaluation alone
    does not assert the input pair is realizable. In particular, an
    actual residual pair must have hx,hy>=0 and z*z<=hx*hy.
    """
    if not isinstance(degree, int) or degree < 0:
        raise ValueError("degree must be a nonnegative integer")
    if not isinstance(D, int) or D < 2:
        raise ValueError("Gegenbauer residual dimension D must be an integer >=2")
    z, s = F(z), F(hx)*F(hy)
    older, previous = F(1), z
    if degree == 0:
        return older
    for ell in range(2, degree+1):
        older, previous = previous, (
            (2*ell+D-4)*z*previous-(ell-1)*s*older
        )/(ell+D-3)
    return previous


def positive_inverse(a):
    """Exact inverse, verifying all positive Schur pivots of a symmetric Q."""
    n = len(a)
    if any(len(row) != n for row in a) or any(a[i][j] != a[j][i]
                                             for i in range(n) for j in range(n)):
        raise ValueError("anchor Gram must be square and symmetric")
    rows = [list(map(F,row))+[F(i == j) for j in range(n)]
            for i,row in enumerate(a)]
    pivots = []
    for i in range(n):
        pivot = rows[i][i]
        if pivot <= 0:
            raise ValueError("anchor Gram is not positive definite; use an independent anchor basis")
        pivots.append(pivot)
        rows[i] = [x/pivot for x in rows[i]]
        for j in range(n):
            if j != i:
                factor = rows[j][i]
                rows[j] = [x-factor*y for x,y in zip(rows[j],rows[i])]
    return [row[n:] for row in rows], pivots


def projected_kernel_table(r, mu, anchor_adjacency, max_degree=6):
    """Rational profile table for a least--3 SRG and a PD anchor.

    The input describes a possible anchor, not proof it occurs. A table
    record with base_schur_admissible=false is an impossible residual
    pair, not a pair that may be assigned a positive count.
    """
    if r <= 0 or mu <= 0 or not isinstance(r,int) or not isinstance(mu,int):
        raise ValueError("r and mu must be positive integers")
    t = len(anchor_adjacency)
    if any(len(row) != t for row in anchor_adjacency):
        raise ValueError("anchor adjacency must be square")
    if any(anchor_adjacency[i][j] not in (0,1) or
           anchor_adjacency[i][j] != anchor_adjacency[j][i]
           for i in range(t) for j in range(t)) or any(anchor_adjacency[i][i] for i in range(t)):
        raise ValueError("anchor adjacency must be a simple graph")
    v = F((mu+2*r)*(mu+3*r+3),mu)
    k, lam = 3*r+mu, r+mu-3
    f = (3*(v-1)-k)/(r+3)
    if v.denominator != 1 or f.denominator != 1:
        raise ValueError("SRG order and multiplicity must be integral")
    D = int(f)-t
    if D < 2:
        raise ValueError("this exporter needs residual dimension at least two")
    alpha, beta = -F(mu,2*(mu+3*r)), F(r,mu+3*r)
    Q = [[F(1) if i == j else beta if anchor_adjacency[i][j] else alpha
          for j in range(t)] for i in range(t)]
    inv,pivots = positive_inverse(Q)
    profiles = list(product((0,1),repeat=t))
    gs = [[beta if bit else alpha for bit in profile] for profile in profiles]
    projection = [[sum((gx[i]*inv[i][j]*gy[j] for i in range(t) for j in range(t)),F(0))
                   for gy in gs] for gx in gs]
    heights = [1-projection[i][i] for i in range(len(profiles))]
    records = []
    for x in range(len(profiles)):
        for y in range(len(profiles)):
            for adjacency in (0,1):
                z = (beta if adjacency else alpha)-projection[x][y]
                hx,hy = heights[x],heights[y]
                records.append({"profile_x":x,"profile_y":y,"adjacency":adjacency,
                    "z":z,"hx":hx,"hy":hy,
                    "base_schur_admissible":hx >= 0 and hy >= 0 and z*z <= hx*hy,
                    "kernels":[kernel(ell,D,z,hx,hy) for ell in range(max_degree+1)]})
        h = heights[x]
        records.append({"profile_x":x,"profile_y":x,"adjacency":"diagonal",
            "z":h,"hx":h,"hy":h,"base_schur_admissible":h >= 0,
            "kernels":[kernel(ell,D,h,h,h) for ell in range(max_degree+1)]})
    return {"parameters":[int(v),k,lam,mu],"r":r,"f":int(f),"D":D,
            "anchor_adjacency":anchor_adjacency,"anchor_Gram":Q,
            "anchor_positive_pivots":pivots,"profiles":profiles,"heights":heights,
            "max_degree":max_degree,"records":records}


def encode(obj):
    if isinstance(obj,F):
        return obj.numerator if obj.denominator == 1 else str(obj)
    if isinstance(obj,dict):
        return {str(k):encode(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)):
        return [encode(v) for v in obj]
    return obj


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r',type=int,required=True)
    parser.add_argument('--mu',type=int,required=True)
    parser.add_argument('--anchor-bits',required=True,
                        help='upper-triangle bits in lexicographic pair order; e.g. 000')
    parser.add_argument('--max-degree',type=int,default=6)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    bits = args.anchor_bits
    t = 1
    while t*(t-1)//2 < len(bits):
        t += 1
    if t*(t-1)//2 != len(bits) or any(bit not in '01' for bit in bits):
        parser.error('anchor-bits must have triangular length and contain only 0/1')
    if args.max_degree < 0:
        parser.error('max-degree must be nonnegative')
    A = [[0]*t for _ in range(t)]
    for (i,j),bit in zip(combinations(range(t),2),bits):
        A[i][j] = A[j][i] = int(bit)
    data = projected_kernel_table(args.r,args.mu,A,args.max_degree)
    args.output.write_text(json.dumps(encode(data),indent=2)+'\n')
    print(f"Wrote degrees 0..{args.max_degree}, {len(data['records'])} exact records to {args.output}")


if __name__ == '__main__':
    main()
