#!/usr/bin/env python3
"""Independent standard-library audit of the pure 405 seven-point count core.

No optimizer or attack7 import. Verify the graph catalogue/mapping by every
vertex permutation, then clear denominators to verify every SRG extension
equation and the triangular elimination behind the saved exact recipes.
"""
from array import array
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import comb, lcm
from pathlib import Path
import json
import sys
import time

ROOT = Path(__file__).resolve().parent
DATA = ROOT/'data/degree2'
V, K, LAM, MU = 405, 132, 63, 33


def edge_bit(i,j):
    if i > j:
        i,j = j,i
    return 1 << (j*(j-1)//2+i)


def adjacency(n,mask):
    return [[int(i != j and bool(mask & edge_bit(i,j)))
             for j in range(n)] for i in range(n)]


def mask_of(a):
    return sum(edge_bit(i,j) for j in range(len(a)) for i in range(j) if a[i][j])


def extend(n,mask,profile):
    return mask | sum(edge_bit(i,n) for i in range(n) if profile[i])


def main():
    started = time.monotonic()
    graphs = [tuple(map(int,line.split())) for line in
              (ROOT/'atlas/graphs.txt').read_text().splitlines()]
    assert [i for i,n,m in graphs] == list(range(1252))
    maps = {}
    for n in range(1,8):
        mapping = array('h')
        mapping.frombytes((ROOT/f'atlas/map{n}.bin').read_bytes())
        if sys.byteorder != 'little':
            mapping.byteswap()
        assert len(mapping) == 1 << comb(n,2)
        maps[n] = mapping
    map_records = []
    for n in range(1,8):
        representatives = [(idx,mask) for idx,nn,mask in graphs if nn == n]
        edges = [(i,j) for j in range(n) for i in range(j)]
        seen = bytearray(len(maps[n]))
        orbit_counts = Counter()
        comparisons = 0
        # Three chunks of seven bits make exhaustive relabelling cheap,
        # while using a direct independently defined edge permutation.
        for perm in permutations(range(n)):
            weights = [edge_bit(perm[i],perm[j]) for i,j in edges]
            chunks = []
            for start in range(0,len(edges),7):
                ww = weights[start:start+7]
                table = [0]*(1 << len(ww))
                for bits in range(1,len(table)):
                    least = bits & -bits
                    table[bits] = table[bits ^ least] | ww[least.bit_length()-1]
                chunks.append(table)
            for idx,mask in representatives:
                image = 0
                for chunk,table in enumerate(chunks):
                    image |= table[(mask >> (7*chunk)) & (len(table)-1)]
                assert maps[n][image] == idx, ('bad atlas map',n,idx,mask,image)
                comparisons += 1
                if not seen[image]:
                    seen[image] = 1
                    orbit_counts[idx] += 1
        assert all(seen), ('uncovered labelled masks',n)
        assert sum(orbit_counts.values()) == len(maps[n])
        map_records.append({'n':n,'classes':len(representatives),
                            'labelled_masks':len(maps[n]),
                            'permutation_checks':comparisons})
        print('PASS exhaustive graph mapping n=',n,flush=True)

    raw = json.loads((DATA/'core_recipes_405.json').read_text())
    cores = raw['cores']
    expected = [idx for idx,n,mask in graphs if n == 1 or
                min(map(sum,adjacency(n,mask))) >= 3]
    assert cores == expected and len(cores) == 174
    recipes = {int(i):{int(j):F(value) for j,value in row.items()}
               for i,row in raw['expressions'].items()}
    assert sorted(recipes) == list(range(1252))
    den = lcm(*(value.denominator for row in recipes.values() for value in row.values()))
    integer = {i:{j:int(value*den) for j,value in row.items()}
               for i,row in recipes.items()}
    assert all(recipes[idx] == {j:F(1)} for j,idx in enumerate(cores))

    def verify_row(row):
        reduced = {}
        for idx,coef in row.items():
            for j,value in integer[idx].items():
                reduced[j] = reduced.get(j,0)+coef*value
        assert not any(reduced.values()), ('nonzero exact reduced row',row,reduced)

    def cn(a,selected):
        if not selected:
            return V-len(a)
        if len(selected) == 1:
            return K-sum(a[selected[0]])
        i,j = selected
        return (LAM if a[i][j] else MU)-sum(x*y for x,y in zip(a[i],a[j]))

    extension_rows = 0
    for idx,n,mask in graphs:
        if n == 7:
            continue
        a = adjacency(n,mask)
        profiles = list(product((0,1),repeat=n))
        target = [maps[n+1][extend(n,mask,p)] for p in profiles]
        subsets = [()]+[(i,) for i in range(n)]+list(combinations(range(n),2))
        for selected in subsets:
            row = {idx:-cn(a,selected)*2**n}
            for p,j in zip(profiles,target):
                if all(p[i] for i in selected):
                    row[j] = row.get(j,0)+(V-n)
            verify_row(row)
            extension_rows += 1
    print('PASS all exact SRG extension rows:',extension_rows,flush=True)

    rank_key = {idx:(n,-bin(mask).count('1')) for idx,n,mask in graphs}
    elimination_rows = 0
    for idx,n,mask in graphs:
        if idx in cores:
            continue
        a = adjacency(n,mask)
        x = min(range(n),key=lambda i:sum(a[i]))
        assert sum(a[x]) <= 2
        keep = [i for i in range(n) if i != x]
        h = [[a[i][j] for j in keep] for i in keep]
        hm = mask_of(h)
        hid = maps[n-1][hm]
        selected = tuple(i for i,j in enumerate(keep) if a[x][j])
        rest = [i for i in range(n-1) if i not in selected]
        row = {hid:-2**(n-1)*cn(h,selected)}
        singled = []
        for bits in product((0,1),repeat=len(rest)):
            p = [int(i in selected) for i in range(n-1)]
            for i,b in zip(rest,bits):
                p[i] = b
            j = maps[n][extend(n-1,hm,p)]
            row[j] = row.get(j,0)+(V-n+1)
            if not any(bits):
                assert j == idx
                singled.append(j)
            else:
                assert rank_key[j] < rank_key[idx]
        assert singled == [idx] and rank_key[hid] < rank_key[idx]
        verify_row(row)
        elimination_rows += 1
    assert elimination_rows == 1252-174
    # All recipes also match the degree-six run, independently of any PSD cut.
    other = ROOT/'data/degree6/core_recipes_405.json'
    assert json.loads(other.read_text()) == raw
    result = {'status':'PASS','parameters':[V,K,LAM,MU],
              'variables':1252,'free_core_variables':len(cores),
              'exhaustive_mapping':map_records,
              'exact_extension_rows':extension_rows,
              'exact_triangular_elimination_rows':elimination_rows,
              'common_recipe_denominator':den,
              'all_extension_rows_reduce_identically_to_zero':True,
              'degree2_degree6_recipes_identical':True,
              'uniform_count_upper_bound':'y_H <= 2^binom(|H|,2)',
              'upper_bound_proof':'The empty extension identity is sum_alpha y_(H,alpha)=2^|H| y_H. Nonnegativity bounds each summand; induction from y_K1=1 proves the displayed bound.',
              'elapsed_seconds':round(time.monotonic()-started,3),
              'scope':'Audit of exact count coordinates and atlas mapping only; no standalone nonexistence claim.'}
    target = Path(__file__).with_suffix('.json')
    target.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS independent exact core audit',target,flush=True)


if __name__ == '__main__':
    main()
