"""Independent rational audit of active PSD and unique-completion sources.

Uses Python standard library only. Does not import the discovery/recovery
code. Graphs are classified by an independent isomorphism backtrack
against the saved catalogue; no binary orbit-map cache is used.
The full-shell theorem is a mathematical prerequisite, audited separately.
Affine elimination, integer-moment rows and final combination are outside
this verifier and are checked by the companion source audit.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import reduce
from hashlib import sha256
from itertools import product
from math import gcd, lcm
from pathlib import Path
import json
import sys
import time

HERE = Path(__file__).resolve().parent
PORTABLE = (HERE / 'atlas_graphs.txt').is_file()
BASE = HERE if PORTABLE else HERE.parent
CATALOGUE = HERE / 'atlas_graphs.txt' if PORTABLE else BASE / 'work/atlas/graphs.txt'
DATA = HERE if PORTABLE else BASE / 'output/paper/support'


def signature(a):
    degrees = list(map(sum, a))
    return len(a), tuple(sorted((degrees[i], tuple(sorted(degrees[j] for j in range(len(a)) if a[i][j])))
                               for i in range(len(a))))


def decode(n, mask):
    a = [[0] * n for _ in range(n)]
    for j in range(n):
        for i in range(j):
            a[i][j] = a[j][i] = (mask >> (j * (j - 1) // 2 + i)) & 1
    return a


def encode(a):
    return sum(a[i][j] << (j * (j - 1) // 2 + i)
               for j in range(len(a)) for i in range(j))


def isomorphic(a, b):
    n = len(a)
    if len(b) != n:
        return False
    da = list(map(sum, a)); db = list(map(sum, b))
    ka = [(da[i], tuple(sorted(da[j] for j in range(n) if a[i][j]))) for i in range(n)]
    kb = [(db[i], tuple(sorted(db[j] for j in range(n) if b[i][j]))) for i in range(n)]
    if Counter(ka) != Counter(kb):
        return False
    candidates = [[j for j in range(n) if ka[i] == kb[j]] for i in range(n)]
    order = sorted(range(n), key=lambda i: (len(candidates[i]), -da[i]))
    assigned = {}
    used = set()
    def search(pos):
        if pos == n:
            return True
        i = order[pos]
        for j in candidates[i]:
            if j in used or any(a[i][ii] != b[j][jj] for ii, jj in assigned.items()):
                continue
            assigned[i] = j; used.add(j)
            if search(pos + 1):
                return True
            used.remove(j); del assigned[i]
        return False
    return search(0)


def primitive(row):
    row = {int(i): F(q) for i, q in row.items() if q}
    if not row:
        return {}
    den = lcm(*(q.denominator for q in row.values()))
    integers = {i: int(q * den) for i, q in row.items()}
    common = reduce(gcd, map(abs, integers.values()))
    if integers[min(integers)] < 0:
        common = -common
    return {i: q // common for i, q in integers.items()}


def norm(a, h):
    return sum(3 * hi * hi for hi in h) + 2 * sum(a[i][j] * h[i] * h[j]
               for j in range(len(h)) for i in range(j))


def gram_vector(a, h):
    return [3 * h[i] + sum(a[i][j] * h[j] for j in range(len(h))) for i in range(len(h))]


def extend(a, alpha, beta=None, edge=0):
    n = len(a); size = n + 1 + int(beta is not None)
    b = [[0] * size for _ in range(size)]
    for i in range(n):
        for j in range(n):
            b[i][j] = a[i][j]
        b[i][n] = b[n][i] = alpha[i]
        if beta is not None:
            b[i][n+1] = b[n+1][i] = beta[i]
    if beta is not None:
        b[n][n+1] = b[n+1][n] = edge
    return b


def positive_inverse(a):
    """Positive Schur pivots prove PD; Gauss--Jordan obtains exact inverse."""
    n = len(a); schur = [[F(x) for x in row] for row in a]
    pivots = []
    for i in range(n):
        pivot = schur[i][i]
        assert pivot > 0, ('anchor is not positive definite', i, pivot)
        pivots.append(pivot)
        for j in range(i+1, n):
            for k in range(i+1, n):
                schur[j][k] -= schur[j][i] * schur[i][k] / pivot
    aug = [[F(x) for x in row] + [F(i == j) for j in range(n)] for i, row in enumerate(a)]
    for i in range(n):
        pivot = aug[i][i]
        aug[i] = [q / pivot for q in aug[i]]
        for j in range(n):
            if i != j:
                multiplier = aug[j][i]
                aug[j] = [q - multiplier * r for q, r in zip(aug[j], aug[i])]
    inv = [row[n:] for row in aug]
    assert all(sum(a[i][k] * inv[k][j] for k in range(n)) == (i == j)
               for i in range(n) for j in range(n))
    return inv, pivots


class Audit:
    def __init__(self, v):
        self.v = v
        self.k, self.f = {375: (102, 34), 441: (88, 44), 575: (112, 46)}[v]
        self.eta = F(self.k + 3, v)
        self.out = DATA / f'rebuild_{v}'
        if not self.out.is_dir() and PORTABLE:
            self.out = HERE.parent / f'support/rebuild_{v}'
        self.graphs = {}
        for line in CATALOGUE.read_text().splitlines():
            i, n, mask = map(int, line.split())
            self.graphs[i] = decode(n, mask)
        assert len(self.graphs) == 1252
        self.catalogue_classes = defaultdict(list)
        for i, a in self.graphs.items():
            key = signature(a)
            assert all(not isomorphic(a, self.graphs[j]) for j in self.catalogue_classes[key]), ('duplicate graph type', i)
            self.catalogue_classes[key].append(i)
        self.checked_lookup = {}
        self.forbidden = set()
        self.geometry_cache = {}

    def identify(self, a):
        key = len(a), encode(a)
        if key not in self.checked_lookup:
            matches = [i for i in self.catalogue_classes[signature(a)] if isomorphic(a, self.graphs[i])]
            assert len(matches) == 1, ('graph classification is not unique', key, matches)
            target = matches[0]
            self.checked_lookup[key] = target
        return self.checked_lookup[key]

    def roots(self):
        records = json.loads((self.out / f'forbidden_{self.v}.json').read_text())
        for w in records:
            a = self.graphs[w['id']]
            assert encode(a) == w['mask'] and len(a) == w['n']
            assert w['root'] is not None and all(isinstance(h, int) for h in w['root'])
            assert norm(a, w['root']) in (1, 2)
            self.forbidden.add(w['id'])
        return len(records)

    def completions(self):
        data = json.loads((self.out / 'completion_constraints.json').read_text())
        equations = [primitive(r) for r in data['eq']]
        covered = set(); types = Counter()
        for w in data['witnesses']:
            a = self.graphs[w['anchor']]; h = w['h']; n = len(a)
            assert n == w['n'] and encode(a) == w['anchor_mask']
            assert len(h) == n and all(isinstance(x, int) for x in h)
            assert sum(h) == 1 and norm(a, h) == 3
            z = gram_vector(a, h)
            assert z == w['profile'] and all(x in (0, 1) for x in z)
            b = extend(a, z)
            if w['type'] == 'unique_profile_completion':
                assert n <= 6
                target = self.identify(b)
                assert target == w['completed_graph']
                row = primitive({w['anchor']: -1, target: F(self.v-n, 2**n)})
            else:
                assert w['type'] == 'seven_card_bijection' and n == 7
                assert encode(b) == w['completed_8_mask']
                null = h + [-1]
                assert w['null'] == null and sum(null) == 0 and gram_vector(b, null) == [0] * 8
                deleted = w['deleted_vertex']
                assert null[deleted] in (-1, 1)
                keep = [i for i in range(8) if i != deleted]
                assert keep == w['card_vertices']
                card = [[b[i][j] for j in keep] for i in keep]
                target = self.identify(card)
                assert target == w['card']
                hh = [-null[i] // null[deleted] for i in keep]
                assert hh == w['card_h'] and sum(hh) == 1 and norm(card, hh) == 3
                assert gram_vector(card, hh) == [b[i][deleted] for i in keep]
                row = primitive({w['anchor']: 1, target: -1})
            assert row == equations[w['row_index']], ('completion row mismatch', w)
            covered.add(w['row_index']); types[w['type']] += 1
        assert covered == set(range(len(equations)))
        return {'equations': len(equations), 'witnesses': len(data['witnesses']), 'types': dict(types)}

    def geometry(self, anchor):
        if anchor in self.geometry_cache:
            return self.geometry_cache[anchor]
        a = self.graphs[anchor]; t = len(a)
        gram = [[F(a[i][j] + 3 * (i == j)) - self.eta for j in range(t)] for i in range(t)]
        inv, pivots = positive_inverse(gram)
        patterns = []; vectors = []; heights = []; one = []
        for alpha in product((0, 1), repeat=t):
            target = self.identify(extend(a, alpha))
            vector = [F(x) - self.eta for x in alpha]
            h = 3 - self.eta - sum(vector[i] * inv[i][j] * vector[j] for i in range(t) for j in range(t))
            if target in self.forbidden or h < 0:
                continue
            patterns.append(alpha); vectors.append(vector); heights.append(h); one.append(target)
        projected = [[sum(x[k] * inv[k][l] * y[l] for k in range(t) for l in range(t))
                      for y in vectors] for x in vectors]
        two = [[[self.identify(extend(a, alpha, beta, edge)) for edge in (0, 1)]
                for beta in patterns] for alpha in patterns]
        result = t, patterns, heights, one, projected, two, pivots
        self.geometry_cache[anchor] = result
        return result

    def psd_rows(self):
        cert = json.loads((self.out / 'recovered_exact_certificate.json').read_text())
        count = 0; degree_counts = Counter()
        for index, item in enumerate(cert['rows']):
            src = item['source']
            if src['kind'] != 'PSD':
                continue
            anchor, degree = src['anchor'], src['degree']
            weights = [F(x) for x in src['direction']]
            t, patterns, heights, one, projected, two, pivots = self.geometry(anchor)
            assert degree in (0, 1, 2) and len(weights) == len(patterns)
            D = self.f - t
            assert D > 1
            result = {}
            def add(target, value):
                if value and target not in self.forbidden:
                    result[target] = result.get(target, F(0)) - value
            offscale = F(self.v-t-1, 2**(t+1))
            for i, wi in enumerate(weights):
                diag = [F(1), heights[i], heights[i]**2 * (1-F(1, D))][degree]
                add(one[i], wi**2 * diag)
                for j, wj in enumerate(weights):
                    for edge in (0, 1):
                        z = F(edge) - self.eta - projected[i][j]
                        kernel = [F(1), z, z*z-heights[i]*heights[j]/D][degree]
                        add(two[i][j][edge], offscale * wi * wj * kernel)
            result = {i: q for i, q in result.items() if q}
            expected = {int(i): F(q) for i, q in item['full_y_row'].items() if F(q)}
            assert result == expected, ('active PSD row mismatch', index, src)
            count += 1; degree_counts[degree] += 1
        return {'active_rows': count, 'degrees': dict(degree_counts),
                'anchors': len(self.geometry_cache),
                'positive_schur_pivots': {str(h): list(map(str, v[-1])) for h, v in self.geometry_cache.items()}}


def main(v):
    begin = time.monotonic(); audit = Audit(v)
    roots = audit.roots(); completion = audit.completions(); psd = audit.psd_rows()
    files = [audit.out / f'forbidden_{v}.json', audit.out / 'completion_constraints.json',
             audit.out / 'recovered_exact_certificate.json', CATALOGUE]
    result = {'v': v, 'status': 'PASS', 'arithmetic': 'standard-library fractions.Fraction',
              'norm_1_or_2_forbidden_witnesses': roots, 'completion': completion, 'PSD': psd,
              'independently_classified_graphs': len(audit.checked_lookup),
              'seconds': time.monotonic()-begin,
              'files_sha256': {str(p.relative_to(BASE) if p.is_relative_to(BASE) else p): sha256(p.read_bytes()).hexdigest() for p in files},
              'boundary': 'Verifies source PSD rows, norm witnesses, completion rows. Full-shell theorem is a separate mathematical prerequisite. Does not verify affine elimination, active integer moments, upper bound, or final weighted margin.'}
    path = DATA / f'independent_psd_completion_audit_{v}.json'
    path.write_text(json.dumps(result, indent=2))
    print(json.dumps({k: q for k, q in result.items() if k not in ('PSD', 'files_sha256')}, indent=2))
    print('PSD summary', {k: q for k, q in psd.items() if k != 'positive_schur_pivots'})


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 575)
