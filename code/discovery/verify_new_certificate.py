#!/usr/bin/env python3
"""Independently verify an SRG certificate in the 2026-10-14 article convention.

Only the Python standard library is required. No solver or discovery code is
imported. N_F always counts labelled induced embeddings, not vertex subsets.
All numerical inputs are converted to Fraction, and every test is exact.

Usage: python3 verify_certificate.py [data.json] [certificate.json]
"""
from argparse import ArgumentParser
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import comb, factorial, isqrt
from pathlib import Path
import json
import time


def check(condition, message):
    if not condition:
        raise ValueError(message)


def matrix(n, bits):
    check(len(bits) == n * (n - 1) // 2, "Invalid graph bit length")
    a = [[0] * n for _ in range(n)]
    for (i, j), bit in zip(combinations(range(n), 2), bits):
        check(bit in "01", "Invalid graph bit")
        a[i][j] = a[j][i] = int(bit)
    return a


def graph_code(a, labels):
    return "".join(str(a[i][j]) for i, j in combinations(labels, 2))


def enumerate_graphs(p):
    """Enumerate all labelled graphs, identify every permutation orbit."""
    graphs, arrays, autos, canonical = [], [], [], {}
    for n in range(1, p + 1):
        orders = list(permutations(range(n)))
        for bits_tuple in product("01", repeat=n * (n - 1) // 2):
            bits = "".join(bits_tuple)
            if (n, bits) in canonical:
                continue
            a = matrix(n, bits)
            orbit = [graph_code(a, labels) for labels in orders]
            check(min(orbit) == bits, "Graph enumeration order failed")
            index = len(graphs)
            for equivalent in orbit:
                canonical[n, equivalent] = index
            aut = orbit.count(bits)
            check(len(set(orbit)) * aut == factorial(n), "Orbit-stabilizer failed")
            graphs.append((n, bits))
            arrays.append(a)
            autos.append(aut)
    return graphs, arrays, autos, canonical


def counting_equations(parameters, p, graphs, arrays, autos, canonical):
    """Regenerate precisely (2.9) and (2.6), with no Aut(H) factors."""
    v, k, lam, mu = parameters
    rows, rhs = [], []
    size = len(graphs)
    for hi, (h, _) in enumerate(graphs):
        if h >= p:
            continue
        H = arrays[hi]
        base_sets = [()] + [(i,) for i in range(h)] + list(combinations(range(h), 2))
        for required in base_sets:
            if not required:
                external = v - h
            elif len(required) == 1:
                external = k - sum(H[required[0]])
            else:
                i, j = required
                external = (lam if H[i][j] else mu) - sum(H[i][w] * H[j][w] for w in range(h))
            row = [Q(0)] * size
            row[hi] = -external
            for profile in product((0, 1), repeat=h):
                if not all(profile[i] for i in required):
                    continue
                extension = [H[i][:] + [profile[i]] for i in range(h)]
                extension.append(list(profile) + [0])
                fi = canonical[h + 1, graph_code(extension, range(h + 1))]
                row[fi] += 1
            rows.append(row)
            rhs.append(Q(0))
    for n in range(1, p + 1):
        row = [Q(1, aut) if order == n else Q(0)
               for (order, _), aut in zip(graphs, autos)]
        rows.append(row)
        rhs.append(Q(comb(v, n)))
    return rows, rhs


def rref(rows, rhs):
    """Exact augmented elimination, retaining and rejecting inconsistent rows."""
    columns = len(rows[0])
    a = [[Q(value) for value in row] + [Q(b)] for row, b in zip(rows, rhs)]
    pivots = []
    for col in range(columns):
        target = len(pivots)
        candidate = next((i for i in range(target, len(a)) if a[i][col]), None)
        if candidate is None:
            continue
        a[target], a[candidate] = a[candidate], a[target]
        pivot = a[target][col]
        a[target] = [value / pivot for value in a[target]]
        for i, row in enumerate(a):
            if i != target and row[col]:
                scale = row[col]
                a[i] = [x - scale * y for x, y in zip(row, a[target])]
        pivots.append(col)
    check(all(any(row[:-1]) or not row[-1] for row in a), "Counting equations inconsistent")
    return a[:len(pivots)], pivots


def verify_affine(rows, rhs, affine, free_indices, rank):
    n, q = len(affine), len(free_indices)
    check(all(len(row) == q + 1 for row in affine), "Affine dimensions disagree")
    check(n - rank == q, "Affine space is not complete: dimension mismatch")
    check(len(set(free_indices)) == q, "Repeated free-coordinate index")
    for j, index in enumerate(free_indices):
        check(0 <= index < n, "Invalid free-coordinate index")
        check(affine[index] == [Q(0)] + [Q(j == k) for k in range(q)],
              "Free variables are not the declared original labelled counts")
    for i, (row, b) in enumerate(zip(rows, rhs)):
        check(len(row) == n, "Counting row dimension mismatch")
        actual = [sum(row[k] * affine[k][j] for k in range(n) if row[k]) for j in range(q + 1)]
        check(actual == [b] + [Q(0)] * q, "Counting/affine identity failed at row %d" % i)


def inverse_positive(a):
    """Symmetric positive definiteness via exact Schur pivots; return inverse."""
    n = len(a)
    schur = [row[:] for row in a]
    positive_pivots = []
    for k in range(n):
        pivot = schur[k][k]
        check(pivot > 0, "Anchor Gram matrix is not positive definite")
        positive_pivots.append(pivot)
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                schur[i][j] -= schur[i][k] * schur[k][j] / pivot
    augmented = [row[:] + [Q(i == j) for j in range(n)] for i, row in enumerate(a)]
    for k in range(n):
        pivot = augmented[k][k]
        augmented[k] = [value / pivot for value in augmented[k]]
        for i in range(n):
            if i != k:
                scale = augmented[i][k]
                augmented[i] = [x - scale * y for x, y in zip(augmented[i], augmented[k])]
    inv = [row[n:] for row in augmented]
    check(all(sum(a[i][k] * inv[k][j] for k in range(n)) == int(i == j)
              for i in range(n) for j in range(n)), "Gram inverse check failed")
    return inv, positive_pivots


def gegenbauer(degree, dimension):
    polynomials = [[Q(1)], [Q(0), Q(1)]]
    for ell in range(2, degree + 1):
        coefficients = [Q(0)] * (ell + 1)
        for power, value in enumerate(polynomials[-1]):
            coefficients[power + 1] += (2 * ell + dimension - 4) * value
        for power, value in enumerate(polynomials[-2]):
            coefficients[power] -= (ell - 1) * value
        coefficients = [value / (ell + dimension - 3) for value in coefficients]
        polynomials.append(coefficients)
    return polynomials[degree]


def kernel_data(H, ell, m, d, alpha, beta):
    anchor = matrix(m, H)
    gram = [[Q(1) if i == j else beta if anchor[i][j] else alpha
             for j in range(m)] for i in range(m)]
    inv, pivots = inverse_positive(gram)
    profiles = list(product((0, 1), repeat=m))
    vectors = [[beta if bit else alpha for bit in profile] for profile in profiles]
    projection = [[sum(u[i] * inv[i][j] * w[j] for i in range(m) for j in range(m))
                   for w in vectors] for u in vectors]
    rho = [1 - projection[i][i] for i in range(2 ** m)]
    coefficients = gegenbauer(ell, d - m)
    pair = []
    for i in range(2 ** m):
        row = []
        for j in range(2 ** m):
            values = []
            for edge in (0, 1):
                t = (beta if edge else alpha) - projection[i][j]
                values.append(sum(coefficient * t ** power * (rho[i] * rho[j]) ** ((ell - power) // 2)
                                  for power, coefficient in enumerate(coefficients) if coefficient))
            row.append(values)
        pair.append(row)
    return {"anchor": anchor, "profiles": profiles, "rho": rho, "pair": pair,
            "diagonal": [value ** ell for value in rho], "pivots": pivots}


def direct_profile_matrices(kernel, m, size, canonical):
    """Build article C independently through labelled extension embeddings.

    Fixing profiles and the external edge describes one labelled extension.
    Its number of embeddings is exactly N_F. This is equivalent to (3.6).
    """
    order = 2 ** m
    C = [[[Q(0) for _ in range(order)] for _ in range(order)] for _ in range(size)]
    H, profiles = kernel["anchor"], kernel["profiles"]
    for i, profile in enumerate(profiles):
        extension = [H[a][:] + [profile[a]] for a in range(m)]
        extension.append(list(profile) + [0])
        fi = canonical[m + 1, graph_code(extension, range(m + 1))]
        C[fi][i][i] += kernel["diagonal"][i]
        for j, other in enumerate(profiles):
            for edge in (0, 1):
                extension = [H[a][:] + [profile[a], other[a]] for a in range(m)]
                extension.append(list(profile) + [0, edge])
                extension.append(list(other) + [edge, 0])
                fi = canonical[m + 2, graph_code(extension, range(m + 2))]
                C[fi][i][j] += kernel["pair"][i][j][edge]
    check(all(A[i][j] == A[j][i] for A in C for i in range(order) for j in range(order)),
          "Direct-profile matrices are not symmetric")
    return C


def enumerate_support_weights(kernel, m, f, graphs, arrays, autos):
    """Literally enumerate (3.6): psi, full-support ordered pairs, /Aut(F)."""
    result = [Q(0)] * len(graphs)
    H = kernel["anchor"]
    for fi, ((n, _), F, aut) in enumerate(zip(graphs, arrays, autos)):
        if n not in (m + 1, m + 2):
            continue
        for psi in permutations(range(n), m):
            if any(H[i][j] != F[psi[i]][psi[j]] for i, j in combinations(range(m), 2)):
                continue
            rest = [a for a in range(n) if a not in psi]
            pairs = [(rest[0], rest[0])] if n == m + 1 else [(rest[0], rest[1]), (rest[1], rest[0])]
            for a, b in pairs:
                i = sum(F[psi[t]][a] * 2 ** (m - 1 - t) for t in range(m))
                j = sum(F[psi[t]][b] * 2 ** (m - 1 - t) for t in range(m))
                kval = kernel["diagonal"][i] if a == b else kernel["pair"][i][j][F[a][b]]
                result[fi] += f[i] * f[j] * kval / aut
    return result


def exact_nested(value):
    return [exact_nested(item) for item in value] if isinstance(value, list) else Q(value)


def validate_spectral_embedding(parameters, d, alpha, beta):
    v, k, lam, mu = parameters
    check(0 < k < v - 1, "Verifier expects a nontrivial SRG parameter set")
    check(k * (k - lam - 1) == (v - k - 1) * mu, "SRG parameter relation failed")
    discriminant = (lam - mu) ** 2 + 4 * (k - mu)
    root = isqrt(discriminant)
    check(root * root == discriminant and root > 0, "Expected rational restricted eigenvalues")
    r, s = Q(lam - mu + root, 2), Q(lam - mu - root, 2)
    multiplicity_r = (-k - s * (v - 1)) / (r - s)
    multiplicity_s = v - 1 - multiplicity_r
    check(all(x.denominator == 1 and x > 0 for x in (multiplicity_r, multiplicity_s)),
          "Invalid restricted multiplicities")
    gram_eigenvalues = [1 + beta * k + alpha * (v - k - 1),
                        1 - alpha + (beta - alpha) * r,
                        1 - alpha + (beta - alpha) * s]
    multiplicities = [Q(1), multiplicity_r, multiplicity_s]
    check(alpha != beta and all(value >= 0 for value in gram_eigenvalues), "Spherical Gram is not PSD")
    check(sum(mult for mult, value in zip(multiplicities, gram_eigenvalues) if value > 0) == d,
          "Spherical Gram rank differs from d")
    return {"restricted_eigenvalues": [str(r), str(s)],
            "restricted_multiplicities": [str(multiplicity_r), str(multiplicity_s)],
            "gram_eigenvalues": list(map(str, gram_eigenvalues))}


def verify(data, certificate):
    started = time.monotonic()
    parameters = data["parameters"]
    d, m, p = data["d"], data["m"], data["p"]
    anchor_sizes = data.get("anchor_sizes", [m])
    alpha, beta = Q(data["alpha"]), Q(data["beta"])
    check(1 <= m <= d - 2 and m + 2 <= p <= parameters[0], "Invalid article dimensions")
    spectral = validate_spectral_embedding(parameters, d, alpha, beta)
    graphs, arrays, autos, canonical = enumerate_graphs(p)
    check(data["graphs"] == [{"n": n, "bits": bits, "aut": aut}
                              for (n, bits), aut in zip(graphs, autos)],
          "Saved graph list or automorphism orders differ from independent enumeration")
    rows, rhs = counting_equations(parameters, p, graphs, arrays, autos, canonical)
    reduced, pivots = rref(rows, rhs)
    affine, free_indices = exact_nested(data["affine"]), data["free_indices"]
    verify_affine(rows, rhs, affine, free_indices, len(pivots))
    saved_A, saved_b = exact_nested(data["A"]), exact_nested(data["b"])
    check(len(saved_A) == len(saved_b) and saved_A, "Saved counting dimensions mismatch")
    saved_reduced, saved_pivots = rref(saved_A, saved_b)
    check(saved_reduced == reduced and saved_pivots == pivots,
          "Saved counting system is not equivalent to the complete article system")
    verify_affine(saved_A, saved_b, affine, free_indices, len(pivots))
    q, size, order = len(free_indices), len(graphs), 2 ** m
    kernels, matrices = {}, {}
    for block in data["blocks"]:
        m=block.get("m",data["m"]); order=2**m
        key = (block["H"], block["ell"])
        check(key not in matrices, "Repeated kernel block")
        check((m, key[0]) in graphs and key[1] in data["degrees"], "Invalid kernel block index")
        kernel = kernel_data(*key, m, d, alpha, beta)
        C = direct_profile_matrices(kernel, m, size, canonical)
        M = [[[sum(affine[fi][j] * C[fi][a][b] for fi in range(size) if C[fi][a][b])
               for b in range(order)] for a in range(order)] for j in range(q + 1)]
        check(exact_nested(block["C"]) == C, "Saved raw coefficient matrices disagree at %s" % (key,))
        check(exact_nested(block["M"]) == M, "Saved reduced matrices disagree at %s" % (key,))
        kernels[key], matrices[key] = kernel, (C, M)
    required = set()
    for n, H in graphs:
        if n not in anchor_sizes:
            continue
        try:
            kernel_data(H, 0, n, d, alpha, beta)
        except ValueError:
            continue
        required.update((H, ell) for ell in data["degrees"])
    check(set(matrices) == required, "Saved SDP omits an admissible article kernel block")
    z = exact_nested(certificate["z"])
    check(len(z) == size and all(value >= 0 for value in z), "Count multiplier is negative or missing")
    combined = [sum(z[fi] * affine[fi][j] for fi in range(size)) for j in range(q + 1)]
    term_audits = []
    for number, term in enumerate(certificate["terms"]):
        key, f, eta = (term["H"], term["ell"]), exact_nested(term["f"]), Q(term["eta"])
        check(key in matrices, "Certificate uses an undeclared block")
        m=len(kernels[key]["anchor"]); order=2**m
        check(len(f) == order and eta >= 0, "Invalid weight function or negative multiplier")
        raw = enumerate_support_weights(kernels[key], m, f, graphs, arrays, autos)
        C, M = matrices[key]
        direct = [sum(f[a] * f[b] * C[fi][a][b] for a in range(order) for b in range(order))
                  for fi in range(size)]
        check(raw == direct, "Direct article support enumeration disagrees at term %d" % number)
        kappa = [sum(raw[fi] * affine[fi][j] for fi in range(size)) for j in range(q + 1)]
        quadratic = [sum(f[a] * f[b] * M[j][a][b] for a in range(order) for b in range(order))
                     for j in range(q + 1)]
        check(kappa == quadratic, "Reduced quadratic evaluation disagrees")
        combined = [x + eta * y for x, y in zip(combined, kappa)]
        term_audits.append({"H": key[0], "ell": key[1], "eta": str(eta),
                            "f": list(map(str, f)), "kappa": list(map(str, kappa))})
    check(combined == [Q(-1)] + [Q(0)] * q,
          "Certificate cancellation FAILED; exact coefficients = " + str(list(map(str, combined))))
    return {"verified": True, "arithmetic": "Python standard-library fractions.Fraction only",
            "normalization": "N_F = labelled induced embeddings; kernel coefficient divisor = Aut(F)",
            "parameters": parameters, "spectral_embedding": spectral,
            "graph_counts_by_order": [sum(n == t for n, _ in graphs) for t in range(1, p + 1)],
            "independent_counting_equations": len(rows), "counting_rank": len(pivots),
            "free_dimension": q, "checked_kernel_blocks": len(matrices),
            "certificate_terms": term_audits, "nonzero_count_multipliers": sum(value > 0 for value in z),
            "exact_identity_coefficients": list(map(str, combined)),
            "elapsed_seconds": round(time.monotonic() - started, 3)}


def main():
    parser = ArgumentParser(description=__doc__)
    folder = Path(__file__).resolve().parent
    parser.add_argument("data", nargs="?", type=Path, default=folder / "data.json")
    parser.add_argument("certificate", nargs="?", type=Path, default=folder / "certificate.json")
    parser.add_argument("--output", type=Path, default=folder / "verification.json")
    args = parser.parse_args()
    result = verify(json.loads(args.data.read_text()), json.loads(args.certificate.read_text()))
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
