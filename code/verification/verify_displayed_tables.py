#!/usr/bin/env python3
"""Check the four certificates as actually printed in the 1014 TeX source.

Python standard library only.  By default the manuscript is the sibling
../paper/srg_counting_gegenbauer_article_1014.tex; --tex accepts another location.
The input functions and contradiction coefficients are parsed from that
manuscript, not taken from the older certificate JSON files.

The legacy reconstructor supplies an exact five-vertex graph catalogue,
count equations and host-embedding enumeration of scalar kernel sums.
Every row and free coordinate is converted to labelled embedding counts
before the displayed weighted identities are checked.
"""

import argparse
from fractions import Fraction as Q
import importlib.util
import json
from pathlib import Path
import re


HERE = Path(__file__).resolve().parent
PARAMETERS = {
    540: (540, 147, 66, 30),
    550: (550, 162, 75, 36),
    703: (703, 182, 81, 35),
    1344: (1344, 221, 88, 26),
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def table(source, label):
    blocks = re.findall(r"\\begin\{table\}.*?\\end\{table\}", source, re.S)
    found = [b for b in blocks if r"\label{" + label + "}" in b]
    require(len(found) == 1, "Missing or duplicated manuscript table: " + label)
    return found[0]


def rational_cell(cell):
    cell = cell.strip().strip("$").strip()
    fraction = re.fullmatch(r"\\(?:dfrac|tfrac|frac)\{(-?\d+)\}\{(\d+)\}", cell)
    if fraction:
        return Q(int(fraction[1]), int(fraction[2]))
    require(re.fullmatch(r"-?\d+", cell) is not None, "Unrecognized exact cell: " + cell)
    return Q(int(cell))


def manuscript_data(source):
    directions = {v: [] for v in PARAMETERS}
    current = None
    for line in table(source, "tab:direct-functions").splitlines():
        fields = line.strip().removesuffix(r"\\").split("&")
        if len(fields) != 6 or not re.fullmatch(r"[01]{3}", fields[1].strip()):
            continue
        if fields[0].strip():
            current = int(fields[0])
        require(current in directions, "Unknown parameter row in function table")
        bits = fields[1].strip()
        degree, denominator = int(fields[2]), int(fields[3])
        vector = fields[4].strip()
        require(vector.startswith("$(") and vector.endswith(")$"), "Invalid vector cell")
        integers = [int(x.strip()) for x in vector[2:-2].split(",")]
        weight = int(fields[5])
        require(len(integers) == 8 and denominator > 0 and weight > 0, "Invalid function row")
        require(degree in (1, 2), "Unexpected harmonic degree")
        directions[current].append([bits, degree, denominator, integers, weight])
    require({v: len(t) for v, t in directions.items()} == {540: 5, 550: 5, 703: 4, 1344: 5},
            "Expected exactly the 19 displayed functions")

    coefficients = {v: {} for v in PARAMETERS}
    for label, pair in (("tab:direct-coeff-a", (540, 1344)),
                        ("tab:direct-coeff-b", (550, 703))):
        block = table(source, label)
        header = re.search(r"&\$v=(\d+)\$&\$v=(\d+)\$", block)
        require(header is not None and tuple(map(int, header.groups())) == pair,
                "Coefficient table columns do not match their expected cases")
        for line in block.splitlines():
            fields = line.strip().removesuffix(r"\\").split("&")
            if len(fields) != 3:
                continue
            index = re.fullmatch(r"\$b_([0-4])\$", fields[0].strip())
            if index is None:
                continue
            j = int(index[1])
            for v, cell in zip(pair, fields[1:]):
                require(j not in coefficients[v], "Duplicate contradiction coefficient")
                coefficients[v][j] = rational_cell(cell)
    require(all(set(c) == set(range(5)) for c in coefficients.values()),
            "A displayed contradiction coefficient is missing")
    return directions, {v: [c[j] for j in range(5)] for v, c in coefficients.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tex", type=Path,
                        default=HERE.parent / "paper" / "srg_counting_gegenbauer_article_1014.tex")
    args = parser.parse_args()
    require(args.tex.is_file(), "Manuscript not found; use --tex PATH: " + str(args.tex))
    directions, coefficients = manuscript_data(args.tex.read_text(encoding="utf-8"))
    dependency = HERE.parent / "data" / "legacy_ten" / "verify_1344_reconstructed.py"
    require(dependency.is_file(), "Missing portable legacy reconstructor: " + str(dependency))
    spec = importlib.util.spec_from_file_location("displayed_tables_reconstructor", dependency)
    engine = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(engine)

    free_types = [(5, "0011111111"), (5, "0111111011"),
                  (5, "0111111111"), (5, "1111111111")]
    require(engine.GRAPHS[-4:] == free_types, "The four free graph types have changed")
    require([engine.AUT[g] for g in free_types] == [4, 8, 12, 120],
            "Unexpected free-graph automorphism factors")

    results = []
    for v, parameters in PARAMETERS.items():
        _, k, lam, mu = parameters
        r = lam - mu + 3
        dimension = Q(3 * (v - 1) - k, r + 3)
        require(dimension.denominator == 1 and k == 3 * r + mu, "Invalid spectral data")
        dimension = int(dimension)
        alpha = -Q(k + 3, dimension * (r + 3))
        beta = Q(v - k - 3, dimension * (r + 3))
        rows, rhs = engine.count_system(parameters)
        affine = engine.solve_counts(rows, rhs)

        # The imported parametrization uses X_F = N_F/(v Aut F).
        # Convert both dependent and free coordinates to labelled counts.
        affine_labelled = [
            [v * engine.AUT[graph] * a[0]]
            + [engine.AUT[graph] * a[j + 1] / engine.AUT[free_types[j]]
               for j in range(4)]
            for graph, a in zip(engine.GRAPHS, affine)
        ]
        total = [Q(0)] * 5
        for term in directions[v]:
            row, diagnostic = engine.kernel_row(parameters, dimension, alpha, beta, term)
            # The imported kernel row is S/(v Aut H) in X coordinates.
            # This gives the scalar coefficients c_{H,ell,f}(F) in N.
            aut_h = engine.AUT[3, term[0]]
            labelled_row = [q * aut_h / engine.AUT[graph]
                            for q, graph in zip(row, engine.GRAPHS)]
            multiplier = Q(term[-1], aut_h)
            require(all(Q(p) > 0 for p in diagnostic["positive_elimination_pivots"]),
                    "A displayed reference Gram matrix is not positive definite")
            for j in range(5):
                total[j] += multiplier * sum(q * a[j] for q, a in zip(labelled_row, affine_labelled))

        expected = [-b for b in coefficients[v]]
        require(total == expected, "Printed labelled-count identity fails for v=" + str(v))
        require(all(b > 0 for b in coefficients[v]), "Contradiction coefficient is not positive")
        # Moving all free counts to the left and dividing by b_0 gives -1.
        normalized = [total[0] / coefficients[v][0]] + [
            (total[j] + coefficients[v][j]) / coefficients[v][0] for j in range(1, 5)
        ]
        require(normalized == [-1, 0, 0, 0, 0], "Normalized dual cancellation failed")
        results.append({"parameters": list(parameters), "functions": len(directions[v]),
                        "labelled_identity_verified": True,
                        "positive_b": list(map(str, coefficients[v])),
                        "normalized_dual_coefficients": list(map(str, normalized))})
    print(json.dumps({"verified": True, "input": args.tex.name,
                      "arithmetic": "Python fractions.Fraction; exact",
                      "source_of_functions_and_coefficients": "the displayed manuscript tables",
                      "cases": results}, indent=2))


if __name__ == "__main__":
    main()
