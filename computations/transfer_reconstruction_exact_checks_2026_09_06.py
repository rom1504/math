"""Independent finite exact checks for the standalone upper-bound reconstruction.

These are regression checks of identities, not substitutes for the asymptotic
proof. No floating arithmetic enters the reported inequalities. The optional
full certificate runs the existing exact continuous-envelope verifier too.
"""

import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
import json
from math import factorial
import random


def permanent(matrix):
    order = len(matrix)
    return sum(
        product_entries(matrix[i][j] for i, j in enumerate(perm))
        for perm in permutations(range(order))
    )


def product_entries(values):
    answer = 1
    for value in values:
        answer *= value
    return answer


def gram(vectors):
    return [[sum(a * b for a, b in zip(x, y)) for y in vectors] for x in vectors]


def paley_twelve():
    q = 11
    residues = {i * i % q for i in range(1, q)}
    chi = lambda z: 0 if z % q == 0 else (1 if z % q in residues else -1)
    rows = [[1] * (q + 1)]
    rows.extend(
        [[-1] + [int(x == y) + chi(x - y) for y in range(q)] for x in range(q)]
    )
    assert all(abs(entry) == 1 for row in rows for entry in row)
    assert gram(rows) == [[12 * int(i == j) for j in range(12)] for i in range(12)]
    return rows


def check_pair_tables():
    # Six positions: two zero coordinates and four coordinates of magnitude one.
    order = 6
    half = order // 2
    table_counts = Counter()
    words = []
    for support in combinations(range(order), 4):
        for signs in product((-1, 1), repeat=4):
            word = [0] * order
            for index, sign in zip(support, signs):
                word[index] = sign
            words.append(tuple(word))
            pairs = Counter(zip(word[:half], word[half:]))
            table_counts[tuple(sorted(pairs.items()))] += 1
    expected = factorial(order) * 2**4 // (factorial(2) * factorial(4))
    assert len(words) == expected == 240
    source = {0: F(1, 3), -1: F(1, 3), 1: F(1, 3)}
    for table, count in table_counts.items():
        assert count == factorial(half) // product_entries(factorial(n) for _, n in table)
        original = {pair: F(n, half) for pair, n in table}
        averaged = Counter()
        for (a, b), mass in original.items():
            for aa, bb in ((a, b), (-a, -b), (b, a), (-b, -a)):
                averaged[(aa, bb)] += mass / 4
        for coordinate in (0, 1):
            marginal = Counter()
            for pair, mass in averaged.items():
                marginal[pair[coordinate]] += mass
            assert dict(marginal) == source
        for sign in (-1, 1):
            before, after = Counter(), Counter()
            for (a, b), mass in original.items():
                before[abs(a + sign * b)] += mass
            for (a, b), mass in averaged.items():
                after[abs(a + sign * b)] += mass
            assert before == after
    return {"words": expected, "pair_tables": len(table_counts)}


def check_permanent_projection(rng):
    tested = 0
    for order in range(1, 7):
        for _ in range(4):
            aa = gram([[rng.randrange(-2, 3) for _ in range(3)] for _ in range(order)])
            bb = gram([[rng.randrange(-2, 3) for _ in range(3)] for _ in range(order)])
            hadamard_product = [[aa[i][j] * bb[i][j] for j in range(order)] for i in range(order)]
            assert permanent(hadamard_product) * factorial(order) >= permanent(aa) * permanent(bb)
            if order > 1:
                split = order // 2
                left = [row[:split] for row in aa[:split]]
                right = [row[split:] for row in aa[split:]]
                assert permanent(aa) * factorial(split) * factorial(order - split) <= (
                    permanent(left) * permanent(right) * factorial(order)
                )
            tested += 1
    return tested


def check_graph_contraction():
    # Four degree-three vertex multisets, including repeated coordinate values.
    # The kernel 1+a^2 b^2 is an exact PSD nonnegative Gram kernel.
    lists = ((0, 1, 2), (1, 1, 2), (0, 0, 3), (1, 2, 3))
    vertices = range(len(lists))
    neighbors = {i: [j for j in vertices if i != j] for i in vertices}
    arrangements = [list(permutations(row)) for row in lists]
    kernel = lambda a, b: 1 + a * a * b * b
    total = 0
    for rows in product(*arrangements):
        entries = {
            (i, j): rows[i][neighbors[i].index(j)]
            for i in vertices for j in neighbors[i]
        }
        total += product_entries(
            kernel(entries[(i, j)], entries[(j, i)])
            for i, j in combinations(vertices, 2)
        )
    average = F(total, 6**4)
    norm_square_product = product_entries(
        F(permanent([[kernel(a, b) for b in row] for a in row]), factorial(len(row)))
        for row in lists
    )
    assert average * average <= norm_square_product
    return {"moment": str(average), "squared_norm_product": str(norm_square_product)}


def check_weave_defects(rng):
    order = 4
    retained = 3
    h = [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]]
    selectors = [rng.sample(range(order), retained) for _ in range(order)]
    bases = []
    for _ in range(order):
        perm = rng.sample(range(order), order)
        signs = [rng.choice((-1, 1)) for _ in range(order)]
        bases.append([[row[perm[j]] * signs[j] for j in range(order)] for row in h])
    sign_matrix = [[0] * order for _ in range(order)]
    for i in range(order):
        for j in range(i, order):
            sign_matrix[i][j] = sign_matrix[j][i] = rng.choice((-1, 1))
    matrix = [
        [
            sign_matrix[i][j] * bases[i][selectors[i][a]][j] * bases[j][selectors[j][b]][i]
            for j in range(order) for b in range(retained)
        ]
        for i in range(order) for a in range(retained)
    ]
    assert all(matrix[i][j] == matrix[j][i] for i in range(12) for j in range(12))
    diagonal = sum(matrix[i][i] for i in range(12))
    full_cap = hollow_cap = 0
    for spins in product((-1, 1), repeat=order * retained):
        spectra = [
            [sum(bases[i][selectors[i][a]][j] * spins[i * retained + a] for a in range(retained))
             for j in range(order)]
            for i in range(order)
        ]
        assert sum(v * v for row in spectra for v in row) == order * order * retained
        energy = sum(matrix[i][j] * spins[i] * spins[j] for i in range(12) for j in range(12))
        assert energy == sum(
            sign_matrix[i][j] * spectra[i][j] * spectra[j][i]
            for i in range(order) for j in range(order)
        )
        for sigma in (-1, 1):
            defect = sum(
                (spectra[i][j] - sigma * sign_matrix[i][j] * spectra[j][i])**2
                for i in range(order) for j in range(order)
            )
            assert defect == 2 * (order * order * retained - sigma * energy)
        full_cap = max(full_cap, abs(energy))
        hollow_cap = max(hollow_cap, abs(energy - diagonal))
    assert hollow_cap <= full_cap + 12
    return {"all_spins": 2**12, "full_matrix_cap": full_cap, "hollow_matrix_cap": hollow_cap}


def independent_log_interval(x, terms=90):
    # The reference interval uses no binary range reduction and is practical
    # here because the test arguments stay in [1/3,3]. Signed atanh remainder.
    x = F(x)
    z = (x - 1) / (x + 1)
    partial = 2 * sum(z ** (2 * j + 1) / (2 * j + 1) for j in range(terms))
    tail = 2 * abs(z) ** (2 * terms + 1) / ((2 * terms + 1) * (1 - z * z))
    return (partial, partial + tail) if z >= 0 else (partial - tail, partial)


def check_interval_primitives():
    from continued_feedback_ternary_latent_interval_certificate_2026_09_06 import (
        log_interval, sqrt_interval,
    )
    checked = 0
    for denominator in range(2, 14):
        for numerator in range(denominator // 2 + 1, 2 * denominator + 1):
            value = F(numerator, denominator)
            lo, hi = log_interval(value)
            reflo, refhi = independent_log_interval(value)
            assert lo <= reflo <= refhi <= hi
            slo, shi = sqrt_interval(value)
            assert slo * slo <= value <= shi * shi
            checked += 1
    p = F(31, 32)
    a = F(91470529542342299, 20460000000000000000)
    _, root_upper = sqrt_interval(p)
    cap = F(1, 2) - a / (8 * root_upper)
    assert cap == F(80459630021641337701, 161102201102367360000)
    assert cap < F(499432220485404, 10**15)
    return {"interval_arguments": checked, "outward_cap": str(cap)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-certificate", action="store_true")
    args = parser.parse_args()
    rng = random.Random(20260906)
    paley_twelve()
    result = {
        "paley_twelve_exact": True,
        "pair_tables": check_pair_tables(),
        "permanent_projection_cases": check_permanent_projection(rng),
        "graph_contraction": check_graph_contraction(),
        "weave": check_weave_defects(rng),
        "interval_primitives": check_interval_primitives(),
    }
    if args.full_certificate:
        from continued_feedback_conditional_variance_exact_certificate_2026_09_06 import certificate
        result["full_certificate"] = certificate()
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
