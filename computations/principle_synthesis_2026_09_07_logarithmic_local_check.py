"""Exact finite checks of the logarithmic local-law Fourier identity.

Run with .venv/bin/python. No solver or ignored input is used. Counts and
Fourier identity checks are rational; logarithmic KL reporting is numerical.
"""

from collections import Counter
from fractions import Fraction
from itertools import product
import json
import math


def popcount(value):
    # The repository environment is Python 3.9, before int.bit_count().
    return bin(value).count("1")


def edge_code(matrix, vertices, switches):
    code = 0
    bit = 0
    for i in range(len(vertices)):
        for j in range(i + 1, len(vertices)):
            value = switches[i] * switches[j] * matrix[vertices[i]][vertices[j]]
            code |= (value == 1) << bit
            bit += 1
    return code


def full_counts(matrix, k):
    n = len(matrix)
    result = Counter()
    for vertices in product(range(n), repeat=k):
        for switches in product((-1, 1), repeat=k):
            result[edge_code(matrix, vertices, switches)] += 1
    return result


def exact_fourth_trace(matrix):
    n = len(matrix)
    gram = [[sum(matrix[i][v] * matrix[j][v] for v in range(n))
             for j in range(n)] for i in range(n)]
    return sum(x * x for row in gram for x in row), gram


def row_identity(matrix, t):
    n = len(matrix)
    total = Fraction(0)
    for old_vertices in product(range(n), repeat=t):
        counts = Counter()
        for fresh in range(n):
            for fresh_switch in (-1, 1):
                pattern = tuple(fresh_switch * matrix[fresh][v] for v in old_vertices)
                counts[pattern] += 1
        chi2 = Fraction((2 ** t) * sum(c * c for c in counts.values()), (2 * n) ** 2) - 1
        total += chi2
    total /= n ** t
    fourth, gram = exact_fourth_trace(matrix)
    formula = sum(
        ((1 + Fraction(g, n)) ** t + (1 - Fraction(g, n)) ** t) / 2 - 1
        for row in gram for g in row
    ) / n ** 2
    assert total == formula, (total, formula)
    bound = Fraction((2 ** (t - 1) - 1) * fourth, n ** 4)
    assert total <= bound
    return {"t": t, "mean_chi2": str(total), "bound": str(bound)}


def run_case(name, matrix, k):
    n = len(matrix)
    counts = full_counts(matrix, k)
    size = (2 * n) ** k
    alphabet_size = 2 ** (k * (k - 1) // 2)
    assert sum(counts.values()) == size
    divergence = sum((c / size) * math.log(c * alphabet_size / size) for c in counts.values())
    tv = 0.5 * (sum(abs(c / size - 1 / alphabet_size) for c in counts.values())
                + (alphabet_size - len(counts)) / alphabet_size)
    fourth, _ = exact_fourth_trace(matrix)
    bound = Fraction((2 ** (k - 1) - k) * fourth, n ** 4)
    assert divergence <= float(bound) + 1e-12
    # Every nonempty one-edge or two-distinct-edge moment vanishes exactly.
    edges = k * (k - 1) // 2
    for mask in range(1, 1 << edges):
        if popcount(mask) <= 2:
            assert sum(c * (-1 if popcount(code & mask) % 2 else 1)
                       for code, c in counts.items()) == 0
    return {"name": name, "n": n, "k": k, "support": len(counts),
            "KL": divergence, "KL_upper_exact": str(bound), "TV": tv,
            "rows": [row_identity(matrix, t) for t in range(1, k)]}


def main():
    hadamard4 = [[(-1) ** (popcount(i & j) % 2) for j in range(4)]
                 for i in range(4)]
    cycle5 = [[1 if i == j or (i - j) % 5 in (1, 4) else -1
               for j in range(5)] for i in range(5)]
    arbitrary3 = [[1, 1, -1], [1, 1, 1], [-1, 1, 1]]
    output = [run_case("symmetric_H4", hadamard4, 4),
              run_case("pentagon_plus_diagonal", cycle5, 4),
              run_case("arbitrary_order3", arbitrary3, 4)]
    print(json.dumps({"status": "PASS", "cases": output}, indent=2))


if __name__ == "__main__":
    main()
