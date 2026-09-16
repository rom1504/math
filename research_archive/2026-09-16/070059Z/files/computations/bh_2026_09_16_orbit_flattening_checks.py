"""Exact finite checks for the permutation-flattening appendix."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
import json

from bh_2026_09_16_power_checks import (
    bit_count,
    energies,
    fwht,
    stored_witness,
)


def main():
    n = 5
    edges = list(combinations(range(n), 2))
    signs = {
        edge: 1 if min((edge[1] - edge[0]) % n,
                       (edge[0] - edge[1]) % n) == 1 else -1
        for edge in edges
    }
    assert all(sum(a for e, a in signs.items() if i in e) == 0
               for i in range(n))
    orbit = [
        [signs[tuple(sorted((p[i], p[j])))] for i, j in edges]
        for p in permutations(range(n))
    ]
    for i, e in enumerate(edges):
        for j, f in enumerate(edges):
            actual = Fraction(sum(v[i] * v[j] for v in orbit), len(orbit))
            common = len(set(e).intersection(f))
            expected = (Fraction(1) if common == 2 else
                        Fraction(-1, n - 2) if common == 1 else
                        Fraction(2, (n - 2) * (n - 3)))
            assert actual == expected

    orbit_values = []
    for x in product((-1, 1), repeat=n):
        actual = Fraction(sum(
            sum(v[j] * x[e[0]] * x[e[1]] for j, e in enumerate(edges)) ** 2
            for v in orbit), len(orbit))
        s = sum(x)
        projected = (Fraction(n * (n - 1), 2) - s * s
                     - Fraction(n, n - 2)
                     + Fraction((s * s - n) ** 2, 2 * (n - 2) * (n - 1)))
        assert actual == Fraction(n - 1, n - 3) * projected
        orbit_values.append(actual)
    hilbert_cap_squared = max(orbit_values)
    assert hilbert_cap_squared == Fraction((n - 1) ** 2 * (n + 1),
                                          2 * (n - 2))

    matrix6, source6, _ = stored_witness(6)
    h6 = energies(matrix6)
    for x, h in enumerate(h6):
        parity = -1 if bit_count(x) % 2 else 1
        assert h * h == 15 + 2 * parity * h

    matrix7, source7, _ = stored_witness(7)
    h7 = energies(matrix7)
    raw7 = fwht([h * h for h in h7])
    assert all(v % 128 == 0 for v in raw7)
    c7 = [v // 128 for v in raw7]
    rows = {}
    for r, expected in ((2, {2: 18, 6: 3}), (4, {2: 30, 6: 5})):
        level = [abs(v) for s, v in enumerate(c7) if bit_count(s) == r]
        histogram = dict(Counter(level))
        assert histogram == expected
        eta = Fraction(sum(level) ** 2,
                       len(level) * sum(v * v for v in level))
        assert eta == Fraction(27, 35)
        rows[r] = {"absolute_histogram": histogram, "eta": str(eta)}

    # Independent-sign auxiliary scalarization: multiply by sqrt(D) so
    # all arithmetic remains integral. Its exact supremum is D.
    aux_edges = list(combinations(range(4), 2))
    aux_cap = max(
        abs(sum(y[j] * x[a] * x[b] for j, (a, b) in enumerate(aux_edges)))
        for x in product((-1, 1), repeat=4)
        for y in product((-1, 1), repeat=len(aux_edges))
    )
    assert aux_cap == len(aux_edges)
    print(json.dumps({
        "status": "PASS",
        "row_balanced_gram": {"n": n, "orbit_size": len(orbit),
                              "hilbert_cap_squared": str(hilbert_cap_squared)},
        "H6_exception_source": source6,
        "H7_nonflat_source": source7,
        "H7_square_levels": rows,
        "auxiliary_scaled_cap": aux_cap,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
