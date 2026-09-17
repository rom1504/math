"""Exact negative-eigenpole checks for the coordinate-matching obstruction.

The asymptotic theorem is proved in the accompanying artifact; this
enumerates the complete p=8 query law and tests moment normalizations.
"""

from fractions import Fraction as F
from itertools import combinations, product
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def main():
    spec = importlib.util.spec_from_file_location(
        "involution", ROOT / "computations/paper_discrepancy_2026_09_17_involution_ground.py")
    source = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    p = 8
    n = p * p
    coordinates = list(product(range(p), repeat=2))
    rows = []
    for pairs in source.matchings(tuple(range(p))):
        pi = [0] * p
        for left, right in pairs:
            pi[left], pi[right] = right, left
        for signs in product((-1, 1), repeat=p//2):
            z = [0] * p
            for (left, right), sign in zip(pairs, signs):
                z[left], z[right] = sign, -sign
            rows.append([(-1)**source.parity(u & pi[v]) * z[v]
                         for u, v in coordinates])
    x = np.asarray(rows, dtype=np.int64)
    count = len(x)
    assert count == 1680
    fourier = np.asarray([[(-1)**(source.parity(u & t)+source.parity(v & s))
                           for s, t in coordinates] for u, v in coordinates], dtype=np.int64)
    assert np.array_equal(fourier @ fourier, n * np.eye(n, dtype=np.int64))
    assert np.array_equal(x @ fourier, -p * x)
    covariance_numerator = x.T @ x
    assert np.array_equal((p-1) * covariance_numerator,
                          count * (p*np.eye(n, dtype=np.int64)-fourier))
    cap = n * (p+1) // 2
    energy = ((x @ fourier) * x).sum(axis=1)//2 - n//2
    assert np.all(energy == -cap)
    cp = F(2*p-1, (p-1)*(p-3))
    rng = np.random.default_rng(202609172125)
    tests = [tuple(sorted(rng.choice(n, 4, replace=False))) for _ in range(10000)]
    for v in range(p):
        tests.extend(tuple(u*p+v for u in indices) for indices in combinations(range(p), 4))
    exceptions = 0
    for indices in tests:
        values = [coordinates[i] for i in indices]
        xor_u = values[0][0] ^ values[1][0] ^ values[2][0] ^ values[3][0]
        exceptional = len({v for u, v in values}) == 1 and xor_u == 0
        numerator = abs(int(np.prod(x[:, indices], axis=1).sum()))
        if exceptional:
            assert numerator == count
            exceptions += 1
        else:
            assert numerator * cp.denominator <= count * cp.numerator
    m = n // 2
    abs_means = [F(sum(math.comb(k, j)*abs(k-2*j) for j in range(k+1)), 2**k)
                 for k in range(m+1)]
    costs = [abs_means[k]+abs_means[m-k] for k in range(m+1)]
    den = math.lcm(*(v.denominator for v in costs))
    lookup = np.asarray([int(v*den) for v in costs], dtype=np.int64)
    mean_responses = []
    for _ in range(128):
        order = rng.permutation(n)
        pairs = order.reshape(m, 2)
        signs = rng.choice((-1, 1), size=m)
        statistic = (x[:, pairs[:, 0]] * x[:, pairs[:, 1]] * signs).sum(axis=1)
        actual_second = F(int(statistic @ statistic), count)
        bound = cp*m*m+(1-cp)*m*p/2
        assert actual_second <= bound
        active = (m+statistic)//2
        mean_responses.append(F(int(sum(lookup[active])), den*count))
    result = {
        "status": "PASS exact negative-eigenpole ground-law checks",
        "p": p, "n": n, "query_law_atoms_with_multiplicity": count,
        "exact_cap": cap, "all_atoms_are_exact_absolute_ground_words": True,
        "exact_covariance_identity": "(pI-F)/(p-1)",
        "fourth_moment_checks": len(tests), "tested_exceptional_quadruples": exceptions,
        "signed_coordinate_matchings_tested": len(mean_responses),
        "minimum_tested_exact_average_response": str(min(mean_responses)),
        "maximum_tested_exact_average_response": str(max(mean_responses)),
        "scope": "Finite normalization test; the uniform all-matching asymptotic theorem is analytical.",
    }
    output = ROOT / "tmp/paper_portfolio_2026_09_17/bernoulli/coordinate_matching_hadamard.json"
    output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
