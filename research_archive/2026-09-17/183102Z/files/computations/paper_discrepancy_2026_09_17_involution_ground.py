"""Exact finite checks for the random-involution near-ground law."""

from __future__ import annotations

from fractions import Fraction
import itertools
import json

import numpy as np


def matchings(labels):
    if not labels:
        yield ()
        return
    first = labels[0]
    for i in range(1, len(labels)):
        second = labels[i]
        remaining = labels[1:i]+labels[i+1:]
        for rest in matchings(remaining):
            yield ((first, second),)+rest


def parity(value):
    return bin(int(value)).count("1") % 2


def main():
    p = 8
    a = 3
    n = p*p
    coordinates = list(itertools.product(range(p), repeat=2))
    rows = []
    for pairs in matchings(tuple(range(p))):
        permutation = [0]*p
        for left, right in pairs:
            permutation[left], permutation[right] = right, left
        for bits in itertools.product((0, 1), repeat=p//2):
            g = [0]*p
            for (left, right), bit in zip(pairs, bits):
                g[left] = g[right] = bit
            rows.append([(-1)**(parity(u & permutation[v])+g[v])
                         for u, v in coordinates])
    x = np.asarray(rows, dtype=np.int64)
    count = len(x)
    assert count == 105*16
    fourier = np.asarray([[(-1)**(parity(s & v)+parity(t & u))
                           for u, v in coordinates] for s, t in coordinates], dtype=np.int64)
    assert np.array_equal(x @ fourier, p*x)
    diagonal = np.asarray([(-1)**(parity(u & v)+parity(u)+parity(v))
                           for u, v in coordinates], dtype=np.int64)
    h = (-1)**a*diagonal[:, None]*fourier*diagonal[None, :]
    y = x*diagonal
    assert np.array_equal(y @ h, (-1)**a*p*y)
    assert np.all(np.diag(h) == (-1)**a)
    # Pair moments, evaluated over the complete finite probability law.
    gram = x.T @ x
    np.fill_diagonal(gram, 0)
    assert int(np.abs(gram).max())*(p-1) <= count
    bound = Fraction(2*p-1, (p-1)*(p-3))
    rng = np.random.default_rng(202609174)
    tests = [tuple(sorted(rng.choice(n, 4, replace=False))) for _ in range(10000)]
    for v in range(p):
        tests.extend(tuple(u*p+v for u in values)
                     for values in itertools.combinations(range(p), 4))
    exceptional_count = 0
    max_nonexception = 0
    for indices in tests:
        values = [coordinates[i] for i in indices]
        xor_u = values[0][0]^values[1][0]^values[2][0]^values[3][0]
        exceptional = len({v for u, v in values}) == 1 and xor_u == 0
        numerator = abs(int(np.prod(x[:, indices], axis=1).sum()))
        if exceptional:
            assert numerator == count
            exceptional_count += 1
        else:
            assert numerator*bound.denominator <= count*bound.numerator
            max_nonexception = max(max_nonexception, numerator)
    for _ in range(100):
        order = rng.permutation(n)
        m = int(rng.integers(1, n//2+1))
        edges = order[:2*m].reshape(m, 2)
        signs = rng.choice([-1, 1], size=m)
        statistic = (y[:, edges[:, 0]]*y[:, edges[:, 1]]*signs).sum(axis=1)
        second_moment_sum = int(statistic @ statistic)
        upper = bound*m*m+(1-bound)*m*p/2
        assert second_moment_sum*upper.denominator <= count*upper.numerator
    mode_cases = []
    for k in (2, 4, 8, 16):
        order = rng.permutation(n)
        switches = rng.choice([-1, 1], size=n)
        grouped = (y[:, order]*switches).reshape(count, k, n//k)
        hadamard = np.asarray([[(-1)**parity(b & d) for d in range(k)]
                               for b in range(k)], dtype=np.int64)
        mode_sums = np.einsum("ab,rbj->raj", hadamard, grouped)
        t = (mode_sums*mode_sums).sum(axis=2)
        statistics = []
        for difference in range(1, k):
            edges = [(b, b ^ difference) for b in range(k) if b < (b ^ difference)]
            statistic = sum((grouped[:, left]*grouped[:, right]).sum(axis=1)
                            for left, right in edges)
            statistics.append(statistic)
        u = np.asarray(statistics).T
        assert np.array_equal((t*t).sum(axis=1), k*n*n+4*k*(u*u).sum(axis=1))
        defect = Fraction(int((t*t).sum()), count*k*n*n)-1
        upper = (k-1)*(bound+(1-bound)/p)
        assert 0 <= defect <= upper
        mode_cases.append({"k": k, "exact_expected_defect": str(defect),
                           "proved_defect_upper_bound": str(upper)})
    print(json.dumps({"status": "all exact finite-law checks passed", "p": p, "n": n,
                      "law_atoms_with_multiplicity": count, "fourth_moment_tests": len(tests),
                      "tested_exceptional_quadruples": exceptional_count,
                      "max_tested_nonexception_fourth_moment": str(Fraction(max_nonexception, count)),
                      "uniform_proved_bound": str(bound), "signed_matchings_tested": 100,
                      "switched_permuted_mode_checks": mode_cases}, indent=2))


if __name__ == "__main__":
    main()
