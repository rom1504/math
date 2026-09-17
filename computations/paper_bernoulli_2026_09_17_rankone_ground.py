"""Exact Gram-moment replay for the whole Boolean rank-one code obstruction."""

from fractions import Fraction
import itertools
import json

import numpy as np

from paper_discrepancy_2026_09_17_involution_ground import matchings, parity


def main():
    p = 8
    n = p * p
    coordinates = list(itertools.product(range(p), repeat=2))
    rows = []
    for pairs in matchings(tuple(range(p))):
        permutation = [0] * p
        for left, right in pairs:
            permutation[left], permutation[right] = right, left
        for bits in itertools.product((0, 1), repeat=p // 2):
            cycle_bits = [0] * p
            for (left, right), bit in zip(pairs, bits):
                cycle_bits[left] = cycle_bits[right] = bit
            rows.append([(-1) ** (parity(u & permutation[v]) + cycle_bits[v]
                                  + parity(u & v) + parity(u) + parity(v))
                         for u, v in coordinates])
    law = np.asarray(rows, dtype=np.int64)
    atoms = len(law)
    cp = Fraction(2 * p - 1, (p - 1) * (p - 3))
    rng = np.random.default_rng(2026091705)
    reports = []
    for k in (2, 3, 4, 8, 16):
        m = n // k
        n0 = k * m
        switching = rng.choice((-1, 1), size=n)
        coordinate_order = rng.permutation(n)
        matrices = (law[:, coordinate_order] * switching)[:, :n0].reshape(atoms, k, m)
        grams = np.einsum("aij,akj->aik", matrices, matrices)
        centered = grams - m * np.eye(k, dtype=np.int64)
        assert np.all(np.diagonal(centered, axis1=1, axis2=2) == 0)
        squared_frobenius = np.sum(centered ** 2, axis=(1, 2))
        mean = Fraction(int(squared_frobenius.sum()), atoms * n0 ** 2)
        bound = (1 - Fraction(1, k)) * cp + Fraction((k - 1) * p, 2 * n0) * (1 - cp)
        assert mean <= bound
        largest_gram_eigenvalues = np.linalg.eigvalsh(grams)[:, -1]
        assert np.all(largest_gram_eigenvalues <= m + np.sqrt(squared_frobenius) + 1e-10)
        report = {"k": k, "m": m, "leftovers": n - n0,
                  "exact_mean_normalized_frobenius_squared": str(mean),
                  "proved_upper_bound": str(bound)}
        if k <= 4:
            # For each c, maximize over v exactly by taking the sum of
            # absolute column fields; then maximize over all c.
            left_signs = np.asarray(list(itertools.product((-1, 1), repeat=k)))
            responses = np.einsum("bk,akj->abj", left_signs, matrices)
            maxima = np.abs(responses).sum(axis=2).max(axis=1)
            assert np.all(maxima ** 2 <= n0 * largest_gram_eigenvalues + 1e-8)
            report["smallest_exact_max_correlation_to_rankone_code"] = int(maxima.min())
        reports.append(report)
    print(json.dumps({"status": "PASS", "ground_law_atoms": atoms,
                      "n": n, "arbitrary_switched_grouping_checks": reports,
                      "scope": "Finite replay of moment proof, not asymptotic simulation."}, indent=2))


if __name__ == "__main__":
    main()
