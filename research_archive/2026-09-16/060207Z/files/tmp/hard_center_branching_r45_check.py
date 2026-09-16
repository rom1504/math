#!/usr/bin/env python3
"""Exact finite audit of the Wave 45 hard common-core branching criterion."""

from __future__ import annotations

import itertools
import math
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9
from high_replica_pressure_r44_explore import geometry


def down_up_kernel(selectors: list[frozenset[int]], n: int, m: int, ell: int) -> np.ndarray:
    """K(S,T): choose an ell-subset of S, then an m-superset of it."""
    denominator = math.comb(m, ell) * math.comb(n - ell, m - ell)
    return np.asarray(
        [
            [Fraction(math.comb(len(s & t), ell), denominator) for t in selectors]
            for s in selectors
        ],
        dtype=object,
    )


def audit(name: str, a: np.ndarray, m: int) -> dict:
    z, selector_tuples, qn, energy, row, deficit, h = geometry(a, m)
    n = len(a)
    selectors = list(map(frozenset, selector_tuples))
    # The old C2 class is a stronger finite row qualification than the
    # asymptotic project cap.  It happens to contain every center here.
    eligible = row <= 2 * n * (n - 1)
    favorable = (deficit <= h + 1e-12) & eligible[:, None]
    assert np.all(eligible)
    center_count, selector_count = favorable.shape
    total_incidence = int(np.sum(favorable))
    assert total_incidence > 0
    degrees = np.sum(favorable, axis=1)
    max_degree = Fraction(int(np.max(degrees)), selector_count)
    first_moment = Fraction(total_incidence, center_count * selector_count)
    second_ratio = Fraction(int(np.sum(degrees * degrees)), selector_count * total_incidence)

    rows = []
    for ell in range(m):
        kernel = down_up_kernel(selectors, n, m, ell)
        # Exact stochasticity and symmetry.
        assert all(sum(kernel[i, :]) == 1 for i in range(selector_count))
        assert np.array_equal(kernel, kernel.T)

        # The common-core square identity, computed independently through
        # f K f and through conditional favorable densities b_z(R)^2.
        collision_numerator = Fraction(0)
        for f in favorable:
            ids = np.flatnonzero(f)
            for i in ids:
                for j in ids:
                    collision_numerator += kernel[i, j]
        rho = collision_numerator / total_incidence

        cores = list(itertools.combinations(range(n), ell))
        core_square = Fraction(0)
        for f in favorable:
            for core_tuple in cores:
                core = frozenset(core_tuple)
                supersets = [j for j, s in enumerate(selectors) if core <= s]
                b = Fraction(sum(bool(f[j]) for j in supersets), len(supersets))
                core_square += b * b
        core_square /= center_count * len(cores)
        incidence_mean = Fraction(total_incidence, center_count * selector_count)
        assert core_square / incidence_mean == rho

        lambda_one = Fraction(ell * (n - m), m * (n - ell)) if ell else Fraction(0)
        # Independent numerical spectrum audit of K=D*D.
        numeric = np.asarray([[float(v) for v in row_] for row_ in kernel], dtype=float)
        eigenvalues = np.linalg.eigvalsh(numeric)
        assert eigenvalues[0] >= -2e-12
        assert abs(eigenvalues[-1] - 1.0) < 2e-12
        assert abs(eigenvalues[-2] - float(lambda_one)) < 2e-12

        # Check max a >= rho^t-lambda^t for several powers.
        for power in range(1, 7):
            assert max_degree >= rho**power - lambda_one**power
        rows.append((ell, rho, lambda_one, rho - lambda_one))

    empty_pairs = sum(
        not np.any(favorable[:, i] & favorable[:, j])
        for i, j in itertools.combinations(range(selector_count), 2)
    )
    empty_triples = sum(
        not np.any(favorable[:, i] & favorable[:, j] & favorable[:, k])
        for i, j, k in itertools.combinations(range(selector_count), 3)
    )
    return {
        "name": name,
        "n": n,
        "m": m,
        "centers": center_count,
        "selectors": selector_count,
        "max_degree": max_degree,
        "first_moment": first_moment,
        "second_ratio": second_ratio,
        "rows": rows,
        "empty_pairs": (empty_pairs, math.comb(selector_count, 2)),
        "empty_triples": (empty_triples, math.comb(selector_count, 3)),
        "max_row": int(np.max(row)),
    }


def main() -> None:
    audits = [audit("A6", A6, 5), audit("A8", A8, 6), audit("A9", A9, 7)]
    for result in audits:
        print(
            f"{result['name']}: n={result['n']} m={result['m']} "
            f"centers={result['centers']} selectors={result['selectors']} "
            f"maxR={result['max_row']} max-a={result['max_degree']} "
            f"J1={result['first_moment']} J2/J1={result['second_ratio']} "
            f"empty-pairs={result['empty_pairs'][0]}/{result['empty_pairs'][1]} "
            f"empty-triples={result['empty_triples'][0]}/{result['empty_triples'][1]}"
        )
        for ell, rho, eigenvalue, gap in result["rows"]:
            print(f"  ell={ell}: rho={rho} lambda={eigenvalue} gap={gap}")

    assert audits[1]["rows"][-1][3] == Fraction(-11, 225)
    assert audits[2]["rows"][-2][3] == Fraction(-1, 21168)
    assert audits[2]["rows"][-1][3] == Fraction(-55, 882)
    print("PASS: exact hard collision, common-core square, PSD spectrum, and branching bounds")


if __name__ == "__main__":
    main()
