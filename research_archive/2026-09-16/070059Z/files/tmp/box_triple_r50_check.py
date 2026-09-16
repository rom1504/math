#!/usr/bin/env python3
"""Exact Wave 50 check: box mass does not force coarea at the same row cap.

The matrix is written explicitly, so the obstruction is reproducible without
calling the MILP sampler that originally found it.  All signings, all
six-subsets, and all relevant completions are exhausted.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter
from fractions import Fraction

import numpy as np


A = np.asarray(
    [
        [0, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [1, 0, -1, 1, -1, 1, 1, 1, -1, -1],
        [1, -1, 0, 1, 1, 1, -1, 1, 1, 1],
        [1, 1, 1, 0, 1, 1, -1, -1, -1, -1],
        [1, -1, 1, 1, 0, -1, 1, -1, -1, 1],
        [1, 1, 1, 1, -1, 0, 1, 1, 1, -1],
        [1, 1, -1, -1, 1, 1, 0, -1, 1, -1],
        [1, 1, 1, -1, -1, 1, -1, 0, -1, 1],
        [1, -1, 1, -1, -1, 1, 1, -1, 0, 1],
        [1, -1, 1, -1, 1, -1, -1, 1, 1, 0],
    ],
    dtype=np.int64,
)


def projective_spins(n: int) -> list[np.ndarray]:
    return [np.asarray((1,) + tail, dtype=np.int64)
            for tail in itertools.product((-1, 1), repeat=n - 1)]


def energy(a: np.ndarray, x: np.ndarray) -> int:
    return int(x @ a @ x)


def row_square(a: np.ndarray, x: np.ndarray) -> int:
    return int(np.sum((a @ x) ** 2))


def falling(x: int, j: int) -> int:
    ans = 1
    for i in range(j):
        ans *= x - i
    return ans


def build_ground_relation(m: int):
    n = len(A)
    centers = projective_spins(n)
    selectors = [tuple(s) for s in itertools.combinations(range(n), m)]
    child_caps = []
    for s in selectors:
        child = A[np.ix_(s, s)]
        child_caps.append(max(abs(energy(child, y)) for y in projective_spins(m)))

    favorable = np.zeros((len(centers), len(selectors)), dtype=bool)
    for iz, z in enumerate(centers):
        for js, s in enumerate(selectors):
            child = A[np.ix_(s, s)]
            favorable[iz, js] = abs(energy(child, z[list(s)])) == child_caps[js]
    rows = np.asarray([row_square(A, z) for z in centers], dtype=np.int64)
    return centers, selectors, child_caps, favorable, rows


def exact_box_minimum(m: int) -> tuple[int, tuple]:
    """Exhaust min V(S,y) over every child-ground incidence."""
    n = len(A)
    best = None
    witness = None
    for s in itertools.combinations(range(n), m):
        t = tuple(i for i in range(n) if i not in s)
        child = A[np.ix_(s, s)]
        words = projective_spins(m)
        vals = [energy(child, y) for y in words]
        q_s = max(map(abs, vals))
        for y, val in zip(words, vals):
            if abs(val) != q_s:
                continue
            for w in itertools.product((-1, 1), repeat=n - m):
                z = np.empty(n, dtype=np.int64)
                z[list(s)] = y
                z[list(t)] = w
                r = row_square(A, z)
                if best is None or r < best:
                    best = r
                    witness = (s, tuple(map(int, y)), t, w, q_s,
                               tuple(map(int, z)), energy(A, z))
    assert best is not None and witness is not None
    return best, witness


def aggregate_retention(
    favorable: np.ndarray,
    rows: np.ndarray,
    selectors: list[tuple[int, ...]],
    cap: int,
    ell: int,
) -> Fraction | None:
    """Triple retention after a^2 size bias, using the overlap-count formula."""
    n = len(A)
    m = len(selectors[0])
    kernel_denom = math.comb(m, ell) * math.comb(n - ell, m - ell)
    numerator = 0
    denominator = 0
    for iz in range(len(rows)):
        ids = np.flatnonzero(favorable[iz])
        r = len(ids)
        if rows[iz] > cap or not r:
            continue
        overlap_moment = sum(
            math.comb(len(set(selectors[i]) & set(selectors[j])), ell)
            for i in ids for j in ids
        )
        numerator += r * overlap_moment
        denominator += r * r
    if not denominator:
        return None
    return Fraction(numerator, kernel_denom * denominator)


def main() -> None:
    n, m = 10, 6
    centers = projective_spins(n)
    parent_energies = [energy(A, z) for z in centers]
    q = max(map(abs, parent_energies))
    assert q == 26  # The proved order-ten optimum is q_10=26.

    centers, selectors, child_caps, favorable, rows = build_ground_relation(m)
    active = favorable.any(axis=1)
    assert min(rows) == n
    # Each coordinate of A z is an odd integer, so R_2(z) >= n=10 for every z.
    assert all(all(int(v) % 2 for v in A @ z) for z in centers)

    box_min, witness = exact_box_minimum(m)
    assert box_min == 10
    expected_witness = (
        (0, 1, 2, 6, 7, 8),
        (1, 1, -1, -1, -1, 1),
        (3, 4, 5, 9),
        (1, -1, -1, 1),
        10,
        (1, 1, -1, 1, -1, -1, -1, -1, 1, 1),
        -10,
    )
    assert witness == expected_witness

    low = np.flatnonzero(rows <= box_min)
    assert len(low) == 2 and all(active[i] for i in low)
    families = [tuple(np.flatnonzero(favorable[i])) for i in low]
    assert families[0] == families[1]
    family_sets = [set(selectors[j]) for j in families[0]]
    assert len(family_sets) == 5
    intersection_histogram = Counter(
        len(s & t) for s in family_sets for t in family_sets
    )
    assert intersection_histogram == Counter({4: 10, 2: 10, 6: 5})

    expected = {
        1: (Fraction(2, 27), Fraction(1, 42), Fraction(9, 28), Fraction(369, 350), Fraction(2, 25)),
        2: (Fraction(1, 6), Fraction(29, 1050), Fraction(29, 175), Fraction(1021, 875), Fraction(32, 175)),
        3: (Fraction(2, 7), Fraction(1, 25), Fraction(7, 50), Fraction(168, 125), Fraction(8, 25)),
        4: (Fraction(4, 9), Fraction(17, 225), Fraction(17, 100), Fraction(208, 125), Fraction(64, 125)),
        5: (Fraction(2, 3), Fraction(1, 5), Fraction(3, 10), Fraction(12, 5), Fraction(4, 5)),
    }
    table = []
    for ell in range(1, m):
        lam1 = Fraction(ell * (n - m), m * (n - ell))
        lam2 = Fraction(
            falling(ell, 2) * falling(n - m, 2),
            falling(m, 2) * falling(n - ell, 2),
        ) if ell >= 2 else Fraction(0)
        kappa = (lam1 - lam2) / (1 - lam1)
        p = aggregate_retention(favorable, rows, selectors, box_min, ell)
        assert p is not None
        triple_ratio = p / lam1
        coarea_ratio = (1 - p) / (1 - lam1)
        assert (lam1, p, triple_ratio, coarea_ratio, kappa) == expected[ell]
        assert p < lam1 and coarea_ratio > 1
        pressure_gap = p - lam1
        table.append((ell, lam1, lam2, p, pressure_gap,
                      triple_ratio, coarea_ratio, kappa))

    # Under the uniform projective-center law the double-incidence mass is
    # exactly two copies of a=(5/210)=1/42.
    d_c = Fraction(2, len(centers)) * Fraction(1, 42) ** 2
    assert d_c == Fraction(1, 451584)
    expected_centers = {
        (1, 1, -1, -1, 1, -1, 1, 1, -1, 1),
        (1, 1, -1, 1, -1, -1, -1, -1, 1, 1),
    }
    assert {tuple(map(int, centers[i])) for i in low} == expected_centers

    # A K_0 admixture is a genuine exception to the positive-core mixture
    # obstruction.  K_0 independently resamples the selector, so p_0=a=1/42
    # on both active centers and all its nonconstant eigenvalues vanish.
    p0 = aggregate_retention(favorable, rows, selectors, box_min, 0)
    assert p0 == Fraction(1, 42)
    w0, w1 = Fraction(3, 4), Fraction(1, 4)
    p_mixed = w0 * p0 + w1 * expected[1][1]
    lambda1_mixed = w1 * expected[1][0]
    lambda2_mixed = Fraction(0)
    normalized_pressure = p_mixed - lambda1_mixed
    kappa_mixed = (lambda1_mixed - lambda2_mixed) / (1 - lambda1_mixed)
    coarea_mixed = (1 - p_mixed) / (1 - lambda1_mixed)
    assert p_mixed == Fraction(1, 42)
    assert lambda1_mixed == Fraction(1, 54)
    assert normalized_pressure == Fraction(1, 189)
    assert lambda2_mixed == 0
    assert kappa_mixed == Fraction(1, 53)
    assert coarea_mixed == Fraction(369, 371) < 1
    assert d_c * normalized_pressure == Fraction(1, 85349376)

    # Locate the first row cap at which each scale passes.  This also checks
    # that cap inflation, not a different scale at cap 10, is necessary here.
    active_caps = sorted(set(map(int, rows[active])))
    first_pass = {}
    for ell in range(1, m):
        lam1 = expected[ell][0]
        first_pass[ell] = next(
            (cap for cap in active_caps
             if aggregate_retention(favorable, rows, selectors, cap, ell) >= lam1),
            None,
        )
    assert active_caps == [10, 42, 74, 106, 138]
    assert first_pass == {1: 106, 2: 138, 3: None, 4: None, 5: None}

    print(f"explicit exact minimizer: n={n}, q={q}, m={m}")
    print(f"box minimum={box_min}; witness={witness}")
    print(f"row-{box_min} class: centers={len(low)}, degree=5/{len(selectors)}=1/42, D_C={d_c}")
    print(f"row-{box_min} center vectors={sorted(expected_centers)}")
    print(f"ground-family ordered intersection histogram={dict(sorted(intersection_histogram.items()))}")
    print("ell | lambda1 | lambda2 | triple p | p-lambda1 | p/lambda1 | coarea ratio | kappa")
    for row in table:
        print(" | ".join(map(str, row)))
    print(
        "Kmix=3/4 K0+1/4 K1: "
        f"p={p_mixed}, lambda1={lambda1_mixed}, "
        f"normalized_pressure={normalized_pressure}, lambda2={lambda2_mixed}, "
        f"kappa={kappa_mixed}, coarea={coarea_mixed}, "
        f"G={d_c * normalized_pressure}"
    )
    print(f"active row caps={active_caps}; first passing caps={first_pass}")
    print("PASS: positive-core mixtures fail at cap 10, while the audited K0 admixture passes")


if __name__ == "__main__":
    main()
