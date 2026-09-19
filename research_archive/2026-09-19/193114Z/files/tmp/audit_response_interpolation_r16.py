#!/usr/bin/env python3
"""Exact convex-interpolation diagnostics for the A9 capture wall."""

from fractions import Fraction
from itertools import product

import numpy as np

from check_bicriteria_recurrence_r15 import A, QMIN, all_minimizers, qnorm


BLOCK = (0, 1, 2, 4, 5, 6)
SINGLETONS = ((3,), (7,), (8,))
SPINS = np.array([(1,) + t for t in product((-1, 1), repeat=8)], dtype=np.int64)


def signed_lines(M0: np.ndarray, M1: np.ndarray):
    """Return (intercept, slope, sigma, projective spin) energy lines."""
    e0 = np.einsum("bi,ij,bj->b", SPINS, M0, SPINS)
    e1 = np.einsum("bi,ij,bj->b", SPINS, M1, SPINS)
    lines = []
    for i, z in enumerate(SPINS):
        for sigma in (-1, 1):
            lines.append((int(sigma * e0[i]), int(sigma * (e1[i] - e0[i])), sigma, z))
    return lines


def exact_envelope_min(lines):
    """Minimize the upper envelope on [0,1] at all rational breakpoints."""
    candidates = {Fraction(0), Fraction(1)}
    for i, (a, b, _, _) in enumerate(lines):
        for c, d, _, _ in lines[i + 1 :]:
            if b != d:
                t = Fraction(c - a, b - d)
                if 0 <= t <= 1:
                    candidates.add(t)
    best_value = None
    best_t = None
    best_active = None
    for t in candidates:
        vals = [Fraction(a) + t * b for a, b, _, _ in lines]
        value = max(vals)
        if best_value is None or value < best_value:
            best_value = value
            best_t = t
            best_active = tuple(i for i, v in enumerate(vals) if v == value)
    return best_t, best_value, best_active


def main():
    C = A.copy()
    C[np.ix_(BLOCK, BLOCK)] = 0
    for singleton in SINGLETONS:
        C[np.ix_(singleton, singleton)] = 0
    D = A - C

    assert qnorm(A) == QMIN[9] == 24
    assert qnorm(C) == 22
    assert qnorm(A[np.ix_(BLOCK, BLOCK)]) == 22

    attaining = []
    for G6 in all_minimizers(6)[1]:
        K = np.zeros_like(A)
        K[np.ix_(BLOCK, BLOCK)] = G6
        H = C + K
        if qnorm(H) == 24:
            attaining.append((G6, K, H))
    assert len(attaining) == 40

    minima = {}
    constant_count = 0
    persistent_counts = {}
    persistent_profiles = {}
    for _, K, H in attaining:
        lines = signed_lines(A, H)
        t, value, active = exact_envelope_min(lines)
        minima[(t, value)] = minima.get((t, value), 0) + 1
        # Convexity and equal endpoint values imply Q <= 24 throughout.
        # Record whether one line is identically 24, which forces equality.
        persistent = [(sigma, tuple(z)) for a, b, sigma, z in lines if a == 24 and b == 0]
        if persistent:
            constant_count += 1
        persistent_counts[len(persistent)] = persistent_counts.get(len(persistent), 0) + 1
        e_c = np.einsum("bi,ij,bj->b", SPINS, C, SPINS)
        e_d = np.einsum("bi,ij,bj->b", SPINS, D, SPINS)
        e_k = np.einsum("bi,ij,bj->b", SPINS, K, SPINS)
        profiles = []
        for sigma, z_tuple in persistent:
            i = next(i for i, z in enumerate(SPINS) if tuple(z) == z_tuple)
            profiles.append((int(sigma * e_c[i]), int(sigma * e_d[i]), int(sigma * e_k[i])))
        key = tuple(sorted(profiles))
        persistent_profiles[key] = persistent_profiles.get(key, 0) + 1
        assert value <= 24

    print("attaining completions:", len(attaining))
    print("interpolation envelope minima (t,value):")
    for key, count in sorted(minima.items()):
        print(key, count)
    print("completions with an identically active +24 line:", constant_count)
    print("persistent signed-projective state count distribution:", sorted(persistent_counts.items()))
    print("persistent (cross, original, replacement) profile distributions:")
    for key, count in sorted(persistent_profiles.items()):
        print(key, count)
    assert constant_count == 40
    assert min(persistent_counts) == 4 and max(persistent_counts) == 13
    assert sum(persistent_counts.values()) == 40
    assert all(
        set(key) <= {(14, 10, 10), (18, 6, 6)}
        for key in persistent_profiles
    )

    # For the first completion, report endpoint active-state response data.
    _, K, H = attaining[0]
    e_c = np.einsum("bi,ij,bj->b", SPINS, C, SPINS)
    e_d = np.einsum("bi,ij,bj->b", SPINS, D, SPINS)
    e_k = np.einsum("bi,ij,bj->b", SPINS, K, SPINS)
    e_h = e_c + e_k
    rows = []
    for i, z in enumerate(SPINS):
        for sigma in (-1, 1):
            if sigma * e_h[i] == 24:
                rows.append((int(sigma * e_c[i]), int(sigma * e_d[i]), int(sigma * e_k[i])))
    print("first completion active (cross, original-internal, replacement):")
    print(sorted(rows))


if __name__ == "__main__":
    main()
