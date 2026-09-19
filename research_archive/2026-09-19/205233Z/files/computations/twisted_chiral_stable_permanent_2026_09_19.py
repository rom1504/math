#!/usr/bin/env python3
"""Exact integer audit of conditional local-stability permanents.

No solver, probabilistic certificate, or asymptotic inference. The numerator
polynomial counts all sector permutations and all matching sign vectors.
"""
import itertools
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def convolve(left, right):
    out = Counter()
    for a, x in left.items():
        for b, y in right.items():
            out[a + b] += x * y
    return out


def polynomial_permanent(rows, cols, margins):
    states = {0: Counter({0: 1})}
    for i in rows:
        new = {}
        for mask, polynomial in states.items():
            for bit, j in enumerate(cols):
                if mask >> bit & 1:
                    continue
                allowed = Counter()
                if margins[i, j] >= 0:
                    allowed[1] = 1
                if margins[i, j] >= 2:
                    allowed[-1] = 1
                if not allowed:
                    continue
                target = new.setdefault(mask | (1 << bit), Counter())
                target.update(convolve(polynomial, allowed))
        states = new
    return states.get((1 << len(cols)) - 1, Counter())


def main():
    rng = np.random.default_rng(2026091917)
    records = []
    checked_permutations = 0
    checked_matchings = 0
    for n in range(2, 7):
        for trial in range(10 if n <= 5 else 6):
            raw = rng.choice([-1, 1], (n, n))
            A = np.triu(raw, 1)
            A = A + A.T
            x, y, z = rng.choice([-1, 1], (3, n))
            if trial == 0:
                y = x.copy()
            if trial == 1:
                y = -x
            signs = x * y
            plus = np.flatnonzero(signs == 1).tolist()
            minus = np.flatnonzero(signs == -1).tolist()
            perm = rng.permutation(n)
            image_plus = sorted(perm[plus].tolist())
            image_minus = sorted(perm[minus].tolist())
            relative = -np.ones(n, dtype=np.int64)
            relative[image_plus] = 1
            w = z * relative
            alpha, gamma = x * (A @ x), y * (A @ y)
            b, c = z * (A @ w), w * (A @ z)
            margins = np.minimum(alpha[:, None] + b[None, :],
                                 -gamma[:, None] + c[None, :])
            assert np.all(margins % 2 == 0)
            predicted = convolve(polynomial_permanent(plus, image_plus, margins),
                                 polynomial_permanent(minus, image_minus, margins))
            base = int((x @ A @ x - y @ A @ y) // 2 + z @ A @ w)
            us = np.array(list(itertools.product([-1, 1], repeat=n)), dtype=np.int64)
            sums = us.sum(axis=1)
            actual = Counter()
            permutation_count = 0
            parent_spin = np.concatenate([x, y])
            for pp in itertools.permutations(image_plus):
                for pm in itertools.permutations(image_minus):
                    p = np.zeros(n, dtype=np.int64)
                    p[plus], p[minus] = pp, pm
                    s = z[p] * x
                    B = A[np.ix_(p, p)] * np.outer(s, s)
                    core = np.block([[A, B], [B, -A]])
                    aligned = parent_spin * (core @ parent_spin)
                    assert np.array_equal(aligned[:n], alpha + b[p])
                    assert np.array_equal(aligned[n:], -gamma + c[p])
                    assert int(parent_spin @ core @ parent_spin // 2) == base
                    stable_direct = np.all(aligned[:n] + us > 0, axis=1)
                    stable_direct &= np.all(aligned[n:] + us > 0, axis=1)
                    stable_formula = np.all(margins[np.arange(n), p] + us > 0, axis=1)
                    assert np.array_equal(stable_direct, stable_formula)
                    for value in sums[stable_direct]:
                        actual[int(value)] += 1
                    # Also explicitly assemble one full parent matching per p.
                    u = us[permutation_count % len(us)]
                    C = B + np.diag(u * signs)
                    D = np.block([[A, C], [C, -A]])
                    assert np.all((parent_spin * (D @ parent_spin)) % 2 == 1)
                    assert int(parent_spin @ D @ parent_spin // 2) == base + int(u.sum())
                    permutation_count += 1
            denominator = math.factorial(len(plus)) * math.factorial(len(minus)) * (1 << n)
            assert permutation_count == math.factorial(len(plus)) * math.factorial(len(minus))
            assert predicted == actual
            assert sum(actual.values()) <= denominator
            checked_permutations += permutation_count
            checked_matchings += permutation_count * (1 << n)
            records.append({
                "n": n, "trial": trial, "matrix": A.tolist(),
                "x": x.tolist(), "y": y.tolist(), "z": z.tolist(), "w": w.tolist(),
                "base_energy": base, "sector_sizes": [len(plus), len(minus)],
                "row_margins": margins.tolist(),
                "stable_matching_sum_numerator": dict(sorted(actual.items())),
                "probability_denominator": denominator,
                "sector_permutations": permutation_count,
                "exact_identity_pass": True,
            })
    result = {"status": "exact finite identity checks; no asymptotic conclusion",
              "seed": 2026091917, "case_count": len(records),
              "sector_permutations_checked": checked_permutations,
              "permutation_matching_pairs_checked": checked_matchings,
              "all_pass": True, "cases": records}
    path = ROOT / "computations/results/twisted_chiral_stable_permanent_2026_09_19.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "cases"}))


if __name__ == "__main__":
    main()
