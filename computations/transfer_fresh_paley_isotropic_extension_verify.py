#!/usr/bin/env python3
"""Exact integer checks for the Paley isotropic-ground extension theorem.

No random searches or optimization certificates are inferred from samples.
The r=3 row-extension calculation exhausts both row and spin cubes modulo
their separate global signs. Larger r checks only the explicit ground law.
"""

import itertools
import json

import numpy as np


def paley_square(r):
    squares = {a * a % r for a in range(1, r)}
    d = next(a for a in range(1, r) if a not in squares)
    q = r * r

    def mul(z, w):
        a, b = z % r, z // r
        c, e = w % r, w // r
        return (a * c + d * b * e) % r + r * ((a * e + b * c) % r)

    field_squares = {mul(z, z) for z in range(1, q)}
    nonsquare = next(z for z in range(1, q) if z not in field_squares)
    C = np.zeros((q + 1, q + 1), dtype=np.int64)
    C[0, 1:] = 1
    C[1:, 0] = 1
    for z in range(q):
        for w in range(q):
            if z != w:
                difference = (z % r - w % r) % r + r * ((z // r - w // r) % r)
                C[z + 1, w + 1] = 1 if difference in field_squares else -1

    positive = []
    for u in sorted(field_squares):
        for plus in itertools.combinations(range(r), (r + 1) // 2):
            f = np.full(r, -1, dtype=np.int64)
            f[list(plus)] = 1
            x = np.ones(q + 1, dtype=np.int64)
            for z in range(q):
                x[z + 1] = f[mul(u, z) // r]
            positive.append(x)
    positive = np.asarray(positive)
    negative = np.empty_like(positive)
    negative[:, 0] = -positive[:, 0]
    for z in range(q):
        negative[:, mul(nonsquare, z) + 1] = positive[:, z + 1]
    return C, positive, negative


def check(r):
    C, positive, negative = paley_square(r)
    n = len(C)
    count = len(positive)
    identity = np.eye(n, dtype=np.int64)
    assert np.array_equal(C @ C, r * r * identity)
    assert np.array_equal(positive @ C, r * positive)
    assert np.array_equal(negative @ C, -r * negative)
    assert np.array_equal(r * (positive.T @ positive), count * (r * identity + C))
    assert np.array_equal(r * (negative.T @ negative), count * (r * identity - C))
    combined = np.vstack((positive, negative))
    assert np.array_equal(combined.T @ combined, len(combined) * identity)
    cap = r * n // 2
    result = {
        "r": r,
        "n": n,
        "ground_cap": cap,
        "positive_law_atoms_with_multiplicity": count,
        "isotropic_identity_verified": True,
        "all_rows_extension_lower_bound": cap + r + 1,
    }
    if r == 3:
        X = np.asarray([(1,) + x for x in itertools.product((-1, 1), repeat=n - 1)], dtype=np.int64)
        energies = np.einsum("bi,ij,bj->b", X, C, X) // 2
        assert int(np.max(np.abs(energies))) == cap
        correlations = np.abs(X @ X.T)
        row_caps = np.max(correlations + np.abs(energies)[None, :], axis=1)
        exact_extension = int(np.min(row_caps))
        ground_test = np.max(np.abs(X @ combined.T), axis=1)
        result.update({
            "all_rows_modulo_global_sign": len(X),
            "all_spins_modulo_global_sign": len(X),
            "exact_best_extension": exact_extension,
            "minimum_ground_only_correlation": int(np.min(ground_test)),
            "optimal_row_count_modulo_global_sign": int(np.sum(row_caps == exact_extension)),
            "one_optimal_row": X[int(np.argmin(row_caps))].tolist(),
        })
        assert exact_extension == cap + r + 1
    return result


def check_prime_paley(q):
    assert q % 4 == 1
    chi = np.asarray([0] + [1 if pow(z, (q - 1) // 2, q) == 1 else -1 for z in range(1, q)], dtype=np.int64)
    C = np.zeros((q + 1, q + 1), dtype=np.int64)
    C[0, 1:] = 1
    C[1:, 0] = 1
    for z in range(q):
        for w in range(q):
            C[z + 1, w + 1] = chi[(z - w) % q]
    n = q + 1
    assert np.array_equal(C @ C, q * np.eye(n, dtype=np.int64))

    translate = [0] + [(z + 1) % q + 1 for z in range(q)]
    assert np.array_equal(C[np.ix_(translate, translate)], C)
    invert = [1, 0] + [pow(z, -1, q) + 1 for z in range(1, q)]
    signs = np.asarray([1, 1] + [int(chi[z]) for z in range(1, q)])
    assert np.array_equal(C[np.ix_(invert, invert)] * np.outer(signs, signs), C)
    nonsquare = next(z for z in range(1, q) if chi[z] == -1)
    scale = [0] + [nonsquare * z % q + 1 for z in range(q)]
    antisigns = np.asarray([-1] + [1] * q)
    assert np.array_equal(C[np.ix_(scale, scale)] * np.outer(antisigns, antisigns), -C)

    X = np.asarray([(1,) + x for x in itertools.product((-1, 1), repeat=n - 1)], dtype=np.int64)
    energies = np.einsum("bi,ij,bj->b", X, C, X) // 2
    cap = int(np.max(np.abs(energies)))
    absolute_grounds = X[np.abs(energies) == cap]
    positive_grounds = X[energies == cap]
    negative_grounds = X[energies == -cap]
    assert len(positive_grounds) == len(negative_grounds)
    assert np.array_equal(absolute_grounds.T @ absolute_grounds,
                          len(absolute_grounds) * np.eye(n, dtype=np.int64))
    D = n * (n - 1) // 2
    assert np.array_equal(D * (positive_grounds.T @ positive_grounds),
                          len(positive_grounds) * (D * np.eye(n, dtype=np.int64) + cap * C))
    result = {
        "q": q,
        "n": n,
        "exact_cap": cap,
        "absolute_ground_count_modulo_global_sign": len(absolute_grounds),
        "symmetries_and_isotropic_ground_law_verified": True,
    }
    if q == 5:
        correlations = np.abs(X @ X.T)
        result["exact_best_extension"] = int(np.min(np.max(correlations + np.abs(energies)[None, :], axis=1)))
        assert result["exact_best_extension"] == 9
    return result


if __name__ == "__main__":
    print(json.dumps({"square_field_checks": [check(r) for r in (3, 5, 7)],
                      "prime_field_checks": [check_prime_paley(q) for q in (5, 13, 17)]}, indent=2))
