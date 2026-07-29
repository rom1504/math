#!/usr/bin/env python3
"""Exact checks for finite_fibre_renormalization.md."""

from itertools import combinations, product

import numpy as np


def boolean_norm(matrix):
    """Return max |x^T M x| and one witness, by exhaustive enumeration."""
    n = matrix.shape[0]
    best = (-1, None, None)
    for spin in product((-1, 1), repeat=n):
        x = np.asarray(spin, dtype=np.int64)
        energy = int(x @ matrix @ x)
        candidate = (abs(energy), spin, energy)
        if candidate[0] > best[0]:
            best = candidate
    return best


def fibre_boolean_norm(base, r, d):
    """Enumerate a three-spin fibre lift without traversing 2^(3k) tuples."""
    states = np.asarray(list(product((-1, 1), repeat=3)), dtype=np.int64)
    lift = np.kron(base, r) + np.kron(
        np.eye(base.shape[0], dtype=np.int64), d
    )
    best = (-1, None, None)
    for labels in product(range(8), repeat=base.shape[0]):
        x = states[list(labels)].reshape(-1)
        energy = int(x @ lift @ x)
        candidate = (abs(energy), labels, energy)
        if candidate[0] > best[0]:
            best = candidate
    return best, lift


def signing_from_edges(n, edge_signs):
    matrix = np.zeros((n, n), dtype=np.int64)
    for (i, j), sign in zip(combinations(range(n), 2), edge_signs):
        matrix[i, j] = matrix[j, i] = sign
    return matrix


def main():
    r = np.ones((3, 3), dtype=np.int64) - 2 * np.eye(3, dtype=np.int64)
    d = np.ones((3, 3), dtype=np.int64) - np.eye(3, dtype=np.int64)

    reps = np.asarray(
        [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1)],
        dtype=np.int64,
    )
    kernel = reps @ r @ reps.T
    internal = np.einsum("bi,ij,bj->b", reps, d, reps)
    expected_kernel = np.asarray(
        [[3, 1, 1, 1], [1, -5, 3, 3], [1, 3, -5, 3], [1, 3, 3, -5]],
        dtype=np.int64,
    )
    assert np.array_equal(kernel, expected_kernel)
    assert tuple(internal) == (6, -2, -2, -2)

    c6 = np.asarray(
        [
            [0, 1, 1, 1, 1, 1],
            [1, 0, 1, -1, -1, 1],
            [1, 1, 0, 1, -1, -1],
            [1, -1, 1, 0, 1, -1],
            [1, -1, -1, 1, 0, 1],
            [1, 1, -1, -1, 1, 0],
        ],
        dtype=np.int64,
    )
    assert np.array_equal(c6 @ c6, 5 * np.eye(6, dtype=np.int64))
    assert boolean_norm(c6)[0] == 10

    (lift_q, lift_witness, lift_energy), lift = fibre_boolean_norm(c6, r, d)
    assert lift_q == 78
    assert lift_witness == (1, 1, 3, 2, 4, 5)
    assert lift_energy == -78
    assert np.all(np.diag(lift) == 0)
    off_diagonal = lift[~np.eye(lift.shape[0], dtype=bool)]
    assert set(off_diagonal) == {-1, 1}

    max_t_by_s = {}
    for tail in product((-1, 1), repeat=5):
        z = np.asarray((1,) + tail, dtype=np.int64)
        switched = c6 * np.outer(z, z)
        s_value = sum(switched[i, j] for i, j in combinations(range(6), 2))
        max_t = -10**9
        for colours in product(range(3), repeat=6):
            t_value = sum(
                switched[i, j]
                for i, j in combinations(range(6), 2)
                if colours[i] == colours[j]
            )
            max_t = max(max_t, t_value)
        max_t_by_s[s_value] = max(max_t_by_s.get(s_value, -10**9), max_t)
    assert max_t_by_s == {-5: 2, -3: 3, 3: 4, 5: 5}
    assert min(6 * s - 16 * t - 12 for s, t in max_t_by_s.items()) == -78

    positive_spectral_bound = 18 * (2 + np.sqrt(5))
    assert positive_spectral_bound < 78

    assert boolean_norm(r)[0] == 5
    assert boolean_norm(np.kron(r, r))[0] == 33

    same_q_lift_values = []
    for edges in [(-1, -1, -1), (-1, -1, 1)]:
        base = signing_from_edges(3, edges)
        base_q = boolean_norm(base)[0]
        lifted_q = fibre_boolean_norm(base, r, d)[0][0]
        same_q_lift_values.append((base_q, lifted_q))
    assert same_q_lift_values == [(6, 24), (6, 36)]

    print("four-state kernel and internal vector: verified")
    print("C6^2=5I and Q(C6)=10: verified")
    print(
        "Q(T(C6))=78; doubled/original normalizations:",
        78 / (18 ** 1.5),
        78 / (2 * 18 ** 1.5),
    )
    print("(S,max T) table:", dict(sorted(max_t_by_s.items())))
    print("Q(R)=5 and Q(R tensor R)=33: verified")
    print("same-Q scalar-recursion counterexample:", same_q_lift_values)


if __name__ == "__main__":
    main()
