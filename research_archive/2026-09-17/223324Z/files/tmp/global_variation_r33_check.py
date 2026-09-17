#!/usr/bin/env python3
"""Exact finite checks for Wave 33 global block variation.

All scratch/output stays in the repository tmp directory.  The program:
  * verifies the arbitrary-edge flip identity and global stability on named
    exact minimizers;
  * exhausts small q_b and the block-deletion inequality;
  * solves the continuous response game and its dual law;
  * verifies the A9 witness-migration facts and the abstract orientation wall.
"""

from __future__ import annotations

from itertools import combinations, product
from math import comb, log, sqrt

import numpy as np
from scipy.optimize import linprog


A4 = np.array([
    [0, -1, 1, 1],
    [-1, 0, 1, -1],
    [1, 1, 0, -1],
    [1, -1, -1, 0],
], dtype=np.int64)

A6 = np.array([
    [0, 1, 1, 1, 1, 1],
    [1, 0, -1, -1, 1, 1],
    [1, -1, 0, 1, -1, 1],
    [1, -1, 1, 0, 1, -1],
    [1, 1, -1, 1, 0, -1],
    [1, 1, 1, -1, -1, 0],
], dtype=np.int64)

A9 = np.array([
    [0, 1, -1, 1, 1, 1, 1, -1, -1],
    [1, 0, 1, -1, -1, -1, 1, -1, -1],
    [-1, 1, 0, 1, 1, 1, 1, 1, -1],
    [1, -1, 1, 0, 1, 1, 1, -1, 1],
    [1, -1, 1, 1, 0, 1, -1, 1, -1],
    [1, -1, 1, 1, 1, 0, -1, -1, -1],
    [1, 1, 1, 1, -1, -1, 0, 1, 1],
    [-1, -1, 1, -1, 1, -1, 1, 0, -1],
    [-1, -1, -1, 1, -1, -1, 1, -1, 0],
], dtype=np.int64)


def states(n: int):
    """The 2^n oriented projective full states (sigma,x), x[0]=1."""
    for sigma in (-1, 1):
        for tail in product((-1, 1), repeat=n - 1):
            yield sigma, np.asarray((1,) + tail, dtype=np.int64)


def qnorm(a: np.ndarray) -> int:
    return max(int(sigma * x @ a @ x) for sigma, x in states(len(a)))


def signing_from_bits(b: int, bits: tuple[int, ...]) -> np.ndarray:
    ans = np.zeros((b, b), dtype=np.int64)
    for bit, (i, j) in zip(bits, combinations(range(b), 2)):
        ans[i, j] = ans[j, i] = bit
    return ans


def exact_q_small(b: int) -> int:
    return min(
        qnorm(signing_from_bits(b, bits))
        for bits in product((-1, 1), repeat=comb(b, 2))
    )


def response_profiles(a: np.ndarray, u: tuple[int, ...]):
    """Return local edge vectors d and exact outside responses Z."""
    n = len(a)
    b = len(u)
    outside = tuple(i for i in range(n) if i not in u)
    edges = tuple(combinations(range(b), 2))
    d_rows = []
    responses = []
    labels = []
    auo = a[np.ix_(u, outside)]
    aoo = a[np.ix_(outside, outside)]
    for sigma in (-1, 1):
        for tail in product((-1, 1), repeat=b - 1):
            x = np.asarray((1,) + tail, dtype=np.int64)
            d_rows.append([sigma * x[i] * x[j] for i, j in edges])
            best = -10**9
            for y_tuple in product((-1, 1), repeat=len(outside)):
                y = np.asarray(y_tuple, dtype=np.int64)
                value = sigma * (y @ aoo @ y + 2 * x @ auo @ y)
                best = max(best, int(value))
            responses.append(best)
            labels.append((sigma, x))
    return edges, np.asarray(d_rows, dtype=float), np.asarray(responses, dtype=float), labels


def solve_response_game(a: np.ndarray, u: tuple[int, ...]):
    edges, d, zresp, labels = response_profiles(a, u)
    nstates, medges = d.shape

    # Primal: min t subject to Z_u + 2 z.d_u <= t, -1 <= z_e <= 1.
    c = np.r_[np.zeros(medges), 1.0]
    aub = np.c_[2.0 * d, -np.ones(nstates)]
    primal = linprog(
        c,
        A_ub=aub,
        b_ub=-zresp,
        bounds=[(-1.0, 1.0)] * medges + [(None, None)],
        method="highs",
    )
    assert primal.success

    # Dual written explicitly: max E Z - 2 sum r_e,
    # r_e >= |E d_e|, mu a probability law.
    cdual = np.r_[-zresp, 2.0 * np.ones(medges)]
    dual_rows = []
    for e in range(medges):
        row = np.zeros(nstates + medges)
        row[:nstates] = d[:, e]
        row[nstates + e] = -1.0
        dual_rows.append(row)
        row = np.zeros(nstates + medges)
        row[:nstates] = -d[:, e]
        row[nstates + e] = -1.0
        dual_rows.append(row)
    dual = linprog(
        cdual,
        A_ub=np.asarray(dual_rows),
        b_ub=np.zeros(2 * medges),
        A_eq=np.asarray([np.r_[np.ones(nstates), np.zeros(medges)]]),
        b_eq=np.asarray([1.0]),
        bounds=[(0.0, None)] * (nstates + medges),
        method="highs",
    )
    assert dual.success
    assert abs(primal.fun + dual.fun) < 1e-8

    mu = dual.x[:nstates]
    means = mu @ d
    actual = np.asarray([a[u[i], u[j]] for i, j in edges], dtype=float)
    local = 2.0 * d @ actual
    full = zresp + local
    q = float(np.max(full))
    v = float(primal.fun)
    eta = sqrt(8.0 * medges * len(u) * log(2.0))
    alignment = sum(max(actual[e] * means[e], 0.0) for e in range(medges))
    mean_slack = float(mu @ (q - full))
    mean_local = float(mu @ local)
    mean_response = float(mu @ zresp)

    assert v <= q + 1e-8
    assert v >= q - eta - 1e-8
    assert 4.0 * alignment <= eta + 1e-7
    assert mean_slack <= eta + 1e-7
    assert mean_local <= eta / 2.0 + 1e-7
    assert mean_response >= q - eta - 1e-7
    return {
        "q": q,
        "V": v,
        "Q0": float(np.max(zresp)),
        "eta": eta,
        "alignment": alignment,
        "mean_slack": mean_slack,
        "mean_local": mean_local,
        "mean_response": mean_response,
        "support": int(np.sum(mu > 1e-8)),
        "labels": labels,
    }


def check_flip_identity(a: np.ndarray, u: tuple[int, ...], expected_q: int) -> None:
    full_states = list(states(len(a)))
    energies = np.asarray([sigma * x @ a @ x for sigma, x in full_states], dtype=int)
    assert int(energies.max()) == expected_q
    delta = expected_q - energies
    edges = tuple(combinations(u, 2))
    signed = np.asarray(
        [[sigma * a[i, j] * x[i] * x[j] for i, j in edges]
         for sigma, x in full_states],
        dtype=int,
    )
    for mask in range(1 << len(edges)):
        fids = [e for e in range(len(edges)) if (mask >> e) & 1]
        changed = a.copy()
        for e in fids:
            i, j = edges[e]
            changed[i, j] *= -1
            changed[j, i] *= -1
        lhs = qnorm(changed) - expected_q
        if fids:
            rhs = int(np.max(-4 * np.sum(signed[:, fids], axis=1) - delta))
        else:
            rhs = int(np.max(-delta))
        assert lhs == rhs
        assert lhs >= 0  # named A is an exact global minimizer

    # The singleton specialization has a <=4-slack opposing witness.
    for e in range(len(edges)):
        assert np.any((signed[:, e] == -1) & (delta <= 4))


def zero_internal(a: np.ndarray, u: tuple[int, ...]) -> np.ndarray:
    ans = a.copy()
    for i, j in combinations(u, 2):
        ans[i, j] = ans[j, i] = 0
    return ans


def check_block_deletion(a: np.ndarray, u: tuple[int, ...], expected_q: int, qb: int) -> None:
    q0 = qnorm(zero_internal(a, u))
    assert expected_q - q0 <= qb


def check_a9_migration() -> None:
    q = qnorm(A9)
    assert q == 24
    full_states = list(states(9))
    energies = np.asarray([sigma * x @ A9 @ x for sigma, x in full_states])

    u = (0, 2, 4, 5, 7)
    internal = A9 - zero_internal(A9, u)
    ground_local = {
        int(sigma * x @ internal @ x)
        for (sigma, x), energy in zip(full_states, energies)
        if energy == q
    }
    assert ground_local == {8}
    assert qnorm(zero_internal(A9, u)) == 24

    u_shield = (0, 1, 2, 6)
    assert qnorm(zero_internal(A9, u_shield)) == 32


def check_abstract_orientation_wall(b: int = 4) -> None:
    medges = comb(b, 2)
    q = 20
    zplus = q - 2 * medges
    zminus = q + 2 * medges
    for beta in product((-1, 1), repeat=medges):
        total = sum(beta)
        scores = (zplus + 2 * total, zminus - 2 * total)
        assert max(scores) >= q
    original_total = medges
    assert zplus + 2 * original_total == q
    assert zminus - 2 * original_total == q


def main() -> None:
    known_q = {2: 2, 3: 6, 4: 8, 5: 8}
    for b, expected in known_q.items():
        got = exact_q_small(b)
        assert got == expected, (b, got, expected)

    assert qnorm(A4) == 8
    assert qnorm(A6) == 10
    assert qnorm(A9) == 24

    # Arbitrary simultaneous flips in representative blocks.
    check_flip_identity(A4, (0, 1, 2, 3), 8)
    check_flip_identity(A6, (0, 1, 2, 3), 10)
    check_flip_identity(A9, (2, 3, 4, 6, 7), 24)

    # The prior exact block-deletion inequality q-q0 <= q_b.
    check_block_deletion(A4, (0, 1, 2, 3), 8, known_q[4])
    check_block_deletion(A6, (0, 1, 2, 3), 10, known_q[4])
    check_block_deletion(A9, (2, 3, 4, 6, 7), 24, known_q[5])
    check_block_deletion(A9, (0, 2, 4, 5, 7), 24, known_q[5])

    audits = {
        "A4/all": solve_response_game(A4, (0, 1, 2, 3)),
        "A6/0123": solve_response_game(A6, (0, 1, 2, 3)),
        "A9/loss4": solve_response_game(A9, (2, 3, 4, 6, 7)),
        "A9/migration": solve_response_game(A9, (0, 2, 4, 5, 7)),
        "A9/shield": solve_response_game(A9, (0, 1, 2, 6)),
    }
    check_a9_migration()
    check_abstract_orientation_wall()

    print("exact q_b:", known_q)
    for name, data in audits.items():
        compact = {key: round(value, 9) if isinstance(value, float) else value
                   for key, value in data.items() if key != "labels"}
        print(name, compact)
    print("global_variation_r33_check: PASS")


if __name__ == "__main__":
    main()

