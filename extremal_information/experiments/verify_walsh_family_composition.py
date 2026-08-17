#!/usr/bin/env python3
"""Exact wind tunnel for repeated Walsh-bridge composition at n=4.

This is intentionally tiny: it enumerates every Boolean spin, every
Maiorana--McFarland truth table, and every word up to a requested depth.
It reports normalized response-message counts and tests candidate quotient
statistics.  The theorem checks added below use integer arithmetic only.
"""

from __future__ import annotations

from itertools import product

import numpy as np


def dot2(a: int, b: int) -> int:
    return bin(a & b).count("1") & 1


def walsh(m: int) -> list[list[int]]:
    pts = list(product(range(1 << m), repeat=2))
    return [
        [(-1) ** (dot2(a, u) ^ dot2(b, v)) for u, v in pts]
        for a, b in pts
    ]


def mat_vec(a: list[list[int]], x: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(aij * xj for aij, xj in zip(row, x)) for row in a)


def bilinear(x: tuple[int, ...], a: list[list[int]], y: tuple[int, ...]) -> int:
    ay = mat_vec(a, y)
    return sum(xi * zi for xi, zi in zip(x, ay))


def signs_for_g(m: int, gmask: int) -> tuple[int, ...]:
    q = 1 << m
    return tuple(
        (-1) ** (dot2(u, v) ^ ((gmask >> v) & 1))
        for u, v in product(range(q), repeat=2)
    )


def child_values(
    w: list[list[int]], spins: list[tuple[int, ...]], s: tuple[int, ...]
) -> tuple[int, ...]:
    # A_g is D_s (D_b W D_b) D_s with its diagonal removed.  Since the
    # diagonal trace is zero, 1/2 x^T D_s D_b W D_b D_s x is exact.
    # Here s=s_g and b=s_0, hence s*b=tau_g.
    b = signs_for_g(1, 0)
    tau = tuple(si * bi for si, bi in zip(s, b))
    vals = []
    for x in spins:
        z = tuple(ti * xi for ti, xi in zip(tau, x))
        vals.append(bilinear(z, w, z) // 2)
    return tuple(vals)


def normalize(message: tuple[int, ...]) -> tuple[int, ...]:
    c = message[0]
    return tuple(v - c for v in message)


def run_n4(depth: int = 12) -> None:
    m = 1
    q = 2
    n = q * q
    w = walsh(m)
    spins = list(product((-1, 1), repeat=n))
    bridge = [[bilinear(x, w, y) for y in spins] for x in spins]
    children = [child_values(w, spins, signs_for_g(m, g)) for g in range(1 << q)]

    states: dict[tuple[int, ...], tuple[int, ...]] = {
        normalize(children[g]): (g,) for g in range(1 << q)
    }
    print(f"depth=1 words={1 << q} projective_states={len(states)}")
    for k in range(2, depth + 1):
        nxt: dict[tuple[int, ...], tuple[int, ...]] = {}
        for msg, word in states.items():
            # A representative projective message is enough because max-plus
            # propagation commutes with additive constants.
            for g, hg in enumerate(children):
                out = tuple(
                    hg[j] + max(msg[i] + bridge[i][j] for i in range(len(spins)))
                    for j in range(len(spins))
                )
                nxt.setdefault(normalize(out), word + (g,))
        states = nxt
        print(f"depth={k} words={4**k} projective_states={len(states)}")

    # Pairwise signed truth-table correlations do not determine even the
    # unqueried three-block path optimum.
    by_pairwise: dict[tuple[int, int, int], dict[int, tuple[int, int, int]]] = {}
    tables = [tuple((-1) ** ((g >> v) & 1) for v in range(q)) for g in range(4)]
    for word in product(range(4), repeat=3):
        key = tuple(
            sum(a * b for a, b in zip(tables[word[i]], tables[word[j]]))
            for i, j in ((0, 1), (1, 2), (0, 2))
        )
        best = -10**9
        for i, x in enumerate(spins):
            for j, y in enumerate(spins):
                xy = children[word[0]][i] + bridge[i][j] + children[word[1]][j]
                best = max(
                    best,
                    max(xy + bridge[j][ell] + children[word[2]][ell]
                        for ell in range(len(spins))),
                )
        by_pairwise.setdefault(key, {})[best] = word
    witnesses = [(key, vals) for key, vals in by_pairwise.items() if len(vals) > 1]
    print(f"pairwise-correlation collision classes with distinct optima={len(witnesses)}")
    if witnesses:
        key, vals = witnesses[0]
        print("first pairwise collision", key, sorted(vals.items()))


def sylvester(q: int) -> np.ndarray:
    return np.array(
        [[(-1) ** dot2(a, u) for u in range(q)] for a in range(q)],
        dtype=np.int64,
    )


def linear_modulation(m: int, a: int) -> np.ndarray:
    q = 1 << m
    return np.array(
        [(-1) ** dot2(a, v) for u, v in product(range(q), repeat=2)],
        dtype=np.int64,
    )


def bent_pole(m: int, a: int) -> np.ndarray:
    q = 1 << m
    return np.array(
        [(-1) ** (dot2(u, v) ^ dot2(a, v))
         for u, v in product(range(q), repeat=2)],
        dtype=np.int64,
    )


def path_adjacency(k: int) -> np.ndarray:
    a = np.zeros((k, k), dtype=np.int64)
    for i in range(k - 1):
        a[i, i + 1] = a[i + 1, i] = 1
    return a


def verify_scalable_theorem() -> None:
    checks = 0
    for m in (2, 3):
        q = 1 << m
        n = q * q
        r = sylvester(q)
        w = np.kron(r, r)
        assert np.array_equal(r @ r, q * np.eye(q, dtype=np.int64))
        assert np.array_equal(w @ w, n * np.eye(n, dtype=np.int64))
        checks += 2

        # a_even has Hamming parity zero; a_odd has Hamming parity one.
        for a, sign in ((3, 1), (1, -1)):
            tau = linear_modulation(m, a)
            c = (tau[:, None] * w) * tau[None, :]
            assert np.array_equal(c @ c, n * np.eye(n, dtype=np.int64))
            assert np.array_equal(w @ c, sign * (c @ w))
            checks += 2

            if sign == 1:
                s = bent_pole(m, a)
                y_num = w @ s
                assert np.all(y_num % q == 0)
                y = y_num // q
                assert set(y.tolist()) == {-1, 1}
                assert np.array_equal(c @ s, q * s)
                assert np.array_equal(c @ y, q * y)
                checks += 4
                for k in (2, 3, 5):
                    xs = [s if i % 2 == 0 else y for i in range(k)]
                    energy2 = sum(int(x @ c @ x) for x in xs)
                    energy2 += 2 * sum(
                        int(xs[i] @ w @ xs[i + 1]) for i in range(k - 1)
                    )
                    expected2 = (3 * k - 2) * n * q
                    assert energy2 == expected2
                    checks += 1
            else:
                for k in (2, 3, 5):
                    adjacency = path_adjacency(k)
                    global_matrix = np.kron(np.eye(k, dtype=np.int64), c)
                    global_matrix += np.kron(adjacency, w)
                    target_square = n * np.kron(
                        np.eye(k, dtype=np.int64) + adjacency @ adjacency,
                        np.eye(n, dtype=np.int64),
                    )
                    assert np.array_equal(global_matrix @ global_matrix, target_square)
                    observed = max(abs(np.linalg.eigvalsh(global_matrix))) / q
                    predicted = (1 + 4 * np.cos(np.pi / (k + 1)) ** 2) ** 0.5
                    assert abs(observed - predicted) < 1e-10
                    even_coeff = 1.5 * k - 1
                    odd_coeff = 0.5 * k * predicted
                    assert even_coeff > odd_coeff
                    checks += 3

        # Exact Kronecker carrier for an arbitrary three-label path.
        k = 3
        tables = [
            np.array([1 if ((mask >> v) & 1) == 0 else -1 for v in range(q)],
                     dtype=np.int64)
            for mask in (1, 3, 5)
        ]
        carrier = np.zeros((k * q, k * q), dtype=np.int64)
        for i, table in enumerate(tables):
            carrier[i*q:(i+1)*q, i*q:(i+1)*q] = (
                table[:, None] * r * table[None, :]
            )
        for i in range(k - 1):
            carrier[i*q:(i+1)*q, (i+1)*q:(i+2)*q] = r
            carrier[(i+1)*q:(i+2)*q, i*q:(i+1)*q] = r
        lifted = np.kron(r, carrier)
        # Audit every entry against its child/bridge formula in (i,u,v) order.
        for u, jv, up, lp in product(range(q), range(k*q), range(q), range(k*q)):
            i, v = divmod(jv, q)
            ell, vp = divmod(lp, q)
            if i == ell:
                expected = r[u, up] * tables[i][v] * r[v, vp] * tables[i][vp]
            elif abs(i - ell) == 1:
                expected = r[u, up] * r[v, vp]
            else:
                expected = 0
            assert lifted[u*k*q+jv, up*k*q+lp] == expected
            checks += 1

    print(f"scalable Walsh carrier/holonomy checks passed: {checks}")


if __name__ == "__main__":
    verify_scalable_theorem()
    run_n4()
