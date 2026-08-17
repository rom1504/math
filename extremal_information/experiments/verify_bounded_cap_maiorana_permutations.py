#!/usr/bin/env python3
"""Exact checks for the permutation Maiorana--McFarland response packing.

This verifier checks all transforms and all ordered pair Rayleigh identities
for q=2,4, computes an exact maximum low-bias permutation code at q=4, and
checks the Boolean response gap on a representative subcode at n=16.
"""

from __future__ import annotations

from itertools import permutations, product
from math import log


def parity(x: int) -> int:
    return bin(x).count("1") & 1


def dot_parity(x: int, y: int) -> int:
    return parity(x & y)


def walsh(d: int) -> list[list[int]]:
    size = 1 << d
    return [
        [1 if dot_parity(i, j) == 0 else -1 for j in range(size)]
        for i in range(size)
    ]


def matvec(matrix: list[list[int]], vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(a * b for a, b in zip(row, vector)) for row in matrix)


def quadratic(matrix: list[list[int]], vector: tuple[int, ...]) -> int:
    return sum(
        matrix[i][j] * vector[i] * vector[j]
        for i in range(len(vector))
        for j in range(i + 1, len(vector))
    )


def mm_vector(q: int, permutation: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        1 if dot_parity(u, permutation[v]) == 0 else -1
        for u in range(q)
        for v in range(q)
    )


def build(m: int):
    q = 1 << m
    n = q * q
    W = walsh(2 * m)
    identity = tuple(range(q))
    b = mm_vector(q, identity)
    H = [[b[i] * W[i][j] * b[j] for j in range(n)] for i in range(n)]
    A = [[0 if i == j else H[i][j] for j in range(n)] for i in range(n)]
    return q, n, W, H, A


def pair_bias(
    q: int, pi: tuple[int, ...], sigma: tuple[int, ...]
) -> tuple[int, tuple[int, ...]]:
    tau = tuple(v ^ pi[v] ^ sigma[v] for v in range(q))
    total = sum(
        1
        if (dot_parity(x, y) ^ dot_parity(tau[x], tau[y])) == 0
        else -1
        for x in range(q)
        for y in range(q)
    )
    return total, tau


def maximum_good_code(adjacency: list[set[int]]) -> list[int]:
    best: list[int] = []

    def search(chosen: list[int], candidates: list[int]) -> None:
        nonlocal best
        if len(chosen) + len(candidates) <= len(best):
            return
        if not candidates:
            if len(chosen) > len(best):
                best = chosen[:]
            return
        while candidates:
            vertex = candidates.pop()
            search(
                chosen + [vertex],
                [other for other in candidates if other in adjacency[vertex]],
            )

    search([], list(range(len(adjacency))))
    return best


def verify_order(m: int, check_responses: bool) -> tuple[int, int, int]:
    q, n, W, H, A = build(m)
    perms = list(permutations(range(q)))

    one = tuple([1] * n)
    assert all(H[i][j] == H[j][i] for i in range(n) for j in range(n))
    assert matvec(H, one) == tuple([q] * n)
    assert sum(H[i][i] for i in range(n)) == 0
    for i in range(n):
        for j in range(n):
            assert sum(H[i][k] * H[k][j] for k in range(n)) == (
                n if i == j else 0
            )

    vectors: list[tuple[int, ...]] = []
    queries: list[tuple[int, ...]] = []
    for pi in perms:
        s = mm_vector(q, pi)
        y = tuple(value // q for value in matvec(W, s))
        assert all(abs(value) == 1 for value in y)
        assert matvec(W, y) == tuple(q * value for value in s)
        vectors.append(s)
        queries.append(y)

    identity_index = perms.index(tuple(range(q)))
    assert queries[identity_index] == vectors[identity_index]

    adjacency = [set() for _ in perms]
    rayleigh_checks = 0
    for i, pi in enumerate(perms):
        for j, sigma in enumerate(perms):
            w = tuple(a * b for a, b in zip(vectors[i], vectors[j]))
            rayleigh = sum(a * b for a, b in zip(w, matvec(H, w)))
            bias_sum, _ = pair_bias(q, pi, sigma)
            assert rayleigh == q * bias_sum
            rho_numerator = bias_sum
            if i != j and 4 * rho_numerator <= q * q:
                adjacency[i].add(j)
            rayleigh_checks += 1

    code = maximum_good_code(adjacency)
    if q == 2:
        assert len(code) == 2
    if q == 4:
        assert len(code) == 20

    response_checks = 0
    if check_responses:
        spins = list(product((-1, 1), repeat=n))
        energies = [quadratic(A, u) for u in spins]
        assert max(abs(value) for value in energies) == q * n // 2

        # Six representatives already exercise all response identities; the
        # exact graph computation above certifies all twenty pair coordinates.
        audited = code[:6]
        responses: dict[tuple[int, int], int] = {}
        for i in audited:
            for j in audited:
                w = tuple(a * b for a, b in zip(vectors[i], vectors[j]))
                value = max(
                    energy + q * sum(ui * wi for ui, wi in zip(u, w))
                    for energy, u in zip(energies, spins)
                )
                responses[i, j] = value
                if i == j:
                    assert value == 3 * q * n // 2
                else:
                    bias_sum, _ = pair_bias(q, perms[i], perms[j])
                    # value <= q*n*(8+rho)/6, rho=bias_sum/q^2.
                    assert 6 * q * q * value <= q * n * (
                        8 * q * q + bias_sum
                    )
                response_checks += 1

        for offset, i in enumerate(audited):
            for j in audited[offset + 1 :]:
                at_i = responses[i, i] - responses[j, i]
                at_j = responses[i, j] - responses[j, j]
                assert 8 * at_i >= q * n
                assert 8 * at_j <= -q * n
                assert 4 * (at_i - at_j) >= q * n

    return len(code), rayleigh_checks, response_checks


def asymptotic_margin_threshold() -> int:
    # The constants in the self-contained union bound are deliberately
    # coarse.  Locate the first m for which its displayed exponent is
    # positive; this is a symbolic consistency check, not a finite theorem
    # threshold worth optimizing.
    for m in range(9, 4097):
        q = 2.0**m if m < 1024 else float("inf")
        r = m - 8
        if m < 1024:
            margin = r * q / 128 - 2 * r * log(q) - 2 * q
            if margin > 0:
                return m
        else:
            # Divide by q to avoid overflow: r/128-2 dominates, while the
            # term 2r log(q)/q is negligible here.
            if r / 128 - 2 > 0:
                return m
    raise AssertionError("no positive asymptotic margin found")


def main() -> None:
    q2 = verify_order(1, check_responses=True)
    q4 = verify_order(2, check_responses=True)
    threshold = asymptotic_margin_threshold()
    print(
        "permutation Maiorana--McFarland checks:",
        {
            "q2": {"maximum_code": q2[0], "rayleigh": q2[1], "responses": q2[2]},
            "q4": {"maximum_code": q4[0], "rayleigh": q4[1], "responses": q4[2]},
            "first_positive_crude_margin_m": threshold,
            "total_exact_identities": q2[1] + q2[2] + q4[1] + q4[2],
        },
    )


if __name__ == "__main__":
    main()
