#!/usr/bin/env python3
"""Exact finite regression checks for the bounded-cap linear-rate theorem.

The theorem's asymptotic inputs are standard concentration inequalities and
are not proved by finite computation.  This script checks every deterministic
identity on which the probabilistic argument rests:

* regular-Walsh construction and exact cap at n=4 and n=16;
* the weighted top-deficit formula BC.1 for all small integer fields at n=4;
* the two-query row reduction (S,T)=(U+V,U-V), including both adversarial
  target signs, for every possible overlap through n=10;
* the projective query-linked exposure inequality for every 4 by 4 sign
  bridge in a fixed nontrivial sample and every ordered query pair.

All arithmetic is integral except the final displayed normalizations.
"""

from __future__ import annotations

from itertools import product


def parity(x: int) -> int:
    return bin(x).count("1") & 1


def dot_parity(x: int, y: int) -> int:
    return parity(x & y)


def sign(x: int) -> int:
    return 1 if x >= 0 else -1


def walsh(d: int) -> list[list[int]]:
    size = 1 << d
    return [
        [1 if dot_parity(i, j) == 0 else -1 for j in range(size)]
        for i in range(size)
    ]


def mm_base(m: int) -> tuple[int, ...]:
    q = 1 << m
    return tuple(
        1 if dot_parity(u, v) == 0 else -1
        for u in range(q)
        for v in range(q)
    )


def matvec(matrix: list[list[int]], vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(a * b for a, b in zip(row, vector)) for row in matrix)


def quadratic(matrix: list[list[int]], vector: tuple[int, ...]) -> int:
    return sum(
        matrix[i][j] * vector[i] * vector[j]
        for i in range(len(vector))
        for j in range(i + 1, len(vector))
    )


def regular_walsh(m: int):
    q = 1 << m
    n = q * q
    W = walsh(2 * m)
    b = mm_base(m)
    H = [[b[i] * W[i][j] * b[j] for j in range(n)] for i in range(n)]
    A = [[0 if i == j else H[i][j] for j in range(n)] for i in range(n)]
    return q, n, W, H, A


def verify_regular_child(m: int, exhaustive: bool) -> int:
    q, n, _, H, A = regular_walsh(m)
    one = tuple([1] * n)
    assert matvec(H, one) == tuple([q] * n)
    assert sum(H[i][i] for i in range(n)) == 0
    for i in range(n):
        for j in range(n):
            assert sum(H[i][k] * H[k][j] for k in range(n)) == (
                n if i == j else 0
            )
    if exhaustive:
        spins = list(product((-1, 1), repeat=n))
        energies = [quadratic(A, u) for u in spins]
        assert max(abs(value) for value in energies) == q * n // 2
        return len(spins)
    # The spectral proof is exact once H^2=nI and H1=q1 are checked.
    return 0


def verify_weighted_deficit() -> int:
    q, n, _, _, A = regular_walsh(1)
    spins = list(product((-1, 1), repeat=n))
    energies = {u: quadratic(A, u) for u in spins}
    top = max(energies.values())
    checks = 0
    for h in product(range(-2, 3), repeat=n):
        direct = top + sum(abs(value) for value in h) - max(
            energies[u] + sum(ui * hi for ui, hi in zip(u, h)) for u in spins
        )
        formula = min(
            top
            - energies[u]
            + 2
            * sum(
                abs(hi)
                for ui, hi in zip(u, h)
                if ui != sign(hi)
            )
            for u in spins
        )
        assert direct == formula
        checks += 1
    assert top == q * n // 2
    return checks


def verify_row_reduction() -> int:
    checks = 0
    for n in range(2, 11, 2):
        y = tuple([1] * n)
        rows = list(product((-1, 1), repeat=n))
        for same in range(n + 1):
            z = tuple([1] * same + [-1] * (n - same))
            for row in rows:
                S = sum(ri * yi for ri, yi in zip(row, y))
                T = sum(ri * zi for ri, zi in zip(row, z))
                U = sum(row[:same])
                V = sum(row[same:])
                assert (S, T) == (U + V, U - V)
                for target in (-1, 1):
                    direct = sign(S) * sign(T) != target
                    reduced = sign(U + V) * sign(U - V) != target
                    assert direct == reduced
                    checks += 1
    return checks


def response(A: list[list[int]], B: list[list[int]], y: tuple[int, ...]) -> int:
    n = len(A)
    field = matvec(B, y)
    return max(
        quadratic(A, x) + sum(xi * hi for xi, hi in zip(x, field))
        for x in product((-1, 1), repeat=n)
    )


def switch_matrix(A: list[list[int]], s: tuple[int, ...]) -> list[list[int]]:
    n = len(A)
    return [[s[i] * A[i][j] * s[j] for j in range(n)] for i in range(n)]


def verify_query_linked_exposure() -> int:
    _, n, W, _, A = regular_walsh(1)
    spins = list(product((-1, 1), repeat=n))
    top = max(quadratic(A, u) for u in spins)
    u_star = next(u for u in spins if quadratic(A, u) == top)

    # W and two deterministic perturbations exercise flat and nonflat fields.
    bridges = [W]
    for mask in (0x1357, 0x5A3C):
        bridges.append(
            [
                [1 if ((mask >> (i * n + j)) & 1) == 0 else -1 for j in range(n)]
                for i in range(n)
            ]
        )

    checks = 0
    for B in bridges:
        linked: dict[tuple[int, ...], tuple[tuple[int, ...], list[list[int]], int]] = {}
        for y in spins:
            By = matvec(B, y)
            s = tuple(us * sign(value) for us, value in zip(u_star, By))
            child = switch_matrix(A, s)
            value = response(child, B, y)
            assert value == top + sum(abs(value) for value in By)
            linked[y] = (s, child, value)
            checks += 1

        for y in spins:
            sy, Ay, Ryy = linked[y]
            By = matvec(B, y)
            for z in spins:
                if y == z:
                    continue
                sz, Az, Rzz = linked[z]
                Bz = matvec(B, z)
                Rzy = response(Az, B, y)
                Ryz = response(Ay, B, z)

                def deficit(h: tuple[int, ...]) -> int:
                    return top + sum(abs(value) for value in h) - max(
                        quadratic(A, u) + sum(ui * hi for ui, hi in zip(u, h))
                        for u in spins
                    )

                cross_y = tuple(si * value for si, value in zip(sz, By))
                cross_z = tuple(si * value for si, value in zip(sy, Bz))
                delta_y = deficit(cross_y)
                delta_z = deficit(cross_z)
                assert Ryy - Rzy == delta_y
                assert Ryz - Rzz == -delta_z
                assert ((Ryy - Rzy) - (Ryz - Rzz)) >= delta_y + delta_z
                checks += 1
    return checks


def main() -> None:
    child_checks = verify_regular_child(1, exhaustive=True)
    verify_regular_child(2, exhaustive=False)
    deficit_checks = verify_weighted_deficit()
    row_checks = verify_row_reduction()
    exposure_checks = verify_query_linked_exposure()
    print(
        "bounded-cap linear-rate exact checks:",
        {
            "child_spins": child_checks,
            "weighted_deficits": deficit_checks,
            "row_reductions": row_checks,
            "query_linked_exposures": exposure_checks,
            "total": child_checks + deficit_checks + row_checks + exposure_checks,
        },
    )


if __name__ == "__main__":
    main()
