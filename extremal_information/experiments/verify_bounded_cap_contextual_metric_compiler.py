#!/usr/bin/env python3
"""Exact finite checks for the bounded-cap contextual anti-pin compiler.

The asymptotic code is probabilistic.  This script checks the deterministic
Hadamard, parent-cap, channel, metric, and ownership normalizations on the
explicit Maiorana--McFarland subcodes at n=4 and n=16.
"""

from __future__ import annotations

from itertools import product


def parity(x: int) -> int:
    return bin(x).count("1") & 1


def walsh(d: int) -> list[list[int]]:
    n = 1 << d
    return [[1 if parity(i & j) == 0 else -1 for j in range(n)]
            for i in range(n)]


def mm_vector(m: int, mask: int) -> tuple[int, ...]:
    q = 1 << m
    return tuple(
        1 if (parity(u & v) ^ ((mask >> v) & 1)) == 0 else -1
        for u in range(q) for v in range(q)
    )


def matvec(matrix: list[list[int]], vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(a * b for a, b in zip(row, vector)) for row in matrix)


def quad_hollow(matrix: list[list[int]], x: tuple[int, ...]) -> int:
    return sum(matrix[i][j] * x[i] * x[j]
               for i in range(len(x)) for j in range(i + 1, len(x)))


def build(m: int):
    q = 1 << m
    n = q * q
    W = walsh(2 * m)
    b = mm_vector(m, 0)
    H = [[b[i] * W[i][j] * b[j] for j in range(n)] for i in range(n)]
    A = [[0 if i == j else H[i][j] for j in range(n)] for i in range(n)]
    return q, n, H, A


def clique_energy(total: int, q: int) -> int:
    return (total * total - q) // 2


def spin_mask(spin: tuple[int, ...]) -> int:
    return sum((value == -1) << i for i, value in enumerate(spin))


def energy_table(A: list[list[int]]) -> list[int]:
    n = len(A)
    table = []
    for mask in range(1 << n):
        value = 0
        for i in range(n):
            xi = -1 if (mask >> i) & 1 else 1
            for j in range(i + 1, n):
                xj = -1 if (mask >> j) & 1 else 1
                value += A[i][j] * xi * xj
        table.append(value)
    return table


def parent_cap(
    energies: list[int], n: int, child_mask: int, query_mask: int, q: int
) -> int:
    """Optimize the exact complete parent, enumerating old spins and y-sums."""

    best = 0
    for x_mask in range(1 << n):
        h = energies[x_mask ^ child_mask]
        distance = bin(x_mask ^ query_mask).count("1")
        a = n - 2 * distance
        for total in range(-q, q + 1, 2):
            value = h + a * total + clique_energy(total, q)
            best = max(best, abs(value))
    return best


def original_distance(
    energies: list[int], n: int, s_mask: int, t_mask: int
) -> int:
    return max(abs(energies[x ^ s_mask] - energies[x ^ t_mask])
               for x in range(1 << n))


def verify(m: int) -> int:
    q, n, H, A = build(m)
    energies = energy_table(A)
    one = tuple([1] * n)
    assert matvec(H, one) == tuple([q] * n)
    assert sum(H[i][i] for i in range(n)) == 0
    for i in range(n):
        for j in range(n):
            assert sum(H[i][r] * H[r][j] for r in range(n)) == (
                n if i == j else 0
            )

    # A solution-hidden finite code: scan the explicit MM switches, retaining
    # only the intrinsic two-sided Rayleigh condition.
    code: list[tuple[int, ...]] = []
    for mask in range(1 << q):
        s = mm_vector(m, mask)
        good = True
        for t in code:
            w = tuple(a * b for a, b in zip(s, t))
            rayleigh = sum(a * b for a, b in zip(w, matvec(H, w)))
            if 4 * abs(rayleigh) > q * n:
                good = False
                break
        if good:
            code.append(s)
    assert len(code) >= 2

    # Four representatives suffice at n=16; all are used at n=4.
    code = code if n == 4 else code[:4]
    masks = [spin_mask(s) for s in code]
    responses: dict[tuple[int, int], int] = {}
    diagonal = 3 * q * n // 2 + q * (q - 1) // 2
    off_upper = 11 * q * n // 8 + q * (q - 1) // 2

    for i, s in enumerate(code):
        for j, t in enumerate(code):
            w = tuple(a * b for a, b in zip(s, t))
            rayleigh = sum(a * b for a, b in zip(w, matvec(H, w)))
            if i != j:
                assert 4 * abs(rayleigh) <= q * n
            response = parent_cap(energies, n, masks[i], masks[j], q)
            responses[(i, j)] = response
            if i == j:
                assert response == diagonal
            else:
                assert response <= off_upper

    checks = 0
    for i in range(len(code)):
        for j in range(i + 1, len(code)):
            d0 = original_distance(energies, n, masks[i], masks[j])
            differences = [responses[(i, r)] - responses[(j, r)]
                           for r in range(len(code))]
            dc_twice = max(differences) - min(differences)
            assert dc_twice >= q * n // 4
            assert dc_twice <= 2 * d0
            assert 4 * dc_twice >= d0  # d_C=dc_twice/2 >= d0/8.
            assert d0 >= 3 * q * n // 8
            checks += 1

    print(
        f"m={m}, q={q}, n={n}, audited code={len(code)}, "
        f"parents={len(code) ** 2}, metric pairs={checks}"
    )
    return checks


def main() -> None:
    total = verify(1) + verify(2)
    print(f"bounded-cap contextual compiler checks: PASS ({total} pairs)")


if __name__ == "__main__":
    main()
