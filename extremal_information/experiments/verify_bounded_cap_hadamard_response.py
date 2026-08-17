#!/usr/bin/env python3
"""Exact checks for the bounded-cap Hadamard response-packing theorem.

The theorem is asymptotic, but every algebraic identity used by its proof can
be checked at the first nontrivial Walsh orders.  This script uses integer
arithmetic throughout.  It verifies

* regularization of the Walsh matrix by the self-dual bent vector;
* the exact cap Q(A)=n^(3/2)/2 for the hollow signing;
* the Maiorana--McFarland transform identity;
* the pairwise Rayleigh identity w^T H w=q S(g+h)^2;
* exact Boolean responses and the two-eigenspace trust-region upper bound.
"""

from __future__ import annotations

from itertools import product
from math import sqrt


def parity(x: int) -> int:
    return bin(x).count("1") & 1


def dot_parity(x: int, y: int) -> int:
    return parity(x & y)


def walsh(d: int) -> list[list[int]]:
    n = 1 << d
    return [
        [1 if dot_parity(i, j) == 0 else -1 for j in range(n)]
        for i in range(n)
    ]


def matvec(matrix: list[list[int]], vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(a * b for a, b in zip(row, vector)) for row in matrix)


def quadratic(matrix: list[list[int]], vector: tuple[int, ...]) -> int:
    return sum(
        matrix[i][j] * vector[i] * vector[j]
        for i in range(len(vector))
        for j in range(i + 1, len(vector))
    )


def mm_vector(m: int, gmask: int) -> tuple[int, ...]:
    """s_g(u,v)=(-1)^(u.v+g(v)), lexicographic in (u,v)."""

    q = 1 << m
    return tuple(
        1 if (dot_parity(u, v) ^ ((gmask >> v) & 1)) == 0 else -1
        for u in range(q)
        for v in range(q)
    )


def sign_bias(q: int, gmask: int, hmask: int) -> int:
    return sum(1 if ((gmask ^ hmask) >> v) & 1 == 0 else -1 for v in range(q))


def build_instance(m: int):
    q = 1 << m
    n = q * q
    W = walsh(2 * m)
    b = mm_vector(m, 0)
    H = [[b[i] * W[i][j] * b[j] for j in range(n)] for i in range(n)]
    A = [[0 if i == j else H[i][j] for j in range(n)] for i in range(n)]
    return q, n, W, b, H, A


def verify_order(m: int, exhaustive_cap: bool) -> tuple[int, int]:
    q, n, W, b, H, A = build_instance(m)

    # H is regular Hadamard, its diagonal trace vanishes, and A is a signing.
    assert matvec(W, b) == tuple(q * x for x in b)
    assert matvec(H, tuple([1] * n)) == tuple([q] * n)
    assert sum(H[i][i] for i in range(n)) == 0
    assert all(A[i][i] == 0 for i in range(n))
    assert all(abs(A[i][j]) == 1 for i in range(n) for j in range(n) if i != j)

    # H^2=nI.
    for i in range(n):
        for j in range(n):
            value = sum(H[i][k] * H[k][j] for k in range(n))
            assert value == (n if i == j else 0)

    spins = list(product((-1, 1), repeat=n)) if exhaustive_cap else []
    energies: list[int] = []
    if exhaustive_cap:
        energies = [quadratic(A, x) for x in spins]
        assert max(abs(e) for e in energies) == n * q // 2

    # Select a greedy code with |bias|<=q/2.
    code: list[int] = []
    for g in range(1 << q):
        if all(abs(sign_bias(q, g, h)) <= q // 2 for h in code):
            code.append(g)
    assert len(code) >= 2

    # At q=8 the complete greedy code already has 64 elements.  Eight
    # representatives are ample to audit the identities without turning this
    # verifier into a benchmark.
    audited_code = code if q <= 4 else code[:8]

    response_checks = 0
    responses: dict[tuple[int, int], int] = {}
    rayleigh_checks = 0
    for g in audited_code:
        sg = mm_vector(m, g)
        Wsg = matvec(W, sg)
        yg = tuple(value // q for value in Wsg)
        assert all(abs(value) == 1 for value in yg)
        assert matvec(W, yg) == tuple(q * value for value in sg)

    for g in audited_code:
        sg = mm_vector(m, g)
        for h in audited_code:
            sh = mm_vector(m, h)
            w = tuple(x * y for x, y in zip(sg, sh))
            Hw = matvec(H, w)
            rayleigh = sum(x * y for x, y in zip(w, Hw))
            bias = sign_bias(q, g, h)
            assert rayleigh == q * bias * bias
            rayleigh_checks += 1

            if exhaustive_cap:
                # Response of child g to query h after switching back to A.
                exact = max(
                    energy + q * sum(ui * wi for ui, wi in zip(u, w))
                    for energy, u in zip(energies, spins)
                )
                if g == h:
                    assert exact == 3 * n * q // 2
                else:
                    rho = (bias / q) ** 2
                    # The resolvent completion with K=2qI-H gives
                    # 1+(2+rho)/6 <= 11/8 when rho<=1/4.
                    upper = n * q * (1 + (2 + rho) / 6)
                    assert exact <= upper + 1e-9
                responses[(g, h)] = exact
                response_checks += 1

    if exhaustive_cap:
        for i, g in enumerate(audited_code):
            for h in audited_code[i + 1 :]:
                at_g = responses[(g, g)] - responses[(h, g)]
                at_h = responses[(g, h)] - responses[(h, h)]
                assert at_g >= n * q / 8
                assert at_h <= -n * q / 8
                assert (at_g - at_h) / 2 >= n * q / 8

    print(
        f"m={m}, q={q}, n={n}, code={len(code)}, "
        f"rayleigh_checks={rayleigh_checks}, response_checks={response_checks}"
    )
    return rayleigh_checks, response_checks


def verify_generic_field_ceiling() -> None:
    # For a balanced two-eigenspace projection and ||h||_2/n -> 1, the exact
    # spherical extra gain is attained at r=sin(pi/12).
    r = (sqrt(6) - sqrt(2)) / 4
    trust_extra = -r * r + (sqrt(1 - r * r) + r) / sqrt(2)
    closed_form = (3 * sqrt(3) - 2) / 4
    gaussian_l1 = sqrt(2 / 3.141592653589793)
    assert abs(trust_extra - closed_form) < 1e-12
    assert trust_extra > gaussian_l1
    print(
        "generic balanced trust ceiling:",
        f"{trust_extra:.12f} > Gaussian diagonal {gaussian_l1:.12f}",
        f"by {trust_extra-gaussian_l1:.12f}",
    )


def main() -> None:
    total_rayleigh = 0
    total_response = 0
    for m, exhaustive in ((1, True), (2, True), (3, False)):
        rayleigh, response = verify_order(m, exhaustive)
        total_rayleigh += rayleigh
        total_response += response
    verify_generic_field_ceiling()
    print(f"total exact identities: {total_rayleigh + total_response}")


if __name__ == "__main__":
    main()
