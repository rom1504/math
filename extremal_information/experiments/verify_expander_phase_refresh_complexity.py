#!/usr/bin/env python3
"""Finite checks for expander_phase_refresh_complexity.md.

The script checks the Walsh gauge separation, verifies ER.10 on exhaustive
small reversible kernels/functions, and reproduces the numerical constants
and state-size scaling in ER.19--ER.21.  It is a verifier, not a proof.
"""

from __future__ import annotations

import itertools
import math

import numpy as np


def walsh_gauge_check() -> None:
    H = np.array(
        [
            [1, 1, 1, 1],
            [1, -1, 1, -1],
            [1, 1, -1, -1],
            [1, -1, -1, 1],
        ],
        dtype=float,
    )
    D = np.diag([1, -1, -1, 1])
    gap = np.linalg.norm((H - D @ H @ D) / 2.0, 2)
    assert abs(gap - 2.0) < 1e-12

    # Tensoring by another normalized Hadamard keeps the norm gap exactly two;
    # hollowing cancels because diagonal switching changes no diagonal entry.
    H2 = np.kron(H, H)
    D2 = np.kron(D, np.eye(4))
    A = H2 - np.diag(np.diag(H2))
    A_switch = D2 @ A @ D2
    gap2 = np.linalg.norm((A - A_switch) / 4.0, 2)
    assert abs(gap2 - 2.0) < 1e-12


def exhaustive_er10_check() -> int:
    """Check ER.10 for small symmetric doubly stochastic kernels."""

    kernels = []
    # P_a is a lazy walk on the complete three-vertex graph.  Its mean-zero
    # eigenvalue is 1-3a, so choose a in [0,1/3].
    for a in (0.05, 0.1, 0.2, 0.3, 1.0 / 3.0):
        P = np.full((3, 3), a)
        np.fill_diagonal(P, 1.0 - 2.0 * a)
        rho = np.linalg.norm(P - np.ones((3, 3)) / 3.0, 2)
        kernels.append((P, rho))

    checked = 0
    pi = np.ones(3) / 3.0
    for (P1, rho1), (P2, rho2) in itertools.product(kernels, repeat=2):
        rho = max(rho1, rho2)
        for vals in itertools.product((0.0, 0.25, 0.75, 1.0), repeat=3):
            g = np.array(vals)
            B = 1.0
            product = P1 @ P2
            mean = float(pi @ g)
            rhs_mix = B * rho**2 / math.sqrt(1.0 / 3.0)
            for x in range(3):
                # Set delta_j to the exact positive one-step violations of
                # g <= P_j g + delta_j.  Iteration must imply ER.10.
                d1 = max(0.0, float(np.max(g - P1 @ g)))
                d2 = max(0.0, float(np.max(g - P2 @ g)))
                lhs = float(g[x] - mean)
                rhs = rhs_mix + d1 + d2
                assert lhs <= rhs + 1e-12
                # Also check the direct iterated inequality numerically.
                assert g[x] <= float((product @ g)[x]) + d1 + d2 + 1e-12
                checked += 1
    return checked


def walsh_constants_and_scaling() -> None:
    c_star = 89.0 / (48.0 * math.sqrt(3.0))
    d_star = c_star - 1.01
    assert abs(c_star - 1.07050362412) < 2e-11
    assert abs(d_star - 0.06050362412) < 2e-11

    rho = 0.5
    lam = math.log(1.0 / rho)
    kappa = 0.01
    B = 2.0

    def log_state_lower(delta: float) -> float:
        return (
            math.log(kappa)
            + lam * d_star / delta
            - 2.0 * lam
            - 2.0 * math.log(2.0 * B / d_star)
        )

    # For N=4^r and delta=1/sqrt(N), the asymptotic bit coefficient is D_*.
    rows = []
    for r in (8, 10, 12, 14):
        sqrt_n = 2.0**r
        delta = 1.0 / sqrt_n
        lower_bits = log_state_lower(delta) / math.log(2.0)
        rows.append((r, int(sqrt_n), lower_bits))
    successive_slopes = [
        (rows[i + 1][2] - rows[i][2])
        / (rows[i + 1][1] - rows[i][1])
        for i in range(len(rows) - 1)
    ]
    for slope in successive_slopes:
        assert abs(slope - d_star) < 1e-12

    print("Walsh c_*:", f"{c_star:.12f}")
    print("Walsh D_*:", f"{d_star:.12f}")
    print("rho=1/2, delta=N^-1/2 state-description lower bounds:")
    for r, sqrt_n, bits in rows:
        print(f"  r={r:2d} sqrt(N)={sqrt_n:5d} log2(S)>={bits:10.3f}")


def main() -> None:
    walsh_gauge_check()
    checked = exhaustive_er10_check()
    walsh_constants_and_scaling()
    print(f"PASS: gauge identities and {checked} ER.10 inequalities verified")


if __name__ == "__main__":
    main()

