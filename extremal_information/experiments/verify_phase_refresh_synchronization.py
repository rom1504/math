#!/usr/bin/env python3
"""Exact arithmetic checks for phase_refresh_synchronization.md.

This is not a proof search.  It checks the two-state sharpness example, the
Walsh obstruction constant, and finite-state instances of the quantitative
minorization inequality.
"""

from fractions import Fraction
from math import sqrt


def two_state_falsifier(r: int) -> None:
    a = Fraction(1, 2**r)
    # T_0=0, T_1=1; P=(1-a)I+a*uniform.
    average_0 = a * Fraction(1, 2)
    average_1 = (1 - a) + a * Fraction(1, 2)
    assert abs(average_0 - 0) == a / 2
    assert abs(average_1 - 1) == a / 2


def finite_minorization_bound() -> None:
    # K_x = alpha*nu + (1-alpha)*rho_x.  At a maximizing phase, the
    # one-sided response inequality must pay alpha*(M-E_nu phi).
    phi = [Fraction(1, 5), Fraction(3, 5), Fraction(7, 5)]
    nu = [Fraction(1, 6), Fraction(1, 3), Fraction(1, 2)]
    alpha = Fraction(2, 7)
    M = max(phi)
    mean = sum(w * z for w, z in zip(nu, phi))
    required = alpha * (M - mean)
    # Choose the residual measure at the maximizing state, which is the
    # most favorable possible residual for avoiding a defect.
    K_at_max = alpha * mean + (1 - alpha) * M
    assert M - K_at_max == required
    mu = min(nu)
    assert M - mean >= mu * (M - min(phi))


def walsh_constant() -> None:
    c = 89 / (48 * sqrt(3))
    gap = c - 1.01
    assert 1.0705 < c < 1.0706
    assert 0.0605 < gap < 0.0606
    # In Q/n^(3/2) normalization the same obstruction is half as large.
    assert abs(gap / 2 - (89 / (96 * sqrt(3)) - 0.505)) < 1e-15


def main() -> None:
    for r in range(1, 40):
        two_state_falsifier(r)
    finite_minorization_bound()
    walsh_constant()
    print("phase refresh synchronization checks: PASS")
    print(f"Walsh operator-response refresh gap: {89/(48*sqrt(3))-1.01:.12f}")


if __name__ == "__main__":
    main()
