#!/usr/bin/env python3
"""Finite checks for the Wave 39 selector--cut Gram transport reduction.

The concentration inequality itself is analytic (Hanson--Wright plus entropy
duality).  This script checks its exact algebraic inputs and the fractional-
cover information identity on the finite minimizers.
"""

from __future__ import annotations

import itertools
import math
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9


def projective_spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def q_value(a: np.ndarray) -> int:
    return max(abs(int(z @ a @ z)) for z in projective_spins(len(a)))


def gram_mean_check(a: np.ndarray, m: int, name: str) -> None:
    n = len(a)
    a2 = a @ a
    a4 = a2 @ a2
    q = q_value(a)
    selectors = list(itertools.combinations(range(n), m))
    p = Fraction(m, n)
    p2 = Fraction(m * (m - 1), n * (n - 1))
    for z in projective_spins(n):
        row = int(z @ a2 @ z)
        values = []
        for s0 in selectors:
            s = np.asarray(s0, dtype=int)
            y = z[s]
            cols = a[:, s] @ y
            values.append(int(cols @ cols))
        exact_mean = Fraction(sum(values), len(values))
        claimed = p2 * row + (p - p2) * n * (n - 1)
        assert exact_mean == claimed

        h = np.diag(z) @ a2 @ np.diag(z)
        one = np.ones(n, dtype=np.int64)
        assert int(one @ h @ one) == row
        assert int((h @ one) @ (h @ one)) == int(z @ a4 @ z)
        assert int(np.sum(h * h)) == int(np.trace(a4))
        # Exact-minimizer spectral consequence used in the analytic proof.
        # Verify the resulting integer inequalities on the finite test cases.
        assert int(z @ a4 @ z) <= 2 * q * row
        assert int(np.trace(a4)) <= 2 * q * n * (n - 1)
    print(name, "Gram mean/norm identities PASS", "selectors", len(selectors))


def cover_information_check(a: np.ndarray, m: int, name: str) -> None:
    """Use unit weights on every oriented state as a toy feasible cover."""
    n = len(a)
    q = q_value(a)
    selectors = list(itertools.combinations(range(n), m))
    states = []
    for x in projective_spins(n):
        raw = int(x @ a @ x)
        for sigma in (-1, 1):
            states.append((x, sigma, sigma * raw))
    inc = np.zeros((len(selectors), len(states)), dtype=bool)
    for si, s0 in enumerate(selectors):
        s = np.asarray(s0, dtype=int)
        child = a[np.ix_(s, s)]
        for di, (x, sigma, energy) in enumerate(states):
            e = sigma * int(x[s] @ child @ x[s])
            inc[si, di] = 2 * e - energy >= q
    assert np.all(np.sum(inc, axis=1) > 0)

    # P(S) uniform and P(D|S) uniform among its incident unit-weight states.
    joint = inc.astype(float)
    joint /= np.sum(joint, axis=1)[:, None]
    joint /= len(selectors)
    ps = np.sum(joint, axis=1)
    pd = np.sum(joint, axis=0)
    assert np.max(np.abs(ps - 1 / len(selectors))) < 1e-14

    mutual = 0.0
    conditional_selector_kl = 0.0
    prior_divergence = 0.0
    prior = np.full(len(states), 1 / len(states))
    kernel_to_prior = 0.0
    for si in range(len(selectors)):
        for di in np.flatnonzero(inc[si]):
            mass = joint[si, di]
            mutual += mass * math.log(mass / (ps[si] * pd[di]))
            kernel = mass / ps[si]
            kernel_to_prior += mass * math.log(kernel / prior[di])
    for di in np.flatnonzero(pd > 0):
        prior_divergence += pd[di] * math.log(pd[di] / prior[di])
        for si in np.flatnonzero(inc[:, di]):
            cond = joint[si, di] / pd[di]
            conditional_selector_kl += joint[si, di] * math.log(
                cond * len(selectors)
            )

    assert abs(conditional_selector_kl - mutual) < 2e-12
    assert abs(kernel_to_prior - mutual - prior_divergence) < 2e-12
    # Unit weights have W=#states and Z_S>=1, hence the general cover bound.
    assert mutual <= math.log(len(states)) + 1e-12
    print(
        name,
        "cover information PASS",
        "I(S;D)", mutual,
        "kernel-to-prior", kernel_to_prior,
    )


def main() -> None:
    gram_mean_check(A6, 5, "A6,m=5")
    gram_mean_check(A8, 6, "A8,m=6")
    gram_mean_check(A9, 7, "A9,m=7")
    cover_information_check(A6, 5, "A6,m=5")
    print("PASS internal_gram_transport_r39_check")


if __name__ == "__main__":
    main()
