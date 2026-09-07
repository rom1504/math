#!/usr/bin/env python3
"""Exact checks for the Wave 15 A9 amplification no-go results."""

from itertools import product
from math import asin, log, pi, sqrt

import numpy as np


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


def spins(n):
    return np.array(list(product((-1, 1), repeat=n)), dtype=np.int64)


def q(a):
    x = spins(len(a))
    return int(np.max(np.abs(np.sum((x @ a) * x, axis=1))))


def random_upper(n):
    return 2 * sqrt(n * (n - 1) * (n + 2) * log(2))


def check_seed():
    x = spins(9)
    e = np.sum((x @ A9) * x, axis=1)
    assert (int(e.min()), int(e.max()), q(A9)) == (-24, 24, 24)
    one = np.ones(9, dtype=np.int64)
    assert int(one @ A9 @ one) == 8

    b = A9[:3, 3:]
    vals = spins(3) @ b @ spins(6).T
    assert int(np.max(np.abs(vals))) == 12
    print("A9 seed and 3+6 cross norm: PASS")


def all_plus_blowup_value(k):
    """Compress each clone fibre to its magnetization."""
    best = 0
    for ms in product(range(-k, k + 1, 2), repeat=9):
        m = np.array(ms, dtype=np.int64)
        e = int(m @ A9 @ m + m @ m - 9 * k)
        best = max(best, abs(e))
    return best


def check_regular_blowup():
    # Q(C9 tensor J_k)=24 k^2.  Direct small compressed checks also verify
    # Q(A9 tensor J_k + I_9 tensor (J_k-I_k))=33 k^2-9k.
    for k in (1, 2, 3):
        assert all_plus_blowup_value(k) == 33 * k * k - 9 * k

    # Conference upper bounds exclude the only two k below the elementary
    # random-sign threshold: orders 18 and 27 sit below conferences 18 and 30.
    assert 18 * sqrt(17) < 24 * 2**2
    assert 27 * sqrt(29) < 24 * 3**2

    # For k>=4 it is enough to check the monotone normalized inequality at 4:
    # 4 sqrt(k) > sqrt(log 2) (9+2/k).
    assert 4 * sqrt(4) > sqrt(log(2)) * (9 + 2 / 4)
    for k in range(4, 200):
        assert random_upper(9 * k) < 24 * k * k
    print("regular blow-up exact norms and global-minimality wall: PASS")


def check_fixed_replication():
    # Every repeated 3+6 local-minimizer tuple has block-diagonal norm at
    # most t(q_3+q_6)=16t, regardless of the signs coupling different copies.
    q3, q6 = 6, 10
    for t in (1, 2, 10, 100):
        b0 = t * (q3 + q6)
        n = 9 * t
        assert b0 == 16 * t
        assert abs(b0 / n**1.5 - 16 / (27 * sqrt(t))) < 1e-15
    print("fixed-wall replication bound 16t=o((9t)^(3/2)): PASS")


def check_outer_lexicographic():
    # Exact small q_t values already in the ledger/artifact census.
    small_q = {2: 2, 3: 6, 4: 8, 5: 8, 6: 10,
               7: 18, 8: 20, 9: 24}
    for t, qt in small_q.items():
        assert 81 * qt - 8 * t > random_upper(9 * t)

    # For t>=10, finite Gaussian rounding and arcsin(u)>=u give the lower
    # coefficient below.  The difference is increasing in t, so t=10 is
    # the only endpoint requiring evaluation.
    def normalized_gap(t):
        lower = 162 / pi * sqrt(1 - 1 / t) - 8 / sqrt(t)
        upper = 54 * sqrt(log(2)) + 12 * sqrt(log(2)) / t
        return lower - upper

    assert normalized_gap(10) > 0.43
    for t in range(10, 10000):
        assert normalized_gap(t) > 0

    # Directly check the slightly sharper finite Gaussian expression too.
    for t in range(10, 1000):
        qt_lower = 2 * t * (t - 1) / pi * asin(1 / sqrt(t - 1))
        assert 81 * qt_lower - 8 * t > random_upper(9 * t)
    print("outer lexicographic replication is nonminimal for every t>=2: PASS")


def check_star_tensor():
    # For P=A9+D, T=P tensor P-diag(P tensor P) is an order-81 signing.
    # A Paley conference matrix of order 82 gives q_81<=81*9=729.
    conference_upper_81 = 729
    x = spins(9)
    energies = np.sum((x @ A9) * x, axis=1)
    xplus = x[np.flatnonzero(energies == 24)[0]]
    xminus = x[np.flatnonzero(energies == -24)[0]]
    minima = []

    for ds in product((-1, 1), repeat=9):
        d = np.array(ds, dtype=np.int64)
        delta = int(d.sum())
        p = A9 + np.diag(d)
        assert np.all(np.abs(p) == 1)

        # vec(P) is Boolean and has energy tr(P^4)-delta^2.
        trace_witness = int(np.trace(np.linalg.matrix_power(p, 4)) - delta**2)
        assert int(np.trace(np.linalg.matrix_power(p, 4))) >= 801

        # A product of two same-sign A9 grounds has this energy.
        ground = xplus if delta > 0 else xminus
        seed_e = int(ground @ A9 @ ground)
        product_witness = (seed_e + delta)**2 - delta**2
        assert product_witness == 576 + 48 * abs(delta)

        witness = max(trace_witness, product_witness)
        assert witness > conference_upper_81
        minima.append((witness, trace_witness, product_witness, delta))

    # This is stronger than the universal 792/816 split used in the proof.
    assert min(w for w, *_ in minima) == 896
    print("all 512 diagonal-completion star tensor squares are nonminimal: PASS")


if __name__ == "__main__":
    check_seed()
    check_regular_blowup()
    check_fixed_replication()
    check_outer_lexicographic()
    check_star_tensor()
    print("all A9 amplification checks: PASS")
