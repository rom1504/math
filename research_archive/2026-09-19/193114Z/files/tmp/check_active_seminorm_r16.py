#!/usr/bin/env python3
"""Exact finite checks for the Wave-16 active-state response functional.

The normalization is the ledger normalization

    Q(M) = max_z |z.T @ M @ z|,

for symmetric zero-diagonal matrices.  We enumerate projective Boolean states
by fixing z[0] = 1, and enumerate both orientations sigma = +/-1.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import product

import numpy as np


A9 = np.array(
    [
        [0, 1, -1, 1, 1, 1, 1, -1, -1],
        [1, 0, 1, -1, -1, -1, 1, -1, -1],
        [-1, 1, 0, 1, 1, 1, 1, 1, -1],
        [1, -1, 1, 0, 1, 1, 1, -1, 1],
        [1, -1, 1, 1, 0, 1, -1, 1, -1],
        [1, -1, 1, 1, 1, 0, -1, -1, -1],
        [1, 1, 1, 1, -1, -1, 0, 1, 1],
        [-1, -1, 1, -1, 1, -1, 1, 0, -1],
        [-1, -1, -1, 1, -1, -1, 1, -1, 0],
    ],
    dtype=np.int64,
)


@lru_cache(None)
def projective_spins(n: int) -> np.ndarray:
    return np.array(
        [(1,) + tail for tail in product((-1, 1), repeat=n - 1)],
        dtype=np.int64,
    )


def energies(M: np.ndarray, spins: np.ndarray | None = None) -> np.ndarray:
    if spins is None:
        spins = projective_spins(len(M))
    return np.einsum("bi,ij,bj->b", spins, M, spins)


def qnorm(M: np.ndarray) -> int:
    return int(np.max(np.abs(energies(M))))


def oriented_data(C: np.ndarray, K: np.ndarray):
    """Return (Q(C), deficits, oriented K energies, exact increment)."""
    spins = projective_spins(len(C))
    c = energies(C, spins)
    k = energies(K, spins)
    signed_c = np.concatenate((c, -c))
    signed_k = np.concatenate((k, -k))
    qc = int(np.max(signed_c))
    deficit = qc - signed_c
    increment = int(np.max(signed_k - deficit))
    assert increment == qnorm(C + K) - qc
    return qc, deficit, signed_k, increment


def active_profile(C: np.ndarray, K: np.ndarray, t: int):
    """Compute h_t, r_t, and the two asymmetric active gauges."""
    qc, deficit, signed_k, increment = oriented_data(C, K)
    active = deficit <= t
    assert np.any(active)
    h = int(np.max(signed_k[active]))
    r = int(np.max(signed_k[active] - deficit[active]))
    plus = max(0, h)
    reverse = max(0, int(np.max(-signed_k[active])))
    return {
        "qc": qc,
        "increment": increment,
        "count": int(np.sum(active)),
        "h": h,
        "r": r,
        "plus": plus,
        "reverse": reverse,
    }


def check_sandwich(C: np.ndarray, K: np.ndarray):
    """Check all three deterministic sandwiches at every deficit breakpoint."""
    qc, deficit, signed_k, increment = oriented_data(C, K)
    qk = qnorm(K)
    # Breakpoints suffice, but include neighboring integers to test genuine
    # strict-tail intervals as well.
    tests = {0, int(np.max(deficit)) + 1}
    for value in map(int, np.unique(deficit)):
        tests.update((max(0, value - 1), value, value + 1))
    for t in sorted(tests):
        active = deficit <= t
        h = int(np.max(signed_k[active]))
        r = int(np.max(signed_k[active] - deficit[active]))
        plus = max(0, h)
        reverse = max(0, int(np.max(-signed_k[active])))

        # Exact deficit-corrected active/tail sandwich.
        assert r <= increment <= max(r, qk - t), (t, r, increment, qk)
        # Signed active-support corollary.
        assert h - t <= increment <= max(h, qk - t), (t, h, increment, qk)
        # Nonnegative asymmetric-seminorm corollary.
        assert -reverse - t <= increment <= max(plus, qk - t)
    return qc, qk, increment


def gauge_signing(n: int, mask: int) -> np.ndarray:
    """First-row-positive representative of one switching class."""
    G = np.zeros((n, n), dtype=np.int64)
    G[0, 1:] = G[1:, 0] = 1
    bit = 0
    for i in range(1, n):
        for j in range(i + 1, n):
            G[i, j] = G[j, i] = 1 if (mask >> bit) & 1 else -1
            bit += 1
    return G


@lru_cache(None)
def all_minimizers(n: int):
    """Enumerate every labelled minimizer, using gauge reps and switchings."""
    best = n * (n - 1)
    representatives = []
    free_edges = (n - 1) * (n - 2) // 2
    for mask in range(1 << free_edges):
        G = gauge_signing(n, mask)
        value = qnorm(G)
        if value < best:
            best, representatives = value, [G]
        elif value == best:
            representatives.append(G)

    labelled = {}
    for G in representatives:
        for tail in product((-1, 1), repeat=n - 1):
            switch = np.array((1,) + tail, dtype=np.int64)
            H = G * np.outer(switch, switch)
            labelled[H.tobytes()] = H
    return best, tuple(labelled.values())


def split_matrix(A: np.ndarray, blocks):
    C = A.copy()
    D = np.zeros_like(A)
    for block in blocks:
        ix = np.ix_(block, block)
        D[ix] = A[ix]
        C[ix] = 0
    return C, D


def oriented_energies(M: np.ndarray) -> np.ndarray:
    e = energies(M)
    return np.concatenate((e, -e))


def check_a9_capture_wall():
    blocks = ((0, 1, 2, 4, 5, 6), (3,), (7,), (8,))
    big = blocks[0]
    C, D = split_matrix(A9, blocks)

    assert qnorm(A9) == 24
    assert qnorm(C) == 22
    assert qnorm(D) == 22
    assert qnorm(A9[np.ix_(big, big)]) == 22
    assert check_sandwich(C, D) == (22, 22, 2)

    expected_D = {
        0: {"count": 3, "h": 2, "r": 2, "plus": 2, "reverse": 6},
        4: {"count": 19, "h": 6, "r": 2, "plus": 6, "reverse": 10},
        8: {"count": 55, "h": 10, "r": 2, "plus": 10, "reverse": 14},
    }
    for t, expected in expected_D.items():
        got = active_profile(C, D, t)
        for key, value in expected.items():
            assert got[key] == value, (t, key, got, expected)

    q6, minimizers = all_minimizers(6)
    assert q6 == 10 and len(minimizers) == 384
    optimal = []
    for G6 in minimizers:
        G = np.zeros_like(A9)
        G[np.ix_(big, big)] = G6
        check_sandwich(C, G)
        if qnorm(C + G) == 24:
            optimal.append(G)
    assert len(optimal) == 40

    persistent_counts = []
    persistent_profiles = set()
    signed_A = oriented_energies(A9)
    signed_C = oriented_energies(C)
    signed_D = oriented_energies(D)
    for G in optimal:
        # Active C states alone miss the positive finite response completely.
        p0 = active_profile(C, G, 0)
        assert p0 == {
            "qc": 22,
            "increment": 2,
            "count": 3,
            "h": -6,
            "r": -6,
            "plus": 0,
            "reverse": 6,
        }
        # Deficit 4 first sees the actual +2 response.  Deficit 8 is the first
        # breakpoint at which r_t and the universal tail Q(G)-t both equal 2.
        p4 = active_profile(C, G, 4)
        p8 = active_profile(C, G, 8)
        assert (p4["h"], p4["r"], p4["plus"], p4["reverse"]) == (6, 2, 6, 10)
        assert (p8["h"], p8["r"], p8["plus"], p8["reverse"]) == (10, 2, 10, 10)
        assert qnorm(G) - 8 == p8["r"] == 2

        # Audit the persistent-ground interpolation wall independently.
        H = C + G
        signed_H = oriented_energies(H)
        signed_G = oriented_energies(G)
        common = (signed_A == 24) & (signed_H == 24)
        count = int(np.sum(common))
        assert 4 <= count <= 13
        persistent_counts.append(count)
        for values in zip(signed_C[common], signed_D[common], signed_G[common]):
            persistent_profiles.add(tuple(map(int, values)))
        # Convexity plus a common ground makes the whole interpolation flat;
        # integer endpoint/midpoint checks guard the explicit implementation.
        assert qnorm(2 * C + D + G) == 48

    assert persistent_profiles == {(14, 10, 10), (18, 6, 6)}
    print("A9 6+1+1+1:")
    print("  Q(C),Q(D),Delta_C(D) = 22,22,2")
    print("  optimal q6 completions = 40")
    print("  every optimal G: (h0,r0,a0+,a0-) = (-6,-6,0,6)")
    print("  every optimal G: r4=2; r8=Q(G)-8=2")
    print(
        "  persistent-ground counts min/max =",
        min(persistent_counts),
        max(persistent_counts),
    )


def check_a9_three_plus_six():
    blocks = ((0, 1, 2), (3, 4, 5, 6, 7, 8))
    C, D = split_matrix(A9, blocks)
    assert qnorm(C) == 24
    assert qnorm(D) == 20
    assert check_sandwich(C, D) == (24, 20, 0)
    p0 = active_profile(C, D, 0)
    assert (p0["count"], p0["h"], p0["r"], p0["plus"], p0["reverse"]) == (
        2,
        0,
        0,
        0,
        0,
    )

    # The explicit best replacement from (10.553).
    G3 = np.array([[0, -1, -1], [-1, 0, -1], [-1, -1, 0]], dtype=np.int64)
    G6 = np.array(
        [
            [0, -1, -1, 1, -1, 1],
            [-1, 0, -1, 1, 1, -1],
            [-1, -1, 0, -1, -1, -1],
            [1, 1, -1, 0, -1, -1],
            [-1, 1, -1, -1, 0, 1],
            [1, -1, -1, -1, 1, 0],
        ],
        dtype=np.int64,
    )
    G = np.zeros_like(A9)
    G[np.ix_(blocks[0], blocks[0])] = G3
    G[np.ix_(blocks[1], blocks[1])] = G6
    assert qnorm(G3) == 6 and qnorm(G6) == 10 and qnorm(G) == 16
    assert check_sandwich(C, G) == (24, 16, 4)
    p0 = active_profile(C, G, 0)
    assert (p0["h"], p0["r"], p0["plus"], p0["reverse"]) == (4, 4, 4, 4)
    print("A9 3+6: original Delta=0; explicit best replacement Delta=4=h0")


def check_all_positive(n: int = 12, s: int = 3):
    """Audit the adaptive witness (10.572) without floating point.

    F=(n-1)E=(n-1)D-(s-1)A is the integer-scaled centered mask.
    """
    assert n % s == 0 and (n // s) % 2 == 0 and s < n // 2
    A = np.ones((n, n), dtype=np.int64) - np.eye(n, dtype=np.int64)
    blocks = tuple(tuple(range(j, j + s)) for j in range(0, n, s))
    C, D = split_matrix(A, blocks)
    F = (n - 1) * D - (s - 1) * A
    scaled_C = (n - 1) * C

    assert qnorm(C) == n * (n - s)
    assert qnorm(F) == n * n * (s - 1)
    # C+E=(1-p)A.  Scaling by n-1 makes this an integer identity.
    assert np.array_equal(scaled_C + F, (n - s) * A)
    assert check_sandwich(scaled_C, F) == (
        (n - 1) * n * (n - s),
        n * n * (s - 1),
        0,
    )
    p0 = active_profile(scaled_C, F, 0)
    assert p0["count"] == 1
    assert (p0["h"], p0["r"], p0["plus"], p0["reverse"]) == (0, 0, 0, 0)

    # Half the blocks have each sign.  This is exactly the adaptive state in
    # (10.572), scaled by n-1.
    z = np.empty(n, dtype=np.int64)
    for j, block in enumerate(blocks):
        z[list(block)] = 1 if j < len(blocks) // 2 else -1
    c_witness = int(z @ C @ z)
    f_witness = int(z @ F @ z)
    assert c_witness == -n * s
    assert f_witness == n * n * (s - 1)
    # With the orientation making F positive, the C-deficit is n^2 before
    # scaling and is therefore far outside every subleading active margin.
    assert qnorm(C) - c_witness == n * n

    internal_sum = sum(qnorm(A[np.ix_(block, block)]) for block in blocks)
    assert internal_sum == n * (s - 1)
    print("all-positive centered witness:")
    print(f"  n={n}, s={s}, sum_i Q(D_i)={internal_sum}")
    print("  a_{C,0}^+(F)=a_{C,0}^+(-F)=0")
    print(f"  Q(F)={qnorm(F)} and Delta_C(F)=0 (integer scale n-1)")
    print(f"  adaptive witness F-energy={f_witness}, unscaled={f_witness}/{n-1}")


def main():
    check_a9_capture_wall()
    check_a9_three_plus_six()
    check_all_positive()
    print("PASS: every exact active-state identity and finite wall check succeeded")


if __name__ == "__main__":
    main()
