#!/usr/bin/env python3
"""Exact finite active-state audit for the A9 6+1+1+1 capture wall.

The checker uses only integer Boolean enumeration.  Radical signs in the
temporal comparisons are certified by explicit integer square inequalities.
Spins are projectivized by fixing their first coordinate to +1.  An oriented
state is a pair (sigma,z), sigma in {+1,-1}; it is active for M when
sigma*z^T M z = Q(M).
"""

from __future__ import annotations

import argparse
from collections import Counter
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

A8 = np.array(
    [
        [0, 1, 1, 1, 1, 1, 1, 1],
        [1, 0, 1, -1, 1, 1, -1, -1],
        [1, 1, 0, 1, -1, 1, -1, -1],
        [1, -1, 1, 0, -1, -1, -1, 1],
        [1, 1, -1, -1, 0, -1, 1, -1],
        [1, 1, 1, -1, -1, 0, 1, 1],
        [1, -1, -1, -1, 1, 1, 0, 1],
        [1, -1, -1, 1, -1, 1, 1, 0],
    ],
    dtype=np.int64,
)

A7 = np.array(
    [
        [0, 1, 1, 1, 1, 1, 1],
        [1, 0, 1, 1, -1, -1, 1],
        [1, 1, 0, 1, -1, 1, -1],
        [1, 1, 1, 0, 1, -1, -1],
        [1, -1, -1, 1, 0, -1, -1],
        [1, -1, 1, -1, -1, 0, -1],
        [1, 1, -1, -1, -1, -1, 0],
    ],
    dtype=np.int64,
)

A5 = np.array(
    [
        [0, -1, 1, -1, 1],
        [-1, 0, -1, 1, 1],
        [1, -1, 0, 1, 1],
        [-1, 1, 1, 0, 1],
        [1, 1, 1, 1, 0],
    ],
    dtype=np.int64,
)

# The wall partition from (10.579), in zero-based indexing.
V = (0, 1, 2, 4, 5, 6)
SINGLETONS = (3, 7, 8)
QMIN = {1: 0, 2: 2, 3: 6, 4: 8, 5: 8, 6: 10, 7: 18, 8: 20, 9: 24}

# Raw upper-triangle sign codes of the forty labelled G6 completions.  Bit k
# is one exactly when edge k in lexicographic (i,j), i<j, order has sign +1.
EXPECTED_CODES = (
    0x069A, 0x0A4E, 0x0C8E, 0x129C, 0x1A07, 0x1B9F, 0x1C0B, 0x1E48,
    0x1E84, 0x222E, 0x2607, 0x279F, 0x2CDF, 0x2EFA, 0x2FAE, 0x308E,
    0x340D, 0x3628, 0x3682, 0x38BF, 0x3AFC, 0x3BCE, 0x3EB3, 0x3ED5,
    0x3F2B, 0x3F4D, 0x460B, 0x4ADF, 0x520D, 0x5ED9, 0x5F8D, 0x62BF,
    0x6E6B, 0x6EA7, 0x76B9, 0x778B, 0x7A6D, 0x7AC7, 0x7CAD, 0x7CCB,
)


@lru_cache(None)
def spins(n: int) -> np.ndarray:
    return np.array([(1,) + t for t in product((-1, 1), repeat=n - 1)], dtype=np.int64)


def energies(M: np.ndarray) -> np.ndarray:
    Z = spins(len(M))
    return np.einsum("bi,ij,bj->b", Z, M, Z)


def qnorm(M: np.ndarray) -> int:
    return int(np.max(np.abs(energies(M))))


def edge_code(M: np.ndarray) -> int:
    bit = 0
    code = 0
    for i in range(len(M)):
        for j in range(i + 1, len(M)):
            if int(M[i, j]) == 1:
                code |= 1 << bit
            bit += 1
    return code


def grounds(M: np.ndarray) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Projective positive and negative grounds; both must attain Q(M)."""
    q = qnorm(M)
    e = energies(M)
    positive = [z.copy() for z, value in zip(spins(len(M)), e) if int(value) == q]
    negative = [z.copy() for z, value in zip(spins(len(M)), e) if int(value) == -q]
    return positive, negative


def max_balanced_root_residual(M: np.ndarray) -> int:
    """At P=N=Q, (10.556) says root residual is |p^T M n|."""
    positive, negative = grounds(M)
    assert positive and negative
    return max(abs(int(p @ M @ n)) for p in positive for n in negative)


def gauge_signing(n: int, mask: int) -> np.ndarray:
    """First-row-positive representative of a switching class."""
    G = np.zeros((n, n), dtype=np.int64)
    G[0, 1:] = G[1:, 0] = 1
    bit = 0
    for i in range(1, n):
        for j in range(i + 1, n):
            G[i, j] = G[j, i] = 1 if (mask >> bit) & 1 else -1
            bit += 1
    return G


@lru_cache(None)
def all_minimizers(n: int) -> tuple[int, tuple[np.ndarray, ...]]:
    """All labelled minimizers: gauge reps followed by all switchings."""
    best = n * (n - 1)
    reps: list[np.ndarray] = []
    free_edges = (n - 1) * (n - 2) // 2
    for mask in range(1 << free_edges):
        G = gauge_signing(n, mask)
        value = qnorm(G)
        if value < best:
            best, reps = value, [G]
        elif value == best:
            reps.append(G)
    labelled: dict[bytes, np.ndarray] = {}
    for G in reps:
        for tail in product((-1, 1), repeat=n - 1):
            s = np.array((1,) + tail, dtype=np.int64)
            H = G * np.outer(s, s)
            labelled[H.tobytes()] = H
    return best, tuple(labelled.values())


def wall_matrices() -> tuple[np.ndarray, np.ndarray]:
    D = A9[np.ix_(V, V)].copy()
    C = A9.copy()
    C[np.ix_(V, V)] = 0
    # The other three diagonal blocks are singletons, hence already zero.
    return C, D


def hybrid(C: np.ndarray, G: np.ndarray) -> np.ndarray:
    H = C.copy()
    H[np.ix_(V, V)] = G
    return H


def active_records(C: np.ndarray, D: np.ndarray, G: np.ndarray):
    """Return all oriented active states and their exact energy profile."""
    H = hybrid(C, G)
    Z = spins(9)
    eh = energies(H)
    ec = energies(C)
    ea = energies(A9)
    out = []
    for i, z in enumerate(Z):
        y = z[list(V)]
        gd = int(y @ G @ y)
        dd = int(y @ D @ y)
        for sigma in (-1, 1):
            if sigma * int(eh[i]) != 24:
                continue
            # c,g,d denote oriented cross, replacement-block, original-block
            # energies.  delta_G and delta_D are the corresponding deficits.
            c = sigma * int(ec[i])
            g = sigma * gd
            d = sigma * dd
            assert c + g == 24
            assert sigma * int(ea[i]) == c + d
            out.append(
                {
                    "sigma": sigma,
                    "z": tuple(map(int, z)),
                    "c": c,
                    "g": g,
                    "d": d,
                    "delta_g": 10 - g,
                    "delta_d": 22 - d,
                    "old_total": c + d,
                    "point_gain": (22 - d) - (10 - g),
                }
            )
    return out


def projective(x: np.ndarray) -> tuple[int, ...]:
    if int(x[0]) < 0:
        x = -x
    return tuple(map(int, x))


def endpoint_records(H: np.ndarray):
    """All projective positive/negative ground pairs and both retained shores."""
    Z = spins(9)
    e = energies(H)
    positives = [z.copy() for z, value in zip(Z, e) if int(value) == 24]
    negatives = [z.copy() for z, value in zip(Z, e) if int(value) == -24]
    records = []
    for p in positives:
        for n in negatives:
            agree = tuple(i for i in range(9) if int(p[i] * n[i]) == 1)
            disagree = tuple(i for i in range(9) if int(p[i] * n[i]) == -1)
            b = int(p @ H @ n)
            # H is balanced (P=N=24), so total root residual is |p^T H n|.
            for X, shore in ((agree, "agree"), (disagree, "disagree")):
                Y = tuple(i for i in range(9) if i not in X)
                HX = H[np.ix_(X, X)]
                qx = qnorm(HX) if X else 0
                # In the p-gauge, h_X is the energy of p restricted to X.
                px = p[list(X)] if X else np.empty(0, dtype=np.int64)
                hx = int(px @ HX @ px) if X else 0
                residual_sum = abs(hx)  # I(H)=0 in (10.534).
                records.append(
                    {
                        "p": projective(p),
                        "n": projective(n),
                        "cross_gram": b,
                        "shore": shore,
                        "X": X,
                        "Y": Y,
                        "m": len(X),
                        "qx": qx,
                        "decrement": 24 - qx,
                        "h": hx,
                        "residual_sum": residual_sum,
                        # Sign of centered demand x = Q(H[X])-(8/9)|X|^(3/2)
                        # is certified without floating point by 81*qx^2 ? 64*m^3.
                        "temporal_sign_poly": 81 * qx * qx - 64 * len(X) ** 3,
                    }
                )
    return positives, negatives, records


def active_accessibility(H: np.ndarray, records: list[dict]) -> Counter:
    """Max endpoint residual available while keeping each active state fixed."""
    positive, negative = grounds(H)
    out = Counter()
    for row in records:
        z = np.array(row["z"], dtype=np.int64)
        opposite = negative if row["sigma"] == 1 else positive
        accessible = max(abs(int(z @ H @ w)) for w in opposite)
        out[accessible] += 1
    return out


def audit_a8_response_wall() -> dict:
    """Falsify Delta_response <= original endpoint residual on exact A8."""
    blocks = ((0, 3, 4, 5), (1, 2, 6, 7))
    assert qnorm(A8) == 20
    block_norms = tuple(qnorm(A8[np.ix_(block, block)]) for block in blocks)
    assert block_norms == (12, 12)
    C = A8.copy()
    for block in blocks:
        C[np.ix_(block, block)] = 0
    assert qnorm(C) == 16

    q4, minimizers4 = all_minimizers(4)
    assert q4 == 8 and len(minimizers4) == 48
    best = 10**9
    best_hybrids: list[np.ndarray] = []
    for G0 in minimizers4:
        for G1 in minimizers4:
            H = C.copy()
            H[np.ix_(blocks[0], blocks[0])] = G0
            H[np.ix_(blocks[1], blocks[1])] = G1
            value = qnorm(H)
            if value < best:
                best, best_hybrids = value, [H]
            elif value == best:
                best_hybrids.append(H)
    assert best == 20 and len(best_hybrids) == 16

    positive, negative = grounds(A8)
    assert (len(positive), len(negative)) == (4, 4)
    cross_grams = {int(p @ A8 @ n) for p in positive for n in negative}
    assert cross_grams == {0}
    for p in positive:
        for n in negative:
            shores = (
                tuple(i for i in range(8) if int(p[i] * n[i]) == 1),
                tuple(i for i in range(8) if int(p[i] * n[i]) == -1),
            )
            for X in shores:
                assert len(X) == 4
                MX = A8[np.ix_(X, X)]
                px = p[list(X)]
                assert qnorm(MX) == 8 and int(px @ MX @ px) == 0

    # Original response Phi_D=20-(12+12)=-4; best response
    # Phi_G=20-(8+8)=4.  Its rise is eight, exactly the captured internal
    # excess of these response blocks; they are not endpoint shores.
    phi_original = 20 - sum(block_norms)
    phi_best = best - 2 * q4
    response_rise = phi_best - phi_original
    captured_excess = sum(value - q4 for value in block_norms)
    original_residual = max_balanced_root_residual(A8)
    assert (phi_original, phi_best, response_rise, captured_excess, original_residual) == (-4, 4, 8, 8, 0)
    assert response_rise > original_residual  # exact counterexample

    # Every best hybrid has gained residual 8 and shares four oriented grounds
    # with A8.  Those are replacement-matrix resources, not A8 resources.
    hybrid_residual_hist = Counter()
    persistent_hist = Counter()
    ea = energies(A8)
    for H in best_hybrids:
        hybrid_residual_hist[max_balanced_root_residual(H)] += 1
        eh = energies(H)
        persistent = sum(
            1
            for vh, va in zip(eh, ea)
            for sigma in (-1, 1)
            if sigma * int(vh) == 20 and sigma * int(va) == 20
        )
        persistent_hist[persistent] += 1
    assert hybrid_residual_hist == Counter({8: 16})
    assert persistent_hist == Counter({4: 16})

    # Exact centered endpoint demand is 8-5*sqrt(2)>0, since 8^2>50.
    assert 8 * 8 > 25 * 2
    return {
        "block_norms": block_norms,
        "cross_norm": qnorm(C),
        "best": best,
        "best_count": len(best_hybrids),
        "response_rise": response_rise,
        "captured_excess": captured_excess,
        "original_residual": original_residual,
        "hybrid_residual_hist": hybrid_residual_hist,
    }


def audit_timing_walls() -> dict:
    """Recheck the A5 prefix, A7 suffix, and A8 neutrality walls exactly."""
    # A5: the displayed positive/negative pair has zero residual; retaining
    # the order-four shore has Q=8 and centered demand
    # 8-64*sqrt(5)/25>0 (equivalently 200>64*sqrt(5)).
    assert qnorm(A5) == 8
    p5 = np.ones(5, dtype=np.int64)
    n5 = np.array((1, 1, 1, 1, -1), dtype=np.int64)
    assert (int(p5 @ A5 @ p5), int(n5 @ A5 @ n5), int(p5 @ A5 @ n5)) == (8, -8, 0)
    child5 = A5[:4, :4]
    assert qnorm(child5) == 8 and int(p5[:4] @ child5 @ p5[:4]) == 0
    assert 200 * 200 > 64 * 64 * 5

    # A7 deterministic 7->6->5 chain.  Deleting vertex 0 leaves a q6
    # minimizer; its negative ground below has the unique field-5 vertex 5,
    # whose deletion leaves a q5 minimizer.  Both orientations of the second
    # endpoint bucket have zero residual.
    assert qnorm(A7) == 18
    n7 = np.array((1, -1, -1, -1, -1, -1, -1), dtype=np.int64)
    assert int(n7 @ A7 @ n7) == -18
    assert tuple(map(int, -n7 * (A7 @ n7))) == (6, 0, 0, 0, 4, 4, 4)
    child6 = A7[1:, 1:]
    assert qnorm(child6) == 10
    n6 = np.array((1, -1, -1, -1, -1, -1), dtype=np.int64)
    assert int(n6 @ child6 @ n6) == -10
    assert tuple(map(int, -n6 * (child6 @ n6))) == (1, 1, 1, 1, 1, 5)
    child5_from7 = child6[:5, :5]
    assert qnorm(child5_from7) == 8
    y = n6[:5]
    h = int(y @ child5_from7 @ y)
    L = abs(int(child6[5, :5] @ y))
    residuals = tuple(max(2 * L + sigma * h - 10, 0) for sigma in (-1, 1))
    assert (h, L, residuals) == (0, 5, (0, 0))

    # Exact signs from (10.544), proved only with integer squaring:
    # x1=10-108*sqrt(42)/49<0 because 490^2<(108^2)*42;
    # x2=(108*sqrt(42)-90*sqrt(35))/49-2>0.  For x2, square twice:
    # 108*sqrt(42)>98+90*sqrt(35) iff
    # 196784>17640*sqrt(35), and the latter square inequality holds.
    # The full sum 8-90*sqrt(35)/49 is negative.
    assert 490 * 490 < 108 * 108 * 42
    assert 196_784 > 0 and 196_784**2 > 17_640**2 * 35
    assert 392 * 392 < 90 * 90 * 35

    return {
        "a5_residual": 0,
        "a5_demand": "8-64*sqrt(5)/25 > 0",
        "a7_residuals": residuals,
        "a7_x1": "10-108*sqrt(42)/49 < 0",
        "a7_x2": "(108*sqrt(42)-90*sqrt(35))/49-2 > 0",
        "a8_endpoint_demand": "8-5*sqrt(2) > 0",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full",
        action="store_true",
        help="print all forty completion fingerprints and summary columns",
    )
    parser.add_argument(
        "--dump-active",
        action="store_true",
        help="print every oriented active spin and its exact energy/deficit profile",
    )
    args = parser.parse_args()

    C, D = wall_matrices()
    assert qnorm(A9) == QMIN[9] == 24
    assert qnorm(C) == 22 and qnorm(D) == 22
    q6, minimizers = all_minimizers(6)
    assert q6 == 10 and len(minimizers) == 384

    completions = []
    for G in minimizers:
        H = hybrid(C, G)
        if qnorm(H) == 24:
            completions.append(G)
    assert len(completions) == 40
    completions.sort(key=edge_code)
    assert tuple(map(edge_code, completions)) == EXPECTED_CODES

    active_hist = Counter()  # orientation suppressed
    orientation_hist = Counter()
    active_count_hist = Counter()
    endpoint_count_hist = Counter()
    endpoint_hist = Counter()  # cross sign and shore name suppressed
    endpoint_pair_hist = Counter()
    persistent_count_hist = Counter()
    persistent_profile_hist = Counter()
    persistent_orientation_hist = Counter()
    accessible_hist = Counter()
    persistent_accessible_hist = Counter()
    completions_without_persistent_high_pair = 0
    neutral_positive_shores = Counter()
    response_min_old_score = 10**9
    response_max_old_score = -10**9
    completion_rows = []
    active_dump = []

    for G in completions:
        H = hybrid(C, G)
        records = active_records(C, D, G)
        if args.dump_active:
            active_dump.extend((edge_code(G), row) for row in records)
        active_count_hist[len(records)] += 1
        orientation_hist.update(r["sigma"] for r in records)
        active_hist.update(
            (r["c"], r["g"], r["delta_g"], r["d"], r["delta_d"], r["old_total"], r["point_gain"])
            for r in records
        )
        old_scores = [r["old_total"] - 22 for r in records]
        response_min_old_score = min(response_min_old_score, *old_scores)
        response_max_old_score = max(response_max_old_score, *old_scores)

        persistent = [r for r in records if r["old_total"] == 24]
        persistent_count_hist[len(persistent)] += 1
        persistent_orientation_hist.update(r["sigma"] for r in persistent)
        persistent_profile_hist.update((r["c"], r["d"], r["g"]) for r in persistent)

        positive, negative = grounds(H)
        access_values = []
        persistent_access_values = []
        for row in records:
            z = np.array(row["z"], dtype=np.int64)
            opposite = negative if row["sigma"] == 1 else positive
            available = max(abs(int(z @ H @ w)) for w in opposite)
            access_values.append(available)
            accessible_hist[available] += 1
            if row["old_total"] == 24:
                persistent_access_values.append(available)
                persistent_accessible_hist[available] += 1
        assert set(access_values) <= {8, 16}
        # This is the exact active-witness failure: Delta_response=12, yet
        # every completion has active (indeed persistent in every completion)
        # states with at most eight units of accessible endpoint residual.
        assert min(persistent_access_values) == 8 < 12
        if max(persistent_access_values) == 8:
            completions_without_persistent_high_pair += 1

        positives, negatives, endpoints = endpoint_records(H)
        endpoint_count_hist[(len(positives), len(negatives), len(endpoints))] += 1
        pair_values = {}
        for row in endpoints:
            pair_values[(row["p"], row["n"])] = abs(row["cross_gram"])
            endpoint_hist[
                (
                    abs(row["cross_gram"]), row["m"], row["qx"],
                    abs(row["h"]), int(np.sign(row["temporal_sign_poly"])),
                )
            ] += 1
        endpoint_pair_hist.update(pair_values.values())
        high_pairs = sum(value == 16 for value in pair_values.values())
        assert high_pairs in (2, 3)
        neutral_pos = sum(
            row["cross_gram"] == 0 and row["temporal_sign_poly"] > 0
            for row in endpoints
        )
        neutral_positive_shores[neutral_pos] += 1
        # Neither the response block V nor its complement is an endpoint shore.
        assert not any(set(row["X"]) in (set(V), set(SINGLETONS)) for row in endpoints)
        completion_rows.append(
            (
                edge_code(G), len(positives), len(negatives), len(records),
                len(persistent), persistent_access_values.count(8),
                persistent_access_values.count(16), high_pairs,
            )
        )

    # Response facts: Phi_D=24-22=2, Phi_G=24-10=14, Delta=12.
    assert (response_min_old_score, response_max_old_score) == (-22, 2)
    assert active_count_hist == Counter({21: 24, 25: 16})
    assert orientation_hist == Counter({1: 496, -1: 408})
    assert persistent_count_hist == Counter({4: 2, 5: 4, 6: 6, 7: 6, 8: 8, 9: 6, 10: 4, 11: 2, 13: 2})
    assert persistent_profile_hist == Counter({(14, 10, 10): 168, (18, 6, 6): 144})
    assert persistent_orientation_hist == Counter({1: 176, -1: 136})
    assert accessible_hist == Counter({8: 712, 16: 192})
    assert persistent_accessible_hist == Counter({8: 248, 16: 64})
    assert completions_without_persistent_high_pair == 8
    assert endpoint_pair_hist == Counter({0: 1944, 4: 1152, 8: 1896, 16: 96})
    assert neutral_positive_shores == Counter({55: 24, 96: 16})
    assert all(
        r["delta_g"] >= 0 and r["delta_d"] >= 0
        for G in completions
        for r in active_records(C, D, G)
    )

    # Each persistent state gives an identically active line on
    # A_t=C+(1-t)D+tG.  Convexity gives Q(A_t)<=24, so Q(A_t)=24 exactly and
    # Phi_t=24-[(1-t)22+t10]=2+12t.  Integer endpoint equality is all that the
    # finite checker must establish.
    assert all(persistent_count_hist.elements())

    # The direct V-retaining deletion from A9 has centered demand
    # 22-(16/3)sqrt(6)>0, but terminal excess 22-q6=12; after paying it the
    # rooted total is 10-(16/3)sqrt(6)<0.  Both signs are exact by squaring.
    assert 66 * 66 > 16 * 16 * 6
    assert 30 * 30 < 16 * 16 * 6

    original_p, original_n, original_endpoints = endpoint_records(A9)
    assert (len(original_p), len(original_n)) == (10, 15)
    original_pair_values = {
        (row["p"], row["n"]): abs(row["cross_gram"])
        for row in original_endpoints
    }
    assert Counter(original_pair_values.values()) == Counter({0: 56, 4: 36, 8: 46, 12: 4, 16: 8})
    assert not any(set(row["X"]) in (set(V), set(SINGLETONS)) for row in original_endpoints)
    original_a9_residual = max(original_pair_values.values())
    assert original_a9_residual == 16 >= 12

    a8 = audit_a8_response_wall()
    timing = audit_timing_walls()

    print("A9 wall: Q(A),Q(C),Q(D),q6,X =", qnorm(A9), qnorm(C), qnorm(D), q6, qnorm(D) - q6)
    print("order-six labelled minimizers / optimal completions =", len(minimizers), len(completions))
    print("active oriented-state count per completion:", sorted(active_count_hist.items()))
    print("active orientation histogram (sigma,count):", sorted(orientation_hist.items()))
    print("persistent-state count per completion:", sorted(persistent_count_hist.items()))
    print("aggregate active profile (cross,G,deltaG,D,deltaD,Aenergy,point_gain):")
    for key, count in sorted(active_hist.items()):
        print(" ", count, key)
    print("persistent (cross,D,G) profile:", sorted(persistent_profile_hist.items()))
    print("old response score range on hybrid-active states:", response_min_old_score, response_max_old_score)
    print("max accessible endpoint residual over all / persistent active states:", sorted(accessible_hist.items()), sorted(persistent_accessible_hist.items()))
    print("completions whose every persistent state has only residual 8:", completions_without_persistent_high_pair)
    print("endpoint (+grounds,-grounds,oriented-shore rows) per completion:", sorted(endpoint_count_hist.items()))
    print("endpoint pair |cross-Gram| histogram:", sorted(endpoint_pair_hist.items()))
    print("original A9 endpoint pair |cross-Gram| histogram:", sorted(Counter(original_pair_values.values()).items()))
    print("aggregate endpoint/deletion profile (|b|,m,QX,|h|,sign(centered)):")
    for key, count in sorted(endpoint_hist.items()):
        print(" ", count, key)
    print("A8 response counterexample (Delta, original residual, captured excess) =", a8["response_rise"], a8["original_residual"], a8["captured_excess"])
    print("A5/A7/A8 timing checks:", timing)

    if args.full:
        print("completion table: code P N active persistent persistent-r8 persistent-r16 |b|=16-pairs")
        for row in completion_rows:
            print(f"  0x{row[0]:04x}", *row[1:])
    if args.dump_active:
        print("active dump: code sigma spin cross G deltaG D deltaD Aenergy point_gain")
        for code, row in active_dump:
            spin = "".join("+" if value == 1 else "-" for value in row["z"])
            print(
                f"  0x{code:04x}", row["sigma"], spin, row["c"], row["g"],
                row["delta_g"], row["d"], row["delta_d"], row["old_total"],
                row["point_gain"],
            )


if __name__ == "__main__":
    main()
