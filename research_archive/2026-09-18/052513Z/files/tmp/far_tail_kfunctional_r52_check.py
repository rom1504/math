#!/usr/bin/env python3
"""Finite algebra audit for the Wave 52 K-functional completion memo.

This exhausts all selectors/local projective states at the listed small
orders.  It checks deterministic identities and probability set inclusions;
it is not asymptotic evidence for the open saved-abundance statement.
"""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from cross_block_dual_r48_check import sample_order_ten
from envelope_block_cover_r27 import A6, A8, A9, projective_spins


def qnorm(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in projective_spins(len(a)))


def exact_kfunctional(beta: np.ndarray, h: int) -> float:
    """K_{1,2}(beta,sqrt(h)) via its elementary dual formula.

    The dual is max beta.z subject to ||z||_infty<=1 and ||z||_2<=sqrt(h).
    After signs and rearrangement, z_i=min(1,b_i/lambda), with lambda chosen
    so sum z_i^2=h, unless the Euclidean constraint is inactive.
    """

    b = np.abs(np.asarray(beta, dtype=np.float64))
    b = b[b > 0]
    if len(b) == 0:
        return 0.0
    if h >= len(b):
        return float(np.sum(b))

    lo = 0.0
    hi = float(np.max(b))
    # At hi the sum can still exceed h, so enlarge until feasible.
    while float(np.sum(np.minimum(1.0, b / hi) ** 2)) > h:
        hi *= 2.0
    for _ in range(100):
        mid = (lo + hi) / 2.0
        value = float(np.sum(np.minimum(1.0, b / mid) ** 2))
        if value > h:
            lo = mid
        else:
            hi = mid
    z = np.minimum(1.0, b / hi)
    assert abs(float(z @ z) - h) < 1e-10
    return float(b @ z)


def profile_proxy(beta: np.ndarray, h: int) -> tuple[float, float, float]:
    b = np.sort(np.abs(np.asarray(beta, dtype=np.float64)))[::-1]
    j = min(h, len(b))
    head = float(np.sum(b[:j]))
    diffuse = math.sqrt(h) * float(np.linalg.norm(b[j:]))
    return head, diffuse, head + diffuse


def audit_matrix(a: np.ndarray, name: str, m: int) -> dict[str, float]:
    n = len(a)
    q = qnorm(a)
    hfun = max(1, min(3, n - m - 1))
    spins_out = list(itertools.product((-1, 1), repeat=n - m))

    states = 0
    completions = 0
    cap_threshold_checks = 0
    union_checks = 0
    max_decomposition_error = 0.0
    max_orthogonality_error = 0.0
    max_l1_minus_q = -math.inf
    min_k_over_proxy = math.inf
    max_k_over_proxy = 0.0

    for s0 in itertools.combinations(range(n), m):
        s = tuple(s0)
        tset = tuple(i for i in range(n) if i not in s)
        child = a[np.ix_(s, s)]
        cross_block = a[np.ix_(tset, s)]
        outside = a[np.ix_(tset, tset)]

        for sigma in (-1, 1):
            for y in projective_spins(m):
                states += 1
                e = sigma * int(y @ child @ y)
                beta = 2 * sigma * (cross_block @ y)
                bmat = sigma * outside

                lvals = []
                qvals = []
                zvals = []
                fullvals = []
                by_spin: dict[tuple[int, ...], tuple[int, int, int]] = {}
                for wt in spins_out:
                    w = np.asarray(wt, dtype=np.int64)
                    lv = int(beta @ w)
                    qv = int(w @ bmat @ w)
                    zv = lv + qv
                    fv = e + zv
                    lvals.append(lv)
                    qvals.append(qv)
                    zvals.append(zv)
                    fullvals.append(fv)
                    by_spin[tuple(int(x) for x in w)] = (lv, qv, zv)

                larr = np.asarray(lvals, dtype=np.float64)
                qarr = np.asarray(qvals, dtype=np.float64)
                zarr = np.asarray(zvals, dtype=np.float64)
                farr = np.asarray(fullvals, dtype=np.float64)
                completions += len(zarr)

                assert np.max(np.abs(farr)) <= q
                max_decomposition_error = max(
                    max_decomposition_error,
                    float(np.max(np.abs(zarr - larr - qarr))),
                )
                assert abs(float(np.mean(qarr))) < 1e-12
                assert abs(float(np.mean(larr))) < 1e-12
                max_orthogonality_error = max(
                    max_orthogonality_error, abs(float(np.mean(larr * qarr)))
                )
                assert abs(float(np.mean(larr**2)) - float(beta @ beta)) < 1e-12
                assert (
                    abs(float(np.mean(qarr**2)) - 2 * (n - m) * (n - m - 1))
                    < 1e-12
                )
                assert abs(float(beta @ beta) - 4 * float((cross_block @ y) @ (cross_block @ y))) < 1e-12

                # Exact antipodal identities and their pointwise cap consequence.
                for wt, (lv, qv, zv) in by_spin.items():
                    neg = tuple(-x for x in wt)
                    lvn, qvn, zvn = by_spin[neg]
                    assert lvn == -lv and qvn == qv
                    assert min(zv, zvn) == qv - abs(lv)
                    assert max(zv, zvn) == qv + abs(lv)
                    assert max(zv, zvn) <= q - e
                    assert min(zv, zvn) <= q - e - 2 * abs(lv)

                # ||beta||_1=max_w L(w)<=q, the circularity audit.
                l1 = float(np.sum(np.abs(beta)))
                assert abs(float(np.max(larr)) - l1) < 1e-12
                max_l1_minus_q = max(max_l1_minus_q, l1 - q)
                assert l1 <= q + 1e-12

                # Exact K-functional sanity checks and Holmstedt-proxy diagnostics.
                kval = exact_kfunctional(beta, hfun)
                head, diffuse, proxy = profile_proxy(beta, hfun)
                assert kval <= l1 + 1e-10
                assert kval <= math.sqrt(hfun) * float(np.linalg.norm(beta)) + 1e-10
                if proxy > 0:
                    ratio = kval / proxy
                    min_k_over_proxy = min(min_k_over_proxy, ratio)
                    max_k_over_proxy = max(max_k_over_proxy, ratio)
                else:
                    assert kval == head == diffuse == 0

                # Probability form of the cap-pairing inequality, for all
                # integer thresholds (event changes occur only on this lattice).
                for r in range(0, 2 * q + 1):
                    lhs = float(np.mean(zarr <= -r))
                    rhs = 0.5 * float(np.mean(np.abs(larr) >= (q - e + r) / 2))
                    assert lhs + 1e-12 >= rhs
                    cap_threshold_checks += 1

                # Dependence-safe union inequality.  For a>=b and r=a-b,
                # {L<=-a,Q<=b} is contained in {L+Q<=-r}; no independence.
                acuts = sorted(set(float(abs(x)) for x in larr))
                bcuts = sorted(set(float(x) for x in qarr))
                for acut in acuts:
                    for bcut in bcuts:
                        if acut + 1e-12 < bcut:
                            continue
                        r = acut - bcut
                        lhs = float(np.mean(zarr <= -r + 1e-12))
                        lower = float(np.mean(larr <= -acut + 1e-12)) - float(
                            np.mean(qarr > bcut + 1e-12)
                        )
                        assert lhs + 1e-12 >= lower
                        union_checks += 1

    out = {
        "states": states,
        "completions": completions,
        "cap_checks": cap_threshold_checks,
        "union_checks": union_checks,
        "decomp_err": max_decomposition_error,
        "orth_err": max_orthogonality_error,
        "l1_minus_q": max_l1_minus_q,
        "min_k_proxy": min_k_over_proxy,
        "max_k_proxy": max_k_over_proxy,
    }
    print(
        f"{name:10s} n={n:2d} m={m:2d} q={q:2d} H={hfun} "
        f"states={states:5d} completions={completions:7d} "
        f"cap-thresholds={cap_threshold_checks:7d} union={union_checks:7d} "
        f"max(||beta||1-q)={max_l1_minus_q:6.1f} "
        f"K/Phi=[{min_k_over_proxy:.6f},{max_k_over_proxy:.6f}]"
    )
    assert max_decomposition_error < 1e-12
    assert max_orthogonality_error < 1e-12
    return out


def exponent_audit() -> None:
    # Exact affine exponent identities used in (R52.7), (R52.10), (R52.11).
    for c in (0.01, 0.10, 0.20, 0.249):
        h_exp = 3 / 4 - c
        frob_exp = 1 + h_exp / 2
        op_exp = 3 / 4 + h_exp
        target_exp = 3 / 2 - c
        variance_threshold_exp = 2 * target_exp - h_exp
        assert abs(op_exp - target_exp) < 1e-12
        assert frob_exp < target_exp
        assert abs(variance_threshold_exp - (9 / 4 - c)) < 1e-12
    print("exponent identities checked for 0<c<1/4")


def main() -> None:
    results = []
    for a, name, m in (
        (A6, "A6", 3),
        (A8, "A8", 4),
        (A9, "A9", 5),
        (sample_order_ten(0), "A10_seed0", 5),
    ):
        results.append(audit_matrix(a, name, m))
    exponent_audit()
    assert sum(int(x["cap_checks"]) for x in results) > 0
    assert sum(int(x["union_checks"]) for x in results) > 0
    print("PASS far-tail K-functional algebra audit")


if __name__ == "__main__":
    main()
