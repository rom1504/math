#!/usr/bin/env python3
"""Exact/reproducible checks for Wave 48C cross-block consistency.

The script checks the Bellman projection identity for completion slacks,
the exact A8 two-block dual incompatibility and ground-layer migration,
some instances of the fractional edge-cover certificate, and a deterministic
random-objective sample of exact order-ten minimizers.
"""

from __future__ import annotations

import itertools
import sys
from collections import Counter

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, milp

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import A6, A8, A9  # noqa: E402
from tropical_decomposition_r46_check import (  # noqa: E402
    edges,
    matrix_from_normalized_code,
    qnorm,
    spins,
)


def restriction_key(sigma: int, x: np.ndarray, S: tuple[int, ...]) -> tuple[int, ...]:
    """Gauge-free oriented restriction: sigma and sigma*x_i*x_j."""
    return (int(sigma),) + tuple(
        int(sigma * x[i] * x[j]) for i, j in itertools.combinations(S, 2)
    )


def restrict_existing_key(
    key: tuple[int, ...], S: tuple[int, ...], R: tuple[int, ...]
) -> tuple[int, ...]:
    sigma = key[0]
    pair_value = dict(zip(itertools.combinations(S, 2), key[1:]))
    return (sigma,) + tuple(pair_value[e] for e in itertools.combinations(R, 2))


def bellman_projection_audit(A: np.ndarray) -> int:
    """Exhaust all R subset S and verify s_R=min_{extensions}s_S."""
    n = len(A)
    X = spins(n).astype(np.int64)
    raw = np.einsum("bi,ij,bj->b", X, A.astype(np.int64), X, optimize=True)
    states = [(sigma, x, int(sigma * e)) for sigma in (-1, 1) for x, e in zip(X, raw)]
    q = max(e for _, _, e in states)

    subsets = [tuple(i for i in range(n) if (mask >> i) & 1) for mask in range(1 << n)]
    lifts: dict[tuple[int, ...], dict[tuple[int, ...], int]] = {}
    for S in subsets:
        profile: dict[tuple[int, ...], int] = {}
        for sigma, x, energy in states:
            key = restriction_key(sigma, x, S)
            profile[key] = max(profile.get(key, -10**9), energy)
        lifts[S] = profile

    checked = 0
    for S in subsets:
        for R in subsets:
            if not set(R).issubset(S):
                continue
            projected: dict[tuple[int, ...], int] = {}
            for key_S, lift_energy in lifts[S].items():
                key_R = restrict_existing_key(key_S, S, R)
                projected[key_R] = max(projected.get(key_R, -10**9), lift_energy)
            assert projected == lifts[R]
            # Equivalent slack form, stated explicitly.
            for key_R, lift_energy in lifts[R].items():
                rhs = min(
                    q - value
                    for key_S, value in lifts[S].items()
                    if restrict_existing_key(key_S, S, R) == key_R
                )
                assert q - lift_energy == rhs
            checked += 1
    return checked


def oriented_ground_features(A: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    n = len(A)
    ee = edges(n)
    X = spins(n).astype(np.int64)
    raw = np.einsum("bi,ij,bj->b", X, A.astype(np.int64), X, optimize=True)
    q = int(np.max(np.abs(raw)))
    ground_ids = np.flatnonzero(np.abs(raw) == q)
    sigmas = np.where(raw[ground_ids] > 0, 1, -1).astype(np.int64)
    features = np.asarray(
        [
            [sigma * A[i, j] * X[g, i] * X[g, j] for i, j in ee]
            for g, sigma in zip(ground_ids, sigmas)
        ],
        dtype=np.int64,
    )
    assert np.all(2 * np.sum(features, axis=1) == q)
    return features, ground_ids, q


def a8_cross_block_certificate() -> dict[str, object]:
    """Audit the two genuinely bad maximal A8 four-blocks exactly."""
    A = A8
    n = len(A)
    ee = edges(n)
    B, ground_ids, q = oriented_ground_features(A)
    assert q == 20 and len(ground_ids) == 8

    selectors = list(itertools.combinations(range(n), 4))
    norms = [qnorm(A[np.ix_(S, S)]) for S in selectors]
    qstar = max(norms)
    maximal = [S for S, value in zip(selectors, norms) if value == qstar]
    bad: list[tuple[int, ...]] = []
    zero_supports: list[set[int]] = []
    for S in maximal:
        ids = [ee.index(e) for e in itertools.combinations(S, 2)]
        child = 2 * np.sum(B[:, ids], axis=1)
        if not np.any(child == qstar):
            bad.append(S)
            assert set(map(int, child)) == {0, 4}
            zero = set(int(i) for i in np.flatnonzero(child == 0))
            zero_supports.append(zero)
            # Uniform mass on the four zero-energy lifts is a zero dual.
            assert len(zero) == 4
            assert np.all(np.sum(B[sorted(zero)][:, ids], axis=0) == 0)

    expected_bad = [(0, 3, 4, 5), (1, 2, 6, 7)]
    assert bad == expected_bad
    assert zero_supports[0].isdisjoint(zero_supports[1])
    assert zero_supports[0] | zero_supports[1] == set(range(len(ground_ids)))

    union_edges = [e for S in bad for e in itertools.combinations(S, 2)]
    union_ids = [ee.index(e) for e in union_edges]
    pressure_edges = [(0, 5), (3, 4), (1, 6), (2, 7)]
    pressure_ids = [ee.index(e) for e in pressure_edges]
    # This four-edge overlap certificate is pointwise on the global ground face.
    assert np.all(np.sum(B[:, pressure_ids], axis=1) == 2)
    # Hence every common law has some union-edge mean >=1/2.  Uniform mass on
    # all eight grounds attains <=1/2 on every union edge, proving equality.
    uniform_means = np.sum(B[:, union_ids], axis=0) / len(B)
    assert np.max(uniform_means) == 0.5

    # Independently solve the finite minimax LP and recover the same value.
    c = np.r_[np.zeros(len(B)), 1.0]
    Aub = np.c_[B[:, union_ids].T, -np.ones(len(union_ids))]
    Aeq = np.zeros((1, len(B) + 1))
    Aeq[0, : len(B)] = 1
    lp = linprog(
        c,
        A_ub=Aub,
        b_ub=np.zeros(len(union_ids)),
        A_eq=Aeq,
        b_eq=[1],
        bounds=[(0, None)] * len(B) + [(None, None)],
        method="highs",
    )
    assert lp.success and abs(lp.fun - 0.5) < 1e-10

    # Flip the pressure edges.  All old grounds fall by eight, but the old
    # slack-eight layer moves up by eight and the exact norm remains 20.
    X = spins(n).astype(np.int64)
    raw = np.einsum("bi,ij,bj->b", X, A.astype(np.int64), X, optimize=True)
    old_ground_count = 0
    migrated_count = 0
    new_values = []
    for sigma in (-1, 1):
        energy = sigma * raw
        flip_sum = np.asarray(
            [sum(sigma * A[e] * x[e[0]] * x[e[1]] for e in pressure_edges) for x in X],
            dtype=np.int64,
        )
        new_energy = energy - 4 * flip_sum
        old = energy == q
        assert np.all(flip_sum[old] == 2)
        assert np.all(new_energy[old] == 12)
        migrated = new_energy == q
        assert np.all(energy[migrated] == 12)
        assert np.all(flip_sum[migrated] == -2)
        old_ground_count += int(np.sum(old))
        migrated_count += int(np.sum(migrated))
        new_values.extend(int(v) for v in new_energy)
    assert max(new_values) == q
    assert old_ground_count == migrated_count == 8

    return {
        "q": q,
        "q_star_m4": qstar,
        "maximal_m4_selectors": len(maximal),
        "bad_selectors": bad,
        "individual_zero_dual_supports": [sorted(z) for z in zero_supports],
        "pressure_edges": pressure_edges,
        "common_dual_minimax": "1/2",
        "old_ground_to_new_energy": 12,
        "migrated_slack_eight_states": migrated_count,
    }


def complement_size_fractional_bound(A: np.ndarray, m: int) -> dict[str, int]:
    """Return gamma_frac >= max(q-2|E(K_n)\E(K_m)|,0)."""
    n = len(A)
    q = qnorm(A)
    outside_edges = n * (n - 1) // 2 - m * (m - 1) // 2
    return {
        "n": n,
        "m": m,
        "q": q,
        "outside_edges": outside_edges,
        "gamma_frac_lower_bound": max(q - 2 * outside_edges, 0),
    }


def edge_uniform_cover(A: np.ndarray, m: int) -> dict[str, object]:
    """Find a fractional edge-uniform cover by maximal m-selectors."""
    n = len(A)
    ee = edges(n)
    selectors = list(itertools.combinations(range(n), m))
    norms = [qnorm(A[np.ix_(S, S)]) for S in selectors]
    qstar = max(norms)
    maximal = [S for S, value in zip(selectors, norms) if value == qstar]
    incidence = np.asarray(
        [[int(i in S and j in S) for S in maximal] for i, j in ee], dtype=float
    )
    rho = m * (m - 1) / (n * (n - 1))
    Aeq = np.r_[np.ones((1, len(maximal))), incidence]
    beq = np.r_[1.0, np.full(len(ee), rho)]
    ans = linprog(
        np.zeros(len(maximal)), A_eq=Aeq, b_eq=beq, bounds=(0, None), method="highs"
    )
    q = qnorm(A)
    certifies = bool(ans.success and rho * q > qstar - 4 + 1e-10)
    return {
        "n": n,
        "m": m,
        "q": q,
        "q_star": qstar,
        "maximal_selectors": len(maximal),
        "uniform_cover": bool(ans.success),
        "rho_q": rho * q,
        "q_star_minus_four": qstar - 4,
        "certifies_tight_pair": certifies,
    }


def a8_m4_uniform_cover_exact_wall() -> list[tuple[int, ...]]:
    """Exact incidence contradiction for a nonzero constant edge load."""
    maximal = [
        S
        for S in itertools.combinations(range(8), 4)
        if qnorm(A8[np.ix_(S, S)]) == 12
    ]
    expected = [
        (0, 1, 2, 5),
        (0, 1, 3, 6),
        (0, 2, 4, 7),
        (0, 3, 4, 5),
        (0, 5, 6, 7),
        (1, 2, 6, 7),
        (1, 3, 4, 7),
        (1, 4, 5, 6),
        (2, 3, 4, 6),
        (2, 3, 5, 7),
    ]
    assert maximal == expected

    # If their weights are l_0,...,l_9 and every edge has the same load r,
    # the pair-incidence equations below force
    # l_1=l_2=l_5=l_7=l_9=:b and
    # l_0=l_3=l_4=l_6=l_8=:a.
    class_a = {0, 3, 4, 6, 8}
    class_b = {1, 2, 5, 7, 9}
    incidence = {
        e: {j for j, S in enumerate(maximal) if e[0] in S and e[1] in S}
        for e in edges(8)
    }
    assert incidence[(0, 1)] == {0, 1}  # load a+b
    assert incidence[(0, 5)] == {0, 3, 4}  # load 3a
    assert incidence[(1, 6)] == {1, 5, 7}  # load 3b
    # Equal loads therefore give a+b=3a=3b, hence b=2a and a=2b,
    # forcing a=b=r=0.  Check the asserted two-color reduction directly:
    for e, selector_ids in incidence.items():
        if len(selector_ids) == 2:
            assert len(selector_ids & class_a) == 1
            assert len(selector_ids & class_b) == 1
        elif len(selector_ids) == 3:
            assert selector_ids in ({0, 3, 4}, {1, 5, 7}, {2, 5, 9}, {3, 6, 8})
        else:
            raise AssertionError((e, selector_ids))
    return maximal


def exhaustive_triangle_audit(max_n: int = 6) -> dict[int, int]:
    """Exhaust switching-normalized signings and all grounds for m=3."""
    counts: dict[int, int] = {}
    for n in range(3, max_n + 1):
        variable_edges = (n - 1) * (n - 2) // 2
        total = 1 << variable_edges
        for code in range(total):
            A = matrix_from_normalized_code(n, code)
            X = spins(n).astype(np.int64)
            raw = np.einsum("bi,ij,bj->b", X, A.astype(np.int64), X, optimize=True)
            q = int(np.max(np.abs(raw)))
            for x, energy in zip(X, raw):
                if abs(int(energy)) != q:
                    continue
                sigma = 1 if energy > 0 else -1
                assert any(
                    all(sigma * A[i, j] * x[i] * x[j] == 1 for i, j in itertools.combinations(S, 2))
                    for S in itertools.combinations(range(n), 3)
                )
        counts[n] = total
    return counts


def sample_order_ten(seed: int) -> np.ndarray:
    """Use the exact M_10=13 cap and a deterministic random linear objective."""
    n = 10
    half_cap = 13
    X = spins(n).astype(np.int64)
    variable_edges = [(i, j) for i in range(1, n) for j in range(i + 1, n)]
    products = np.asarray(
        [[x[i] * x[j] for i, j in variable_edges] for x in X], dtype=float
    )
    fixed = np.sum(X[:, 1:], axis=1)
    constant = fixed - np.sum(products, axis=1)
    constraints = LinearConstraint(2 * products, -half_cap - constant, half_cap - constant)
    objective = np.random.default_rng(seed).normal(size=len(variable_edges))
    ans = milp(
        objective,
        integrality=np.ones(len(variable_edges)),
        bounds=Bounds(np.zeros(len(variable_edges)), np.ones(len(variable_edges))),
        constraints=constraints,
        options={"time_limit": 30.0, "mip_rel_gap": 0.0},
    )
    assert ans.success
    bits = np.rint(ans.x).astype(np.int64)
    assert np.max(np.abs(ans.x - bits)) < 1e-7
    A = np.zeros((n, n), dtype=np.int64)
    A[0, 1:] = A[1:, 0] = 1
    for bit, (i, j) in zip(bits, variable_edges):
        A[i, j] = A[j, i] = 2 * bit - 1
    assert qnorm(A) == 26
    return A


def common_parent_report(A: np.ndarray) -> tuple[tuple[int, ...], int, int]:
    """Return maximal-principal profile and all-size common-parent count."""
    n = len(A)
    ee = edges(n)
    X = spins(n).astype(np.int64)
    products = np.asarray([[x[i] * x[j] for i, j in ee] for x in X], dtype=np.int64)
    avec = np.asarray([A[e] for e in ee], dtype=np.int64)
    full = 2 * (products @ avec)
    q = int(np.max(np.abs(full)))
    ground_ids = np.flatnonzero(np.abs(full) == q)
    sigmas = np.where(full[ground_ids] > 0, 1, -1).astype(np.int64)

    masks = np.arange(1 << n, dtype=np.int64)
    sizes = np.asarray([bin(int(mask)).count("1") for mask in masks], dtype=np.int64)
    edge_masks = np.asarray(
        [
            [int(((int(mask) >> i) & 1) and ((int(mask) >> j) & 1)) for i, j in ee]
            for mask in masks
        ],
        dtype=np.int64,
    )
    internal = 2 * (products @ (edge_masks * avec[None, :]).T)
    qmask = np.max(np.abs(internal), axis=0)
    qstar = np.asarray([np.max(qmask[sizes == m]) for m in range(n + 1)], dtype=np.int64)

    all_size = 0
    sizewise_hit = np.zeros(n + 1, dtype=bool)
    for ground, sigma in zip(ground_ids, sigmas):
        good = (qmask == qstar[sizes]) & (sigma * internal[ground] == qstar[sizes])
        by_size = np.asarray([np.any(good & (sizes == m)) for m in range(n + 1)])
        sizewise_hit |= by_size
        all_size += int(np.all(by_size))
    assert np.all(sizewise_hit)
    return tuple(int(z) for z in qstar), len(ground_ids), all_size


def order_ten_sample(seeds: int = 12) -> dict[str, object]:
    profiles: Counter[tuple[int, ...]] = Counter()
    common_counts: Counter[int] = Counter()
    codes: set[tuple[int, ...]] = set()
    for seed in range(seeds):
        A = sample_order_ten(seed)
        profile, grounds, common = common_parent_report(A)
        assert grounds == 40
        assert common > 0
        profiles[profile] += 1
        common_counts[common] += 1
        codes.add(tuple(int(A[i, j]) for i in range(1, 10) for j in range(i + 1, 10)))
    return {
        "samples": seeds,
        "distinct_normalized_signings": len(codes),
        "qstar_profiles": {str(key): value for key, value in sorted(profiles.items())},
        "all_size_common_parent_counts": dict(sorted(common_counts.items())),
    }


def main() -> None:
    bellman_pairs = bellman_projection_audit(A8)
    print("Bellman nested-pair checks on A8", bellman_pairs)

    a8 = a8_cross_block_certificate()
    print("A8 cross-block certificate", a8)

    complement_bounds = [
        complement_size_fractional_bound(A8, 4),
        complement_size_fractional_bound(A8, 7),
        complement_size_fractional_bound(A9, 8),
    ]
    assert complement_bounds[0]["gamma_frac_lower_bound"] == 0
    assert complement_bounds[1]["gamma_frac_lower_bound"] == 6
    assert complement_bounds[2]["gamma_frac_lower_bound"] == 8
    print("complement-size fractional bounds", complement_bounds)

    covers = [
        edge_uniform_cover(A6, 5),
        edge_uniform_cover(A8, 7),
        edge_uniform_cover(A9, 8),
        edge_uniform_cover(A8, 4),
    ]
    assert covers[0]["certifies_tight_pair"]
    assert covers[1]["certifies_tight_pair"]
    assert not covers[2]["certifies_tight_pair"]
    assert not covers[3]["uniform_cover"]
    exact_wall = a8_m4_uniform_cover_exact_wall()
    print("fractional edge-cover cases", covers)
    print("A8,m4 exact no-constant-load selector order", exact_wall)

    triangles = exhaustive_triangle_audit()
    print("all-signing triangle audit", triangles)

    ten = order_ten_sample()
    print("order-ten deterministic MILP sample", ten)
    print("PASS cross_block_dual_r48_check")


if __name__ == "__main__":
    main()
