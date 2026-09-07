#!/usr/bin/env python3
"""Wave 23: nonuniform weighting and LP duals of the exact star cover."""

from __future__ import annotations

import itertools
import math

import numpy as np
from scipy.optimize import linprog

from check_response_dual_r16 import A8, A9, spins


def states(A: np.ndarray):
    n = len(A)
    X0 = spins(n)
    X = np.tile(X0, (2, 1))
    sigma = np.repeat(np.array((-1, 1), dtype=np.int64), len(X0))
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    signs = np.column_stack(
        [sigma * A[i, j] * X[:, i] * X[:, j] for i, j in edges]
    )
    energy = 2 * signs.sum(axis=1)
    q = int(energy.max())
    deficit = q - energy
    return edges, signs, deficit, q


def patterns(m: int) -> np.ndarray:
    return np.array(list(itertools.product((-1, 1), repeat=m)), dtype=np.int64)


def transform_minimax(m: int, beta: float) -> tuple[float, float, float]:
    """min_p max_y F_p(y)/exp(-2 beta sum y)."""
    Y = patterns(m)
    # Subsets S are encoded by a_S=1 off S and -1 on S.  As S ranges over
    # the cube, a_S ranges over Y.
    ratio = np.exp(2 * beta * (Y @ Y.T))
    # Variables are p_a and C.  ratio @ p <= C.
    c = np.r_[np.zeros(len(Y)), 1.0]
    A_ub = np.c_[ratio, -np.ones(len(Y))]
    b_ub = np.zeros(len(Y))
    A_eq = np.zeros((1, len(Y) + 1))
    A_eq[0, : len(Y)] = 1
    result = linprog(
        c,
        A_ub=A_ub,
        b_ub=b_ub,
        A_eq=A_eq,
        b_eq=np.ones(1),
        bounds=[(0, None)] * len(Y) + [(None, None)],
        method="highs",
    )
    assert result.success
    expected = math.cosh(2 * beta) ** m
    uniform_value = float(np.max(ratio @ np.full(len(Y), 1 / len(Y))))
    assert math.isclose(result.fun, expected, rel_tol=2e-8, abs_tol=2e-8)
    assert math.isclose(uniform_value, expected, rel_tol=2e-10)
    return float(result.fun), expected, float(np.max(np.abs(result.x[: len(Y)] - 1 / len(Y))))


def abstract_upper_lp(m: int, beta: float, D: float) -> tuple[float, float]:
    """Strongest BC upper bound from all moment lower bounds plus pair law.

    Variables are the normalized star-pattern masses nu(y).  We impose all
    D_S/D >= 1/D inequalities and exact coordinate-flip detailed balance.
    """
    Y = patterns(m)
    count = len(Y)
    h = Y.sum(axis=1)
    # f_S(y)=exp(-4 beta sum_{e in S} y_e).  Boolean subset masks are
    # represented directly below.
    subsets = np.array(list(itertools.product((0, 1), repeat=m)), dtype=np.int64)
    F = np.exp(-4 * beta * (subsets @ Y.T))
    bc_vector = np.exp(-2 * beta * h)

    # Maximize BC, represented as minimizing its negative.
    c = -bc_vector
    A_ub = -F
    b_ub = -np.full(len(subsets), 1 / D)
    eq_rows = [np.ones(count)]
    eq_rhs = [1.0]
    lookup = {tuple(map(int, y)): idx for idx, y in enumerate(Y)}
    for idx, y in enumerate(Y):
        hi = int(h[idx])
        j = lookup[tuple(map(int, -y))]
        if hi < 0 or (hi == 0 and idx > j):
            continue
        row = np.zeros(count)
        # nu(-y)=nu(y) exp(-4 beta h(y)).
        row[j] = 1
        row[idx] = -math.exp(-4 * beta * hi)
        eq_rows.append(row)
        eq_rhs.append(0.0)
    result = linprog(
        c,
        A_ub=A_ub,
        b_ub=b_ub,
        A_eq=np.asarray(eq_rows),
        b_eq=np.asarray(eq_rhs),
        bounds=(0, None),
        method="highs",
    )
    assert result.success
    optimum = -float(result.fun)
    parity = 1.0 if m % 2 == 0 else 1 / math.cosh(2 * beta)
    assert math.isclose(optimum, parity, rel_tol=2e-8, abs_tol=2e-8)

    # Independently construct the two-center cover/moment wall.  Choose a
    # minimum-absolute-field pattern y, and its antipode with the exact Gibbs
    # pair weights.  Every perturbation moment is at least one.
    epsilon = m % 2
    y = np.r_[np.ones((m + epsilon) // 2, dtype=int), -np.ones((m - epsilon) // 2, dtype=int)]
    a_values = subsets @ y
    p_plus = math.exp(2 * beta * epsilon) / (2 * math.cosh(2 * beta * epsilon))
    p_minus = math.exp(-2 * beta * epsilon) / (2 * math.cosh(2 * beta * epsilon))
    moments = p_plus * np.exp(-4 * beta * a_values) + p_minus * np.exp(4 * beta * a_values)
    assert moments.min() >= 1 - 1e-10
    wall_bc = p_plus * math.exp(-2 * beta * epsilon) + p_minus * math.exp(2 * beta * epsilon)
    assert math.isclose(wall_bc, parity, rel_tol=1e-12)

    # The same pair is an exact abstract restricted-cover wall.  Give y
    # deficit zero and -y deficit 4 epsilon.  Both restricted radii equal
    # floor(m/2), and their antipodal balls cover the entire m-cube.
    rho = (m - epsilon) // 2
    target_signs = 1 - 2 * subsets
    distance_plus = ((target_signs != y).sum(axis=1))
    distance_minus = m - distance_plus
    assert np.all(np.minimum(distance_plus, distance_minus) <= rho)
    assert 0 + 2 * epsilon == 2 * m - 4 * rho
    assert 4 * epsilon - 2 * epsilon == 2 * m - 4 * rho
    return optimum, float(moments.min())


def cover_dual(A: np.ndarray, name: str, beta: float) -> None:
    n = len(A)
    m = n - 1
    edges, signs, deficit, q = states(A)
    edge_index = {edge: k for k, edge in enumerate(edges)}
    weights = np.exp(-beta * deficit)
    D = float(weights.sum())
    target_bits = np.array(list(itertools.product((0, 1), repeat=m)), dtype=np.int64)

    tau_values = []
    G_values = []
    bc_values = []
    max_dual_slack = []
    min_cert_mass = []
    for i in range(n):
        H = np.array(
            [edge_index[tuple(sorted((i, j)))] for j in range(n) if j != i],
            dtype=int,
        )
        local = signs[:, H]
        selected_sums = target_bits @ local.T
        shifted = deficit[None, :] + 4 * selected_sums
        cover = shifted <= 0
        assert np.all(cover.any(axis=1))

        h = local.sum(axis=1)
        geometric = np.exp(-beta * (deficit + 2 * h))
        G = float(geometric.sum())
        bc = G / D
        assert bc <= 1 + 1e-10

        # Fractional cover dual: max sum lambda_S subject to
        # sum_{S certified by omega} lambda_S <= geometric_omega.
        result = linprog(
            -np.ones(len(target_bits)),
            A_ub=cover.T.astype(float),
            b_ub=geometric,
            bounds=(0, None),
            method="highs",
        )
        assert result.success
        tau = -float(result.fun)
        assert tau <= G + 1e-8
        tau_values.append(tau)
        G_values.append(G)
        bc_values.append(bc)
        max_dual_slack.append(float(np.max(cover.T @ result.x - geometric)))
        cert_mass = cover @ (weights / D)
        min_cert_mass.append(float(cert_mass.min()))

    transform_value, expected_transform, p_deviation = transform_minimax(m, beta)
    upper_lp, wall_min_moment = abstract_upper_lp(m, beta, D)
    parity = 1.0 if m % 2 == 0 else 1 / math.cosh(2 * beta)
    print(
        name,
        f"beta={beta}",
        {
            "Q": q,
            "D": D,
            "BC_range": (min(bc_values), max(bc_values)),
            "kappa_range": (
                min(-math.log(v) / beta for v in bc_values),
                max(-math.log(v) / beta for v in bc_values),
            ),
            "weighted_cover_tau_range": (min(tau_values), max(tau_values)),
            "tau_over_G_range": (
                min(t / g for t, g in zip(tau_values, G_values)),
                max(t / g for t, g in zip(tau_values, G_values)),
            ),
            "min_Gibbs_mass_of_a_witness_set": min(min_cert_mass),
            "transform_minimax": transform_value,
            "transform_uniform_expected": expected_transform,
            "minimax_p_max_deviation_from_uniform": p_deviation,
            "all_moment_upper_LP": upper_lp,
            "parity_ceiling": parity,
            "two_center_wall_min_moment": wall_min_moment,
            "dual_feasibility_max_slack": max(max_dual_slack),
        },
    )


def main() -> None:
    for beta in (0.5, 1.0):
        for A, name in ((A8, "A8"), (A9, "A9")):
            cover_dual(A, name, beta)
    print("PASS: nonuniform transform, weighted cover duals, and sharp abstract upper LP")


if __name__ == "__main__":
    main()
