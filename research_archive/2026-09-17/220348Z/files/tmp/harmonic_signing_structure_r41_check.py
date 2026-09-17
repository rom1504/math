#!/usr/bin/env python3
"""Wave 41: exact one-deletion completion/Gibbs constraints and OR audit."""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import A4
from harmonic_isotonic_r40_check import audit
from harmonic_resonance_r34_check import Endpoint


def signings(n: int):
    pairs = list(itertools.combinations(range(n), 2))
    for signs in itertools.product((-1, 1), repeat=len(pairs)):
        A = np.zeros((n, n), dtype=np.int64)
        for (i, j), value in zip(pairs, signs):
            A[i, j] = A[j, i] = value
        yield A


def projective_spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def cap(A: np.ndarray) -> int:
    return max(abs(int(x @ A @ x)) for x in projective_spins(len(A)))


def row_fields(A: np.ndarray, x: np.ndarray) -> np.ndarray:
    return x * (A @ x)


def one_deletion_U(A: np.ndarray, x: np.ndarray, beta: float) -> float:
    h = row_fields(A, x)
    return float(np.mean(1.0 / np.cosh(2.0 * beta * h)))


def exact_identity_audit(A: np.ndarray, beta: float) -> dict[str, float]:
    n = len(A)
    Q = cap(A)
    max_component_error = 0.0
    max_odds_error = 0.0
    max_edge_bound_error = 0.0
    max_l1_over_2Q = 0.0
    max_floor_error = 0.0
    max_crossing_odds = 0.0
    for x in projective_spins(n):
        E = int(x @ A @ x)
        h = row_fields(A, x)
        U = one_deletion_U(A, x, beta)
        sum_abs = float(np.sum(np.abs(h)))
        max_l1_over_2Q = max(max_l1_over_2Q, sum_abs / (2.0 * Q))
        floor = math.exp(-2.0 * beta * sum_abs / n)
        max_floor_error = max(max_floor_error, floor - U)
        for i in range(n):
            xf = x.copy()
            xf[i] *= -1
            Ef = int(xf @ A @ xf)
            Uf = one_deletion_U(A, xf, beta)
            omega = beta * (Ef - E)
            chi = math.log(Uf / U)
            component = 1.0 / math.cosh(2.0 * beta * float(h[i]))
            max_odds_error = max(max_odds_error, abs(omega + 4.0 * beta * h[i]))
            max_component_error = max(
                max_component_error,
                abs(component - 1.0 / math.cosh(omega / 2.0)),
            )
            rhs = math.log(n) + math.log(math.cosh(abs(omega) / 2.0))
            max_edge_bound_error = max(max_edge_bound_error, abs(chi) - rhs)
            if omega * (omega + chi) <= 0.0 and abs(chi) > 1e-14:
                max_crossing_odds = max(max_crossing_odds, abs(omega), abs(chi))
    assert max_component_error < 1e-13
    assert max_odds_error < 1e-13
    assert max_edge_bound_error < 1e-13
    assert max_l1_over_2Q <= 1.0 + 1e-13
    assert max_floor_error < 1e-13
    assert max_crossing_odds <= 2.0 * math.log(n) + 1e-12
    return {
        "Q": float(Q),
        "component_odds_error": max_component_error,
        "edge_odds_error": max_odds_error,
        "edge_bound_error": max_edge_bound_error,
        "max_sum_abs_over_2Q": max_l1_over_2Q,
        "floor_error": max_floor_error,
        "max_crossing_odds_or_gap": max_crossing_odds,
        "osc_log_U": max(
            math.log(one_deletion_U(A, x, beta)) for x in projective_spins(n)
        )
        - min(math.log(one_deletion_U(A, x, beta)) for x in projective_spins(n)),
        "osc_bound": 4.0 * beta * Q / n,
    }


def single_component_or(A: np.ndarray, omitted: int, beta: float) -> list[float]:
    """Gauge the omitted star positive and list b_i on the child cut cube."""
    n = len(A)
    gauge = np.ones(n, dtype=np.int64)
    for j in range(n):
        if j != omitted:
            gauge[j] = A[omitted, j]
    B = gauge[:, None] * A * gauge[None, :]
    retained = [j for j in range(n) if j != omitted]
    anchor = retained[0]
    others = retained[1:]
    values = []
    for tail in itertools.product((-1, 1), repeat=len(others)):
        x = np.ones(n, dtype=np.int64)
        x[anchor] = 1
        x[others] = tail
        h = int(sum(B[omitted, j] * x[j] for j in retained))
        values.append(1.0 / math.cosh(2.0 * beta * h))
    return values


def order_four_enumeration(beta: float = 1.0) -> dict[str, object]:
    all_data = [(cap(A), A) for A in signings(4)]
    q4 = min(Q for Q, _ in all_data)
    minimizers = [A for Q, A in all_data if Q == q4]
    max_large_fields = 0
    max_minimizer_ratio = 0.0
    p = 1.0 / math.cosh(2.0 * beta)
    q = 1.0 / math.cosh(6.0 * beta)
    for A in minimizers:
        vals = []
        for x in projective_spins(4):
            h = np.abs(row_fields(A, x))
            assert set(map(int, h)).issubset({1, 3})
            large = int(np.sum(h == 3))
            max_large_fields = max(max_large_fields, large)
            vals.append(one_deletion_U(A, x, beta))
        max_minimizer_ratio = max(max_minimizer_ratio, max(vals) / min(vals))
    assert q4 == 8
    assert max_large_fields == 2
    assert max_minimizer_ratio <= 2.0 + 1e-13

    all_plus = np.ones((4, 4), dtype=np.int64) - np.eye(4, dtype=np.int64)
    plus_vals = [one_deletion_U(all_plus, x, beta) for x in projective_spins(4)]
    assert cap(all_plus) == 12
    assert abs(min(plus_vals) - q) < 1e-13
    assert abs(max(plus_vals) - p) < 1e-13

    component = single_component_or(A4, 0, beta)
    rounded = sorted(round(v, 14) for v in component)
    assert abs(rounded[0] - q) < 1e-13
    assert all(abs(v - p) < 1e-13 for v in rounded[1:])

    # Four-point determinant for G=1/b on the gauged two-bit child chart.
    G = [1.0 / value for value in component]
    determinant = G[0] * G[3] - G[1] * G[2]
    assert abs(abs(determinant) - math.sinh(4.0 * beta) ** 2) < 2e-11
    return {
        "q4": q4,
        "number_minimizers": len(minimizers),
        "max_number_abs_field_3_on_minimizer": max_large_fields,
        "max_minimizer_U_ratio": max_minimizer_ratio,
        "all_plus_Q": cap(all_plus),
        "all_plus_U_ratio": max(plus_vals) / min(plus_vals),
        "single_component_values": component,
        "four_point_determinant": determinant,
    }


def numerical_transport(beta: float) -> dict[str, dict[str, float | int]]:
    all_plus = np.ones((4, 4), dtype=np.int64) - np.eye(4, dtype=np.int64)
    return {
        "exact_A4": audit(Endpoint(A4, beta, 3), order=96),
        "nonminimal_all_plus": audit(Endpoint(all_plus, beta, 3), order=96),
    }


def main() -> None:
    all_plus = np.ones((4, 4), dtype=np.int64) - np.eye(4, dtype=np.int64)
    print("A4 identity", exact_identity_audit(A4, 1.0))
    print("all-plus identity", exact_identity_audit(all_plus, 1.0))
    print("order-four", order_four_enumeration())
    for beta in (1.0, 2.0, 4.0):
        print("transport beta", beta, numerical_transport(beta))
    print("PASS harmonic_signing_structure_r41_check")


if __name__ == "__main__":
    main()
