#!/usr/bin/env python3
"""Wave 40 finite audit of scalar curvature regression and antitone residual."""

from __future__ import annotations

import math
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit, logsumexp

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import A6, A8, A9
from harmonic_resonance_r34_check import Endpoint
from harmonic_curvature_transport_r39_check import weighted_isotonic_residual


def law(H: Endpoint, t: float) -> np.ndarray:
    z = H.lognu + t * H.g
    z -= float(logsumexp(z))
    return np.exp(z)


def curvature_and_derivative(H: Endpoint, s: float) -> tuple[np.ndarray, np.ndarray]:
    curvature = np.zeros(H.N)
    derivative = np.zeros(H.N)
    for edges in H.edges[1:]:
        x, y = edges[:, 0], edges[:, 1]
        dg = H.g[y] - H.g[x]
        odds = H.lognu[y] - H.lognu[x] + s * dg
        p = expit(odds)
        v = p * (1.0 - p) * dg**2
        dv = p * (1.0 - p) * (1.0 - 2.0 * p) * dg**3
        curvature[x] += v
        curvature[y] += v
        derivative[x] += dv
        derivative[y] += dv
    return curvature, derivative


def groups(H: Endpoint) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    keys = np.round(H.g, 11)
    values = np.unique(keys)
    indices = [np.flatnonzero(keys == value) for value in values]
    return keys, values, indices


def regression(
    H: Endpoint, statistic: np.ndarray, indices: list[np.ndarray]
) -> np.ndarray:
    ans = np.zeros(len(indices))
    for j, idx in enumerate(indices):
        weights = H.nu[idx]
        ans[j] = float(weights @ statistic[idx] / np.sum(weights))
    return ans


def scalar_laws(
    H: Endpoint, indices: list[np.ndarray], s: float
) -> tuple[np.ndarray, np.ndarray]:
    mus = law(H, s)
    mu1 = law(H, 1.0)
    ps = np.asarray([np.sum(mus[idx]) for idx in indices], dtype=float)
    p1 = np.asarray([np.sum(mu1[idx]) for idx in indices], dtype=float)
    ps /= np.sum(ps)
    p1 /= np.sum(p1)
    return ps, p1


def endpoint_cost(H: Endpoint) -> float:
    mu1 = law(H, 1.0)
    load = np.zeros(H.N)
    for edges in H.edges[1:]:
        x, y = edges[:, 0], edges[:, 1]
        u0 = H.lognu[y] - H.lognu[x]
        dg = H.g[y] - H.g[x]
        u1 = u0 + dg
        p = expit(u1)
        k = p * dg - (np.logaddexp(0.0, u1) - np.logaddexp(0.0, u0))
        load[x] += k
        load[y] += k
    return float(mu1 @ load)


def identity_check(H: Endpoint, s: float = 0.43) -> dict[str, float]:
    curvature, derivative = curvature_and_derivative(H, s)
    h = 2e-6
    plus, _ = curvature_and_derivative(H, s + h)
    minus, _ = curvature_and_derivative(H, s - h)
    derivative_error = float(np.max(np.abs((plus - minus) / (2 * h) - derivative)))

    # Matched edge identity: the score jump is the exclusion log-ratio.
    exclusion_error = 0.0
    replica_error = 0.0
    for bit, edges in enumerate(H.edges[1:], start=1):
        for x, y in edges:
            dg = float(H.g[y] - H.g[x])
            exclusion_error = max(
                exclusion_error,
                abs(dg - math.log(float(H.excl[bit, x] / H.excl[bit, y]))),
            )
            odds = float(H.lognu[y] - H.lognu[x] + s * dg)
            p = float(expit(odds))
            v = p * (1.0 - p) * dg * dg
            two_replica = 0.5 * (
                2.0 * p * (1.0 - p) * dg * dg
            )
            replica_error = max(replica_error, abs(v - two_replica))
    assert derivative_error < 2e-6
    assert exclusion_error < 3e-8
    assert replica_error < 1e-14
    return {
        "s_derivative_error": derivative_error,
        "matched_exclusion_error": exclusion_error,
        "two_replica_error": replica_error,
    }


def audit(H: Endpoint, order: int = 48) -> dict[str, float | int]:
    _, values, indices = groups(H)
    nodes, weights = leggauss(order)
    K = 0.0
    max_drop = 0.0
    nonmonotone_nodes = 0
    max_residual = 0.0
    parent_entropy = 0.0
    for node, weight0 in zip(nodes, weights):
        s = float((node + 1.0) / 2.0)
        weight0 = float(weight0 / 2.0)
        curvature, _ = curvature_and_derivative(H, s)
        m = regression(H, curvature, indices)
        ps, p1 = scalar_laws(H, indices, s)
        mixture = 0.5 * (ps + p1)
        residual = weighted_isotonic_residual(m, mixture)
        K += weight0 * math.sqrt(s * (1.0 - s) * residual)
        max_residual = max(max_residual, residual)
        drops = m[:-1] - m[1:]
        local_drop = float(np.max(drops)) if len(drops) else 0.0
        max_drop = max(max_drop, local_drop)
        if local_drop > 1e-10:
            nonmonotone_nodes += 1
        mus = law(H, s)
        mean = float(mus @ H.g)
        parent_entropy += weight0 * s * float(mus @ (H.g - mean) ** 2)

    c1 = endpoint_cost(H)
    return {
        "levels": len(values),
        "nonmonotone_nodes": nonmonotone_nodes,
        "quadrature_nodes": order,
        "max_adjacent_drop": max_drop,
        "max_isotonic_residual": max_residual,
        "K_squared": K * K,
        "endpoint_cost": c1,
        "K2_over_endpoint": K * K / c1 if c1 > 1e-30 else 0.0,
        "parent_entropy": parent_entropy,
    }


def abstract_endpoint_wall(L: float, order: int = 160) -> dict[str, float]:
    """Exact square-family wall to K^2 <= C C_V(1).

    The base law is uniform on {00,01,10,11}, while g is zero at 00 and L
    at the other three states.  Only the two coordinate edges incident to
    00 have a nonzero score gap.
    """
    # Matched fixed-size erasure lift on three bits.  The score is the same
    # two-bit OR score and is independent of the third bit.  For k=n-m=1,
    # the posterior omission probabilities and unnormalized omitted-component
    # likelihoods are
    #   r_1=r_2=eps exp(-g), r_3=1-2 eps exp(-g), b_i=exp(g) r_i.
    # Thus sum_i r_i=1, b_i is independent of bit i, and every edge obeys
    # Delta_i g = log(r_i(x)/r_i(y)).
    eps = 0.2
    raw_g = np.asarray(
        [L if ((state & 1) or (state & 2)) else 0.0 for state in range(8)]
    )
    r = np.vstack(
        [
            eps * np.exp(-raw_g),
            eps * np.exp(-raw_g),
            1.0 - 2.0 * eps * np.exp(-raw_g),
        ]
    )
    components = np.exp(raw_g)[None, :] * r
    assert np.max(np.abs(np.sum(r, axis=0) - 1.0)) < 1e-14
    assert np.max(
        np.abs(np.sum(components, axis=0) - np.exp(raw_g)) / np.exp(raw_g)
    ) < 1e-12
    exclusion_error = 0.0
    invariance_error = 0.0
    for bit in range(3):
        for state in range(8):
            neighbor = state ^ (1 << bit)
            exclusion_error = max(
                exclusion_error,
                abs(
                    (raw_g[neighbor] - raw_g[state])
                    - math.log(r[bit, state] / r[bit, neighbor])
                ),
            )
            invariance_error = max(
                invariance_error,
                abs(components[bit, state] - components[bit, neighbor]),
            )
    assert exclusion_error < 1e-12
    assert invariance_error < 1e-12

    nodes, weights = leggauss(order)
    K = 0.0
    signed_migration = 0.0
    energy_integral = 0.0
    max_restoring = 0.0
    for node, weight0 in zip(nodes, weights):
        s = float((node + 1.0) / 2.0)
        weight0 = float(weight0 / 2.0)
        exp_minus_sL = math.exp(-s * L)
        q = 1.0 / (1.0 + exp_minus_sL)
        v = L * L * exp_minus_sL / (1.0 + exp_minus_sL) ** 2
        lambda_s_low = exp_minus_sL / (exp_minus_sL + 3.0)
        lambda_1_low = math.exp(-L) / (math.exp(-L) + 3.0)
        mixture_low = 0.5 * (lambda_s_low + lambda_1_low)
        residual = (
            mixture_low
            * (1.0 - mixture_low)
            * (4.0 * v / 3.0) ** 2
        )
        K += weight0 * math.sqrt(s * (1.0 - s) * residual)
        rs = math.exp(s * L)
        r1 = math.exp(L)
        contemporaneous_v = 2.0 * v * (1.0 + rs) / (1.0 + 3.0 * rs)
        endpoint_v = 2.0 * v * (1.0 + r1) / (1.0 + 3.0 * r1)
        signed_migration += weight0 * s * (endpoint_v - contemporaneous_v)
        energy_integral += weight0 * s * contemporaneous_v
        variance = 3.0 * L * L * rs / (1.0 + 3.0 * rs) ** 2
        max_restoring = max(max_restoring, variance / contemporaneous_v)

    q1 = float(expit(L))
    binary_kl = q1 * L - (float(np.logaddexp(0.0, L)) - math.log(2.0))
    endpoint_context_mass = (1.0 + math.exp(L)) / (1.0 + 3.0 * math.exp(L))
    endpoint_cost = 2.0 * endpoint_context_mass * binary_kl

    endpoint_low = 1.0 / (1.0 + 3.0 * math.exp(L))
    endpoint_high = (1.0 - endpoint_low) / 3.0
    entropy = endpoint_low * math.log(4.0 * endpoint_low) + 3.0 * endpoint_high * math.log(
        4.0 * endpoint_high
    )
    return {
        "L": L,
        "K_squared": K * K,
        "K2_over_L": K * K / L,
        "endpoint_cost": endpoint_cost,
        "K2_over_endpoint": K * K / endpoint_cost,
        "parent_entropy": entropy,
        "signed_migration": signed_migration,
        "energy_integral": energy_integral,
        "max_restoring_ratio": max_restoring,
        "matched_exclusion_error": exclusion_error,
        "component_invariance_error": invariance_error,
    }


def main() -> None:
    beta0 = math.log(2.0) / 2.0
    cases = [
        ("A6", A6, beta0, 3),
        ("A6", A6, beta0, 4),
        ("A6", A6, beta0, 5),
        ("A6", A6, 2.0, 3),
        ("A8", A8, 0.5, 3),
        ("A8", A8, 2.0, 4),
        ("A8", A8, 2.0, 6),
        ("A8", A8, 2.0, 7),
        ("A9", A9, 0.5, 3),
        ("A9", A9, 2.0, 3),
        ("A9", A9, 2.0, 7),
        ("A9", A9, 2.0, 8),
    ]
    checked = set()
    for name, matrix, beta, m in cases:
        H = Endpoint(matrix, beta, m)
        key = (name, beta, m)
        identities = identity_check(H) if key not in checked else {}
        checked.add(key)
        print(name, "beta", beta, "m", m, audit(H), identities)
    for L in (4.0, 8.0, 16.0, 32.0, 64.0):
        wall = abstract_endpoint_wall(L)
        print("abstract_square_wall", wall)
    wall64 = abstract_endpoint_wall(64.0)
    assert wall64["K2_over_L"] > 1e-3
    assert abs(wall64["endpoint_cost"] - 2.0 * math.log(2.0) / 3.0) < 1e-12
    assert abs(wall64["parent_entropy"] - math.log(4.0 / 3.0)) < 1e-12
    print("PASS harmonic_isotonic_r40_check")


if __name__ == "__main__":
    main()
