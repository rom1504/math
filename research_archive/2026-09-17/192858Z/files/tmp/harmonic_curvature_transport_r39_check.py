#!/usr/bin/env python3
"""Finite checks for the Wave 39 curvature-transport identity.

All computations are numerical finite-state audits.  They verify exact
algebraic identities after quadrature; displayed sizes are not asymptotic
evidence.
"""

from __future__ import annotations

import sys

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit, logsumexp

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import A4, A6, A8, A9
from harmonic_resonance_r34_check import Endpoint


def law(H: Endpoint, t: float) -> np.ndarray:
    z = H.lognu + t * H.g
    z -= float(logsumexp(z))
    return np.exp(z)


def load_and_curvature(H: Endpoint, t: float) -> tuple[np.ndarray, np.ndarray]:
    """Statewise vertex KL load L_t and instantaneous curvature V_t."""
    load = np.zeros(H.N)
    curvature = np.zeros(H.N)
    for edges in H.edges[1:]:
        x, y = edges[:, 0], edges[:, 1]
        u0 = H.lognu[y] - H.lognu[x]
        dg = H.g[y] - H.g[x]
        ut = u0 + t * dg
        p = expit(ut)
        cost = p * (ut - u0) - (
            np.logaddexp(0.0, ut) - np.logaddexp(0.0, u0)
        )
        local_curvature = p * (1.0 - p) * dg**2
        load[x] += cost
        load[y] += cost
        curvature[x] += local_curvature
        curvature[y] += local_curvature
    return load, curvature


def regression_groups(
    H: Endpoint, h: np.ndarray, mu_s: np.ndarray, mu_1: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Regression E[h|g] and the two scalar g laws.

    Endpoint construction uses floating log-sum-exp values.  Rounding g to
    11 decimals merges only symmetry copies split by roundoff.  The resulting
    audit tolerance is numerical, while the mathematical identity conditions
    on exact g levels.
    """
    keys = np.round(H.g, 11)
    values, inverse = np.unique(keys, return_inverse=True)
    m = len(values)
    regression = np.zeros(m)
    ps = np.zeros(m)
    p1 = np.zeros(m)
    for j in range(m):
        idx = inverse == j
        # Conditional law given g is independent of interpolation time.  Use
        # the base law to avoid privileging either endpoint numerically.
        base = H.nu[idx]
        regression[j] = float(base @ h[idx] / np.sum(base))
        ps[j] = float(np.sum(mu_s[idx]))
        p1[j] = float(np.sum(mu_1[idx]))
    ps /= np.sum(ps)
    p1 /= np.sum(p1)
    return values, regression, ps, p1


def monotone_quantile_transport(
    values: np.ndarray,
    regression: np.ndarray,
    ps: np.ndarray,
    p1: np.ndarray,
) -> tuple[float, float, float]:
    """Signed and downward regression changes under the quantile coupling."""
    # The endpoint law is an increasing exponential tilt of the time-s law.
    assert np.max(np.cumsum(p1) - np.cumsum(ps)) < 2e-9
    i = j = 0
    ai = float(ps[0])
    bj = float(p1[0])
    signed = downward = order_violation = 0.0
    while i < len(ps) and j < len(p1):
        w = min(ai, bj)
        delta = float(regression[j] - regression[i])
        signed += w * delta
        downward += w * max(-delta, 0.0)
        order_violation = max(order_violation, float(values[i] - values[j]))
        ai -= w
        bj -= w
        if ai < 2e-15:
            i += 1
            if i < len(ps):
                ai = float(ps[i])
        if bj < 2e-15:
            j += 1
            if j < len(p1):
                bj = float(p1[j])
    assert order_violation < 2e-9
    direct = float(p1 @ regression - ps @ regression)
    assert abs(signed - direct) < 3e-9
    return signed, downward, order_violation


def weighted_isotonic_residual(
    values: np.ndarray, weights: np.ndarray
) -> float:
    """Weighted L2 residual from the nondecreasing cone (PAV algorithm)."""
    blocks: list[list[float | int]] = []
    for index, (value, weight) in enumerate(zip(values, weights)):
        blocks.append([index, index + 1, float(weight), float(weight * value)])
        while len(blocks) >= 2:
            left, right = blocks[-2], blocks[-1]
            if float(left[3]) / float(left[2]) <= float(right[3]) / float(right[2]):
                break
            blocks[-2:] = [[
                int(left[0]),
                int(right[1]),
                float(left[2]) + float(right[2]),
                float(left[3]) + float(right[3]),
            ]]
    fitted = np.zeros_like(values)
    for start, stop, weight, total in blocks:
        fitted[int(start):int(stop)] = float(total) / float(weight)
    return float(weights @ (values - fitted) ** 2)


def audit(H: Endpoint, order: int = 96) -> dict[str, float]:
    nodes, weights = leggauss(order)
    mu1 = law(H, 1.0)
    signed_covariance = 0.0
    curvature_transport = 0.0
    weighted_energy = 0.0
    endpoint_curvature_cost = 0.0
    quantile_signed = 0.0
    quantile_downward = 0.0
    mutual_information_bound = 0.0
    isotonic_information_bound = 0.0
    curvature_isotonic_scale = 0.0
    parent_entropy = 0.0
    max_regression_error = 0.0

    for node, weight0 in zip(nodes, weights):
        s = float((node + 1.0) / 2.0)
        weight = float(weight0 / 2.0)
        mus = law(H, s)
        load, curvature = load_and_curvature(H, s)
        mean_g = float(mus @ H.g)
        covariance = float(mus @ ((H.g - mean_g) * load))
        transport = float(mu1 @ curvature - mus @ curvature)
        parent_entropy += weight * s * float(mus @ (H.g - mean_g) ** 2)

        signed_covariance += weight * covariance
        curvature_transport += weight * s * transport
        weighted_energy += weight * s * float(mus @ curvature)
        endpoint_curvature_cost += weight * s * float(mu1 @ curvature)

        values, regression, ps, p1 = regression_groups(H, curvature, mus, mu1)
        direct_regression = float(p1 @ regression - ps @ regression)
        max_regression_error = max(max_regression_error, abs(direct_regression - transport))
        q_signed, q_down, _ = monotone_quantile_transport(
            values, regression, ps, p1
        )
        quantile_signed += weight * s * q_signed
        quantile_downward += weight * s * q_down
        mixture = 0.5 * (ps + p1)
        js = 0.5 * float(
            np.sum(ps * np.log(ps / mixture))
            + np.sum(p1 * np.log(p1 / mixture))
        )
        mean_regression = float(mixture @ regression)
        variance_regression = float(
            mixture @ (regression - mean_regression) ** 2
        )
        information_bound = float(np.sqrt(8.0 * js * variance_regression))
        assert abs(direct_regression) <= information_bound + 2e-10
        mutual_information_bound += weight * s * information_bound
        isotonic_residual = weighted_isotonic_residual(regression, mixture)
        isotonic_bound = float(np.sqrt(8.0 * js * isotonic_residual))
        assert max(-direct_regression, 0.0) <= isotonic_bound + 2e-10
        isotonic_information_bound += weight * s * isotonic_bound
        curvature_isotonic_scale += weight * float(
            np.sqrt(s * (1.0 - s) * isotonic_residual)
        )

    endpoint_load, _ = load_and_curvature(H, 1.0)
    endpoint_cost_direct = float(mu1 @ endpoint_load)
    errors = {
        "covariance_vs_curvature_transport": abs(
            signed_covariance - curvature_transport
        ),
        "endpoint_bregman": abs(endpoint_cost_direct - endpoint_curvature_cost),
        "energy_boundary": abs(
            weighted_energy - (endpoint_cost_direct - signed_covariance)
        ),
        "quantile_signed": abs(quantile_signed - curvature_transport),
        "regression": max_regression_error,
    }
    assert max(errors.values()) < 3e-8, errors
    assert max(-signed_covariance, 0.0) <= quantile_downward + 3e-9
    bootstrap_bound = float(
        np.sqrt(2.0 * parent_entropy) * curvature_isotonic_scale
    )
    assert max(-signed_covariance, 0.0) <= bootstrap_bound + 3e-9
    return {
        "signed_covariance": signed_covariance,
        "curvature_transport": curvature_transport,
        "negative_signed": max(-signed_covariance, 0.0),
        "quantile_downward": quantile_downward,
        "mutual_information_bound": mutual_information_bound,
        "isotonic_information_bound": isotonic_information_bound,
        "curvature_isotonic_scale_squared": curvature_isotonic_scale**2,
        "parent_entropy": parent_entropy,
        "bootstrap_bound": bootstrap_bound,
        "weighted_energy": weighted_energy,
        "endpoint_cost": endpoint_cost_direct,
        "max_error": max(errors.values()),
    }


def main() -> None:
    for name, matrix, beta, m in [
        ("A4", A4, 0.5, 3),
        ("A6", A6, 0.5, 3),
        ("A8", A8, 2.0, 5),
        ("A9", A9, 2.0, 7),
    ]:
        print(name, "beta", beta, "m", m, audit(Endpoint(matrix, beta, m)))
    print("PASS harmonic_curvature_transport_r39_check")


if __name__ == "__main__":
    main()
