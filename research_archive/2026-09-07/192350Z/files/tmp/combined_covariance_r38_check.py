#!/usr/bin/env python3
"""Wave 38 checks for aggregate matched harmonic covariance identities.

All output and code stay under /home/math/quadra/tmp.  The checks are
finite-state numerical audits of exact algebraic identities; decimal sizes
are not asymptotic evidence.
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


def vertex_state(H: Endpoint, t: float) -> dict[str, object]:
    """Return all vertex-coordinate ingredients at interpolation time t."""
    mu = law(H, t)
    mean_g = float(mu @ H.g)
    load = np.zeros(H.N)
    coordinate_covariances = []
    coordinate_energies = []
    mixture_weights = []
    mixture_scores = []
    mixture_costs = []
    split_b = []
    split_shape = []

    # In the oriented-cut chart bit 0 is global orientation; bits 1,...,n-1
    # are vertex flips.  Keep orientation out of this audit.
    for bit, edges in enumerate(H.edges[1:], start=1):
        x, y = edges[:, 0], edges[:, 1]
        mass = mu[x] + mu[y]
        u0 = H.lognu[y] - H.lognu[x]
        ut = u0 + t * (H.g[y] - H.g[x])
        p = expit(ut)
        score = (1.0 - p) * H.g[x] + p * H.g[y]
        cost = p * (ut - u0) - (np.logaddexp(0.0, ut) - np.logaddexp(0.0, u0))
        variance = p * (1.0 - p) * (H.g[y] - H.g[x]) ** 2

        coordinate_covariances.append(float(mass @ ((score - mean_g) * cost)))
        coordinate_energies.append(float(mass @ variance))
        mixture_weights.append(mass / (H.n - 1))
        mixture_scores.append(score)
        mixture_costs.append(cost)

        # k_i(D_-i) is constant at the two endpoints of each i-edge.
        load[x] += cost
        load[y] += cost

        # Exact matched split A' = log(B_i/E U) + H'_i.
        normalized_b_x = H.f[x] * H.excl[bit, x]
        normalized_b_y = H.f[y] * H.excl[bit, y]
        assert np.max(np.abs(normalized_b_x - normalized_b_y)) < 3e-9
        log_b = np.log(normalized_b_x)
        shape_score = score - log_b
        split_b.append(float(mass @ ((log_b - mass @ log_b) * cost)))
        split_shape.append(float(mass @ ((shape_score - mass @ shape_score) * cost)))

    coordinate_covariances = np.asarray(coordinate_covariances)
    coordinate_energies = np.asarray(coordinate_energies)
    mixture_weights = np.concatenate(mixture_weights)
    mixture_scores = np.concatenate(mixture_scores)
    mixture_costs = np.concatenate(mixture_costs)

    total_context_cov = float(np.sum(coordinate_covariances))
    global_cov = float(mu @ ((H.g - mean_g) * load))
    mix_mean_score = float(mixture_weights @ mixture_scores)
    mix_mean_cost = float(mixture_weights @ mixture_costs)
    mixture_cov = float(
        mixture_weights
        @ ((mixture_scores - mix_mean_score) * (mixture_costs - mix_mean_cost))
    )

    # Independent-copy symmetrization of the coordinate-context mixture.
    pair_cov = 0.5 * float(
        np.sum(
            mixture_weights[:, None]
            * mixture_weights[None, :]
            * (mixture_scores[:, None] - mixture_scores[None, :])
            * (mixture_costs[:, None] - mixture_costs[None, :])
        )
    )

    normalized_B = H.f[None, :] * H.excl
    fixed_size_error = float(
        np.max(np.abs(np.sum(normalized_B, axis=0) - (H.n - H.m) * H.f))
    )

    errors = {
        "state_load": abs(total_context_cov - global_cov),
        "coordinate_mixture": abs(total_context_cov - (H.n - 1) * mixture_cov),
        "pair_symmetrization": abs(mixture_cov - pair_cov),
        "matched_split": float(
            np.max(
                np.abs(
                    coordinate_covariances - np.asarray(split_b) - np.asarray(split_shape)
                )
            )
        ),
        "fixed_size_external_sum": fixed_size_error,
        "mean_score": abs(mix_mean_score - mean_g),
    }
    assert max(errors.values()) < 4e-8, errors
    return {
        "mu": mu,
        "load": load,
        "covariances": coordinate_covariances,
        "energy": coordinate_energies,
        "total_cov": total_context_cov,
        "global_cov": global_cov,
        "mixture_cov": mixture_cov,
        "errors": errors,
    }


def integrated_audit(H: Endpoint, order: int = 96) -> dict[str, float]:
    nodes, weights = leggauss(order)
    weighted_energy = 0.0
    signed_covariance = 0.0
    coordinatewise_adverse = 0.0
    aggregate_instantaneous_adverse = 0.0
    max_error = 0.0
    for node, weight0 in zip(nodes, weights):
        t = float((node + 1.0) / 2.0)
        weight = float(weight0 / 2.0)
        z = vertex_state(H, t)
        weighted_energy += weight * t * float(np.sum(z["energy"]))
        signed_covariance += weight * float(z["total_cov"])
        coordinatewise_adverse += weight * float(np.sum(np.maximum(-z["covariances"], 0.0)))
        aggregate_instantaneous_adverse += weight * max(-float(z["total_cov"]), 0.0)
        max_error = max(max_error, max(z["errors"].values()))

    endpoint = vertex_state(H, 1.0)
    mu1 = endpoint["mu"]
    endpoint_cost = float(mu1 @ endpoint["load"])
    boundary_error = abs(weighted_energy - (endpoint_cost - signed_covariance))
    assert boundary_error < 5e-8, boundary_error
    return {
        "weighted_energy": weighted_energy,
        "endpoint_cost": endpoint_cost,
        "signed_covariance": signed_covariance,
        "negative_signed_integral": max(-signed_covariance, 0.0),
        "aggregate_instantaneous_adverse": aggregate_instantaneous_adverse,
        "coordinatewise_adverse": coordinatewise_adverse,
        "boundary_error": boundary_error,
        "max_identity_error": max_error,
    }


def main() -> None:
    cases = [
        ("A4", A4, 0.5, 3),
        ("A6", A6, 0.5, 3),
        ("A8", A8, 2.0, 5),
        ("A9", A9, 2.0, 7),
    ]
    for name, matrix, beta, m in cases:
        result = integrated_audit(Endpoint(matrix, beta, m))
        print(name, "beta", beta, "m", m, result)
    print("PASS combined_covariance_r38_check")


if __name__ == "__main__":
    main()
