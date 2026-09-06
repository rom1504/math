"""Independent complete-seed and Gaussian-quadrature check of sine formulas.

Checks exact comparison identities, not a finite error bound for feedback.
The endpoint test includes zero and negative cosine factors, for which
prefix/suffix products must not be replaced by division.
"""

import itertools
import json

import numpy as np
from numpy.polynomial.hermite import hermgauss

from continued_director_steiner_signing_2026_09_06 import make_signing


def check(a, theta_diagonal):
    n = len(a)
    bmat = a / np.sqrt(n - 1)
    q = bmat @ bmat
    kappa, height = 0.6, 0.4
    first = kappa * np.exp(-0.5)
    frequency = theta_diagonal / first
    tmat = first**2 * bmat @ (np.sinh(q) - q) @ bmat
    sigma2 = np.diag(tmat)
    assert min(sigma2) > -1e-12
    sigma = np.sqrt(np.maximum(sigma2, 0))
    damping = height * np.exp(-frequency**2 * sigma2 / 2)
    angles = frequency * first * q
    products = np.prod(np.cos(angles), axis=1)
    inserted = np.empty((n, n))
    for i in range(n):
        for k in range(n):
            inserted[i, k] = np.sin(angles[i, k]) * np.prod(
                np.cos(np.delete(angles[i], k))
            )
    projection = damping[:, None] * inserted
    derivative = frequency * damping * products
    cross_kernel = np.empty((n, n))
    for i in range(n):
        for j in range(n):
            cross_kernel[i, j] = kappa * damping[i] / 2 * (
                np.prod(np.cos(angles[i] - bmat[j]))
                - np.prod(np.cos(angles[i] + bmat[j]))
            )
    cross_kernel -= first * projection @ bmat
    predicted_gain = np.mean(damping * (
        first * np.sum(q * inserted, axis=1)
        + frequency * sigma2 * products
    ))

    seeds = np.array(list(itertools.product((-1.0, 1.0), repeat=n)))
    g, v = seeds @ bmat, seeds @ q
    coherent = damping[None, :] * np.sin(frequency * first * v)
    residual = kappa * np.sin(g) - first * g
    actual_projection = coherent.T @ seeds / len(seeds)
    actual_kernel = coherent.T @ residual / len(seeds)
    actual_derivative = frequency * damping * np.mean(
        np.cos(frequency * first * v), axis=0
    )
    exact_cross = np.mean(np.sum(
        (coherent @ bmat) * derivative[None, :] * (residual @ bmat), axis=1
    )) / n
    matrix_cross = np.sum(
        ((bmat * derivative[None, :]) @ bmat) * cross_kernel
    ) / n

    nodes, weights = hermgauss(80)
    exact_gain = 0.0
    for node, weight in zip(np.sqrt(2) * nodes, weights / np.sqrt(np.pi)):
        returned = first * v + sigma[None, :] * node
        exact_gain += weight * np.mean(
            height * returned * np.sin(frequency * returned)
        )
    errors = {
        "first_Walsh": float(np.max(np.abs(actual_projection - projection))),
        "bare_star_kernel": float(np.max(np.abs(actual_kernel - cross_kernel))),
        "mean_derivative": float(np.max(np.abs(actual_derivative - derivative))),
        "cross_energy_orientation": float(abs(exact_cross - matrix_cross)),
        "tested_gain": float(abs(exact_gain - predicted_gain)),
    }
    assert max(errors.values()) < 1e-11, errors
    return {
        "order": n,
        "all_seed_count": len(seeds),
        "tb": theta_diagonal,
        "gain": float(predicted_gain),
        "max_errors": errors,
    }


if __name__ == "__main__":
    steiner, _ = make_signing(2)
    rng = np.random.default_rng(9060859)
    upper = np.triu(rng.choice((-1, 1), size=(10, 10)), 1)
    matrices = (steiner, upper + upper.T)
    print(json.dumps({
        "status": "complete Boolean enumeration plus 80-node Gaussian quadrature; floating identity checks only",
        "checks": [check(a, theta) for a in matrices for theta in (0.4, np.pi / 2, 2.4)],
    }, indent=2))
