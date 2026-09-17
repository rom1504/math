"""Finite replay for actual hot quadratic feature laws; no asymptotic proof claim."""
import itertools
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.hermite import hermgauss


def cube(n):
    return np.asarray(list(itertools.product((-1.0, 1.0), repeat=n)))


def hadamard(n):
    h = np.ones((1, 1))
    while h.shape[0] < n:
        h = np.block([[h, h], [h, -h]])
    assert h.shape[0] == n
    return h


def replay_frame(u, v=2.0, order=40):
    n, r = u.shape
    assert np.max(np.abs(u.T @ u - np.eye(r))) < 1e-12
    h = cube(n)
    a = 1 - 1 / v
    feature = h @ u
    energy = np.sum(feature * feature, axis=1)
    density = np.exp(a * energy / 2)
    z = float(np.mean(density))
    probability = density / np.sum(density)
    cov = h.T @ (probability[:, None] * h)
    assert np.max(np.abs(probability @ h)) < 1e-12
    assert np.max(np.abs(np.diag(cov) - 1)) < 1e-12

    nodes, weights = hermgauss(order)
    indices = np.asarray(list(itertools.product(range(order), repeat=r)))
    g = np.sqrt(2 * v) * nodes[indices]
    gh_weights = np.prod(weights[indices] / math.sqrt(math.pi), axis=1)
    latent = math.sqrt(a) * g @ u.T
    logcosh = np.logaddexp(latent, -latent) - math.log(2)
    remainder = np.sum(latent * latent / 2 - logcosh, axis=1)
    assert np.min(remainder) > -1e-10
    assert np.max(remainder - np.sum(latent ** 4, axis=1) / 12) < 1e-7
    tilt = np.exp(-remainder)
    c = float(gh_weights @ tilt)
    latent_probability = gh_weights * tilt / c
    means = np.tanh(latent)
    latent_cov = means.T @ (latent_probability[:, None] * means)
    np.fill_diagonal(latent_cov, 1.0)
    partition_error = abs(z / (v ** (r / 2)) - c)
    covariance_identity_error = float(np.max(np.abs(cov - latent_cov)))
    assert partition_error < 2e-7, partition_error
    assert covariance_identity_error < 2e-7, covariance_identity_error

    projection = u @ u.T
    leverage = np.diag(projection)
    leading_cov = np.eye(n) + (v - 1) * (
        projection - np.diag(leverage)
    )
    delta2 = float(leverage @ leverage)
    assert 1 - z / (v ** (r / 2)) <= a * a * v * v * delta2 / 4 + 1e-12
    feature_second_moment = float(probability @ energy)
    assert feature_second_moment <= v * r + 1e-12
    divergence = a * feature_second_moment / 2 - math.log(z)
    entropy_reference = r * (v - 1 - math.log(v)) / 2
    # Exhaustive query check only at these small dimensions.
    responses = np.abs(h @ h.T) @ probability / math.sqrt(n)
    leading = math.sqrt(2 / math.pi) * np.sqrt(1 + (v - 1) * energy / n)
    return {
        "n": n, "r": r, "v": v, "delta2": delta2,
        "relative_partition": z / (v ** (r / 2)),
        "quadrature_partition_error": partition_error,
        "quadrature_covariance_error": covariance_identity_error,
        "covariance_frobenius_error": float(np.linalg.norm(cov - leading_cov)),
        "covariance_operator_error": float(np.linalg.norm(cov - leading_cov, ord=2)),
        "maximum_leverage": float(np.max(leverage)),
        "feature_second_moment": feature_second_moment,
        "exact_hot_upper_second_moment": v * r,
        "physical_relative_entropy": divergence,
        "gaussian_entropy_reference": entropy_reference,
        "all_query_response_max_error": float(np.max(np.abs(responses - leading))),
        "queries": len(h),
    }


def rank_one_sequence(n, v=2.0):
    # Constant-leverage rank one: exact binomial sums, no cube enumeration.
    a = 1 - 1 / v
    sums = np.arange(n + 1, dtype=float) * 2 - n
    log_mass = np.asarray([
        math.lgamma(n + 1) - math.lgamma(j + 1) - math.lgamma(n - j + 1)
        - n * math.log(2) for j in range(n + 1)
    ])
    weighted = np.exp(log_mass + a * sums * sums / (2 * n))
    z = float(np.sum(weighted))
    law = weighted / z
    response = float(law @ np.abs(sums) / math.sqrt(n))
    offdiag = float((law @ (sums * sums) - n) / (n * (n - 1)))
    return {
        "n": n, "relative_partition": z / math.sqrt(v),
        "aligned_response": response,
        "aligned_response_limit": math.sqrt(2 * v / math.pi),
        "covariance_frobenius_error": math.sqrt(n * (n - 1))
        * abs(offdiag - (v - 1) / n),
    }


def main():
    rng = np.random.RandomState(91726)
    frames = [hadamard(8)[:, :r] / math.sqrt(8) for r in (1, 2, 3)]
    for n, r in ((6, 1), (8, 2), (10, 3)):
        q, _ = np.linalg.qr(rng.normal(size=(n, r)))
        frames.append(q)
    results = [replay_frame(u) for u in frames]
    sequences = [rank_one_sequence(n) for n in (16, 32, 64, 128, 256, 512, 1024)]
    out = {
        "status": "PASS",
        "scope": "Floating exact-cube sums and Gaussian quadrature; proof is in the artifact.",
        "frames": results, "rank_one_sequence": sequences,
        "total_complete_queries": sum(item["queries"] for item in results),
    }
    target = Path("tmp/paper_portfolio_2026_09_17/bernoulli/hot_feature_laws.json")
    target.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
