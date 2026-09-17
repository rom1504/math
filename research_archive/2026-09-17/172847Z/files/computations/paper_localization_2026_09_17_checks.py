#!/usr/bin/env python3
"""Reproducible finite checks for the localization transfer obstruction.

The analytic theorem does not depend on these numerical checks. Exact small
enumerations check switching invariance and a family of finite observation
channels. Random posterior checks stress the pointwise entropy inequality.
Gaussian posterior sampling checks the primary decomposition's normalization.

Run with the repository virtual environment. All generated output is retained
under tmp/paper_portfolio_2026_09_17/localization by default.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


def cube(n: int) -> np.ndarray:
    return np.array(list(itertools.product((-1, 1), repeat=n)), dtype=np.int64)


def sylvester(n: int) -> np.ndarray:
    result = np.ones((1, 1), dtype=np.int64)
    while result.shape[0] < n:
        result = np.block([[result, result], [result, -result]])
    if result.shape[0] != n:
        raise ValueError("order must be a power of two")
    return result


def cap(a: np.ndarray, signs: np.ndarray) -> int:
    return int(np.max(np.abs(np.einsum("bi,ij,bj->b", signs, a, signs) // 2)))


def phi(r: np.ndarray) -> np.ndarray:
    """Mutual information of two unbiased signs, safely at endpoints."""
    out = np.zeros_like(r, dtype=float)
    for t in (1 + r, 1 - r):
        positive = t > 0
        out[positive] += t[positive] * np.log(t[positive]) / 2
    return out


def posterior_checks(rng: np.random.Generator, n: int, samples: int) -> dict:
    # Global-flip classes have a unique representative with first spin +1.
    representatives = np.column_stack((np.ones(2 ** (n - 1), dtype=np.int64), cube(n - 1)))
    edges = np.triu_indices(n, 1)
    feature = representatives[:, edges[0]] * representatives[:, edges[1]]
    worst_ratio = 0.0
    smallest_slack = math.inf
    smallest_row_slack = math.inf
    records = []
    for trial in range(samples):
        concentration = 10 ** rng.uniform(-2, 2)
        p = rng.dirichlet(np.full(len(representatives), concentration))
        positive = p > 0
        divergence = math.log(len(p)) + float(np.sum(p[positive] * np.log(p[positive])))
        correlations = p @ feature
        squared_bias = float(np.dot(correlations, correlations))
        ratio = squared_bias / (n * divergence) if divergence > 1e-12 else 0.0
        worst_ratio = max(worst_ratio, ratio)
        smallest_slack = min(smallest_slack, n * divergence - squared_bias)
        mutual_info = phi(correlations)
        row_sums = np.zeros(n)
        np.add.at(row_sums, edges[0], mutual_info)
        np.add.at(row_sums, edges[1], mutual_info)
        smallest_row_slack = min(smallest_row_slack, divergence - float(np.max(row_sums)))
        if trial < 3:
            records.append({"divergence": divergence, "squared_bias": squared_bias, "ratio": ratio})

    # These observations reveal relative spins xi_j*xi_1 for j<=k+1.
    revelation = []
    for k in range(n):
        information = k * math.log(2)
        squared_bias = k * (k + 1) / 2
        revelation.append({"revealed_relative_spins": k, "information_nats": information,
                           "squared_bias_exact": squared_bias,
                           "bound": n * information})
        assert squared_bias <= n * information + 1e-10
    assert smallest_slack >= -1e-9
    assert smallest_row_slack >= -1e-9
    return {"n": n, "random_posteriors": samples,
            "max_bias_to_bound_ratio": worst_ratio,
            "min_entropy_inequality_slack": smallest_slack,
            "min_row_mutual_information_slack": smallest_row_slack,
            "sample_records": records, "exact_revelation_channels": revelation}


def biased_absolute_checks(rng: np.random.Generator, max_k: int, samples: int) -> dict:
    minimum_slack = math.inf
    worst = None
    kappa = math.sqrt(2 / math.pi)
    for k in range(1, max_k + 1):
        signs = cube(k)
        for _ in range(samples):
            b = rng.uniform(-1, 1, k)
            shift = float(rng.normal() * math.sqrt(k))
            probabilities = np.prod((1 + signs * b) / 2, axis=1)
            # Center the Bernoulli sum and apply an arbitrary external mean.
            values = signs.sum(axis=1) - b.sum() + shift
            expectation = float(probabilities @ np.abs(values))
            variance = float(np.sum(1 - b * b))
            rhs = kappa * math.sqrt(variance) - 2 * k ** (1 / 3)
            slack = expectation - rhs
            if slack < minimum_slack:
                minimum_slack = slack
                worst = {"k": k, "bias": b.tolist(), "shift": shift,
                         "exact_absolute_expectation": expectation, "lower_bound": rhs}
    assert minimum_slack >= -1e-9
    u = np.linspace(0, 5, 100001)
    max_third = float(np.max(3 * u / (1 + u * u) ** 2.5))
    assert max_third < 1
    return {"max_k": max_k, "samples_per_k": samples,
            "minimum_uniform_greedy_row_slack": minimum_slack, "worst_record": worst,
            "max_smoothed_abs_scaled_third_derivative": max_third}


def hadamard_checks(rng: np.random.Generator) -> list[dict]:
    result = []
    for n in (2, 4, 8, 16):
        h = sylvester(n)
        a = h.copy()
        np.fill_diagonal(a, 0)
        signs = cube(n)
        original_cap = cap(a, signs)
        switch_caps = []
        for _ in range(8):
            switch = rng.choice((-1, 1), n)
            switch_caps.append(cap(a * switch[:, None] * switch[None, :], signs))
        assert all(value == original_cap for value in switch_caps)
        assert np.array_equal(h @ h, n * np.eye(n, dtype=np.int64))
        assert np.trace(h) == 0
        assert original_cap <= n ** 1.5 / 2 + 1e-8
        result.append({"n": n, "exact_cap": original_cap,
                       "normalized_cap": original_cap / n ** 1.5,
                       "all_switch_caps_match": True})
    return result


def gaussian_posterior_checks(rng: np.random.Generator, samples: int) -> dict:
    n = 4
    support = cube(n).astype(float)
    a = sylvester(n).astype(float)
    np.fill_diagonal(a, 0)
    energies = np.einsum("bi,ij,bj->b", support, a, support) / 2
    weights = np.exp(0.37 * energies / math.sqrt(n))
    weights /= weights.sum()
    mean = weights @ support
    sigma = support.T @ (weights[:, None] * support) - np.outer(mean, mean)
    orthogonal, _ = np.linalg.qr(rng.normal(size=(n, n)))
    precision = orthogonal @ np.diag((0.2, 0.6, 1.1, 2.0)) @ orthogonal.T
    eigenvalues, eigenvectors = np.linalg.eigh(precision)
    square_root = (eigenvectors * np.sqrt(eigenvalues)) @ eigenvectors.T
    inverse = np.linalg.inv(precision)
    # Equivalent whitened channel sqrt(tau)*L^(1/2) X + standard Gaussian.
    ec = np.zeros((n, n))
    eclc = np.zeros((n, n))
    info = 0.0
    energy_decomposition_error = 0.0
    log_weights = np.log(weights)
    quadratic = np.einsum("bi,ij,bj->b", support, precision, support)
    for _ in range(samples):
        x = support[rng.choice(len(support), p=weights)]
        tau = rng.uniform(1, 2)
        y = math.sqrt(tau) * square_root @ x + rng.normal(size=n)
        logits = log_weights + math.sqrt(tau) * (support @ square_root @ y) - tau * quadratic / 2
        posterior = np.exp(logits - logits.max())
        posterior /= posterior.sum()
        m = posterior @ support
        c = support.T @ (posterior[:, None] * support) - np.outer(m, m)
        ec += c / samples
        eclc += c @ precision @ c / samples
        positive = posterior > 0
        info += float(np.sum(posterior[positive] * (np.log(posterior[positive]) - log_weights[positive]))) / samples
        discrepancy = float(posterior @ energies - m @ a @ m / 2 - np.trace(a @ c) / 2)
        energy_decomposition_error = max(energy_decomposition_error, abs(discrepancy))
    info_bound = np.linalg.slogdet(np.eye(n) + 2 * square_root @ sigma @ square_root)[1] / 2
    # Matrix inequalities are Monte Carlo diagnostics; errors are reported, not
    # declared exact certificates. The energy identity holds samplewise.
    return {"dimension": n, "samples": samples,
            "min_eigenvalue_Linv_minus_Ecov": float(np.linalg.eigvalsh(inverse - ec).min()),
            "min_eigenvalue_Sigma_minus_ECLC": float(np.linalg.eigvalsh(sigma - eclc).min()),
            "estimated_mutual_information": info, "information_upper_bound": float(info_bound),
            "max_samplewise_energy_identity_error": energy_decomposition_error,
            "precision": precision.tolist()}


def greedy_simulation(rng: np.random.Generator, n: int, samples: int) -> dict:
    """Numerical only; exact cap is not computed at these larger dimensions."""
    totals = []
    for _ in range(samples):
        x = np.ones(n, dtype=np.int64)
        total = 0
        for i in range(1, n):
            row = rng.choice((-1, 1), i)
            field = int(row @ x[:i])
            x[i] = 1 if field >= 0 else -1
            total += abs(field)
        totals.append(total / n ** 1.5)
    return {"n": n, "samples": samples, "normalized_greedy_mean": float(np.mean(totals)),
            "normalized_greedy_standard_error": float(np.std(totals, ddof=1) / math.sqrt(samples))}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20260917)
    parser.add_argument("--posterior-samples", type=int, default=300)
    parser.add_argument("--gaussian-samples", type=int, default=20000)
    parser.add_argument("--greedy-samples", type=int, default=200)
    parser.add_argument("--output", type=Path,
                        default=Path("tmp/paper_portfolio_2026_09_17/localization/checks.json"))
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    gamma = 2 * math.sqrt(2 / math.pi) / 3
    result = {"seed": args.seed, "gamma": gamma, "gap_above_one_half": gamma - 0.5,
              "simple_required_information_per_spin": ((gamma - 0.5) / math.sqrt(2 / math.pi)) ** 2,
              "sharper_required_information_per_spin": 0.5 * (1 - 0.5 / gamma) ** (4 / 3),
              "posterior_entropy_checks": [posterior_checks(rng, n, args.posterior_samples) for n in (4, 6, 8)],
              "biased_row_checks": biased_absolute_checks(rng, 10, 50),
              "hadamard_switching_checks": hadamard_checks(rng),
              "gaussian_localization_diagnostic": gaussian_posterior_checks(rng, args.gaussian_samples),
              "unbiased_greedy_simulation": [greedy_simulation(rng, n, args.greedy_samples) for n in (64, 256, 1024)]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"output": str(args.output), "gamma": gamma,
                      "sharper_required_information_per_spin": result["sharper_required_information_per_spin"],
                      "posterior_entropy_ratios": [r["max_bias_to_bound_ratio"] for r in result["posterior_entropy_checks"]],
                      "gaussian_localization_diagnostic": result["gaussian_localization_diagnostic"],
                      "unbiased_greedy_simulation": result["unbiased_greedy_simulation"]}, indent=2))


if __name__ == "__main__":
    main()
