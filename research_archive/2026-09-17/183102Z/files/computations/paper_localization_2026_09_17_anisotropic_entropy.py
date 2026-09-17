#!/usr/bin/env python3
"""Finite channel audit for flat-mode anisotropic Boolean entropy.

All Boolean words are enumerated in the small cases. Mode-cutoff membership,
Gaussian capacity, and normal tails are evaluated numerically; none is a
substitute for the analytic finite theorem.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import norm


def hadamard(k: int) -> np.ndarray:
    h = np.ones((1, 1), dtype=np.int64)
    while len(h) < k:
        h = np.block([[h, h], [h, -h]])
    return h


def entropy(t: float) -> float:
    if t <= 0 or t >= 1:
        return 0.0
    return -t * math.log(t) - (1 - t) * math.log1p(-t)


def finite_case(k: int, p: int, theta: float, zeta: Fraction,
                delta: float) -> dict:
    n, h = k * p, hadamard(k)
    words = np.array(list(itertools.product((-1, 1), repeat=n)), dtype=np.int64)
    transformed = np.einsum("ab,ibp->iap", h, words.reshape(-1, k, p))
    energies = np.sum(transformed ** 2, axis=2)
    effective = np.sqrt(energies).sum(axis=1) ** 2 / (k * n)
    keep = effective / k <= theta + 1e-12
    words, energies = words[keep], energies[keep]
    assert len(words) > 0
    labels = energies * zeta.denominator // (n * zeta.numerator)
    groups = defaultdict(list)
    for i, label in enumerate(labels):
        groups[tuple(int(a) for a in label)].append(i)
    max_histogram_sum = k * zeta.denominator // zeta.numerator
    label_count_upper = math.comb(max_histogram_sum + k, k)
    assert len(groups) <= label_count_upper
    projections = [np.kron(np.outer(row, row) / k, np.eye(p)) for row in h]
    histogram_entropy, average_capacity, average_trace = 0.0, 0.0, 0.0
    max_diagonal_error, max_trace_bound_excess = 0.0, -math.inf
    for label, ids in groups.items():
        probability = len(ids) / len(words)
        histogram_entropy -= probability * math.log(probability)
        upper = (np.array(label, dtype=float) + 1) * float(zeta) / k
        root_sum = float(np.sqrt(upper).sum())
        variances = delta * k * np.sqrt(upper) / root_sum
        gamma = sum((variances[a] * projections[a] for a in range(k)),
                    np.zeros((n, n)))
        diagonal_error = float(np.abs(np.diag(gamma) - delta).max())
        assert diagonal_error < 1e-12
        max_diagonal_error = max(max_diagonal_error, diagonal_error)
        samples = words[ids].astype(float)
        mean = samples.mean(axis=0)
        covariance = samples.T @ samples / len(samples) - np.outer(mean, mean)
        sign1, logdet1 = np.linalg.slogdet(gamma + covariance)
        sign0, logdet0 = np.linalg.slogdet(gamma)
        assert sign1 > 0 and sign0 > 0
        capacity = float((logdet1 - logdet0) / 2)
        trace = float(np.trace(np.linalg.solve(gamma, covariance)) / 2)
        trace_bound = n * root_sum ** 2 / (2 * delta * k)
        assert capacity <= trace + 1e-9
        assert trace <= trace_bound + 1e-9
        max_trace_bound_excess = max(max_trace_bound_excess, trace - trace_bound)
        average_capacity += probability * capacity
        average_trace += probability * trace
    error = float(norm.cdf(-1 / math.sqrt(delta)))
    actual_entropy = math.log(len(words))
    capacity_decoder_bound = histogram_entropy + average_capacity + n * entropy(error)
    finite_theorem_bound = (math.log(label_count_upper)
                            + n * (math.sqrt(theta) + math.sqrt(float(zeta))) ** 2 / (2 * delta)
                            + n * entropy(error))
    assert actual_entropy <= capacity_decoder_bound + 1e-8
    assert actual_entropy <= finite_theorem_bound + 1e-8
    return {"k": k, "p": p, "n": n, "theta": theta,
            "zeta": str(zeta), "delta": delta, "code_size": len(words),
            "realized_histograms": len(groups), "histogram_count_upper_bound": label_count_upper,
            "actual_code_entropy": actual_entropy,
            "histogram_entropy": histogram_entropy,
            "average_gaussian_capacity_upper_bound": average_capacity,
            "average_trace_information_upper_bound": average_trace,
            "exact_marginal_decoder_error": error,
            "capacity_plus_decoder_entropy_bound": capacity_decoder_bound,
            "finite_theorem_bound": finite_theorem_bound,
            "maximum_noise_diagonal_error": max_diagonal_error,
            "maximum_trace_bound_excess": max_trace_bound_excess}


def envelope(theta: float) -> dict:
    ll = math.log(1 / theta)
    objective = lambda t: theta * t + entropy(float(norm.cdf(-math.sqrt(2 * t))))
    result = minimize_scalar(objective, bounds=(1e-10, 2 * ll + 30), method="bounded",
                             options={"xatol": 1e-10})
    explicit_t = ll + 2 * math.log(ll)
    explicit_bound = theta * explicit_t + theta / ll ** 2 * (explicit_t + 1)
    assert result.fun <= explicit_bound + 1e-12
    return {"theta": theta, "optimizing_t_numeric": float(result.x),
            "Psi_numeric": float(result.fun),
            "ratio_to_theta_log_inverse": float(result.fun / (theta * ll)),
            "explicit_analytic_upper_bound": explicit_bound}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=Path("tmp/paper_portfolio_2026_09_17/localization/anisotropic_entropy_audit.json"))
    args = parser.parse_args()
    cases = [finite_case(k, p, theta, zeta, delta)
             for k, p, theta in [(4, 2, .25), (4, 2, .5), (4, 2, 1.0),
                                 (4, 3, .25), (4, 3, .5), (4, 3, .75)]
             for zeta, delta in [(Fraction(1, 5), .2), (Fraction(1, 2), .7)]]
    result = {"status": "PASS", "finite_channel_cases": cases,
              "envelope_diagnostics": [envelope(10.0 ** (-j)) for j in [2, 3, 4, 6, 8, 10, 12]],
              "scope": "Finite numerical capacity/decoder diagnostics; histogram and entropy proof are analytic."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "finite_channel_cases": len(cases),
                      "max_noise_diagonal_error": max(v["maximum_noise_diagonal_error"] for v in cases),
                      "envelope_diagnostics": result["envelope_diagnostics"]}, indent=2))


if __name__ == "__main__":
    main()
