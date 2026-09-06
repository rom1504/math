"""Replay the actual-Steiner nonzero-first covariance obstruction.

Integer signing identities are exact. Displayed Gaussian/binomial quantities
are floating diagnostics; the artifact contains the rational asymptotic gap.
No random sampling, solver, or stored-output overwrite is needed.
"""

import json
import math

import numpy as np
from scipy.special import ndtr
from scipy.stats import binom

from continued_director_steiner_signing_2026_09_06 import make_signing


def main():
    cutoff = 4.0
    t = math.pi / 4
    phi = math.exp(-cutoff**2 / 2) / math.sqrt(2 * math.pi)
    delta = math.sqrt(2 * ((1 + cutoff**2) * ndtr(-cutoff) - cutoff * phi))
    records = []
    for exponent in [2, 3, 4, 5]:
        a, record = make_signing(exponent)
        n = len(a)
        m = n - 1
        b = a / math.sqrt(m)
        q = (a @ a) / m
        gamma = (record["r"] + 3) / math.sqrt(m)
        assert np.max(np.abs(q - np.eye(n) - gamma * b)) < 1e-12
        # Exact Rademacher characteristic-function identity, evaluated in float.
        actual_linear_diagonal = 0.5 * (1 - np.prod(np.cos(2 * t * q), axis=1))
        assert np.max(np.abs(actual_linear_diagonal - 0.5)) < 1e-14
        gaussian_diagonal = 0.5 * (1 - math.exp(-2 * t*t * (1 + gamma*gamma)))
        k = np.arange(m + 1)
        row_values = (2 * k - m) / math.sqrt(m)
        excess = np.maximum(np.abs(row_values) - cutoff, 0)
        delta_n = math.sqrt(float(np.dot(binom.pmf(k, m, 0.5), excess**2)))
        op = record["positive_eigenvalue"] / math.sqrt(m)
        bounded_gap_lower = 0.25 * (0.5 - gaussian_diagonal - t * op * (delta_n + delta))
        records.append({
            "order": n,
            "gamma_squared": gamma*gamma,
            "actual_linear_sine_squared_diagonal": float(np.mean(actual_linear_diagonal)),
            "matched_gaussian_sine_squared_diagonal": gaussian_diagonal,
            "finite_binomial_clipping_l2": delta_n,
            "bounded_feasible_covariance_trace_gap_lower_float": bounded_gap_lower,
        })
    print(json.dumps({
        "normal_clipping_l2": delta,
        "asymptotic_unclipped_gap": 0.5 * math.exp(-3 * math.pi**2 / 16),
        "asymptotic_bounded_gap_lower_using_cap_2": 0.25 * (0.5 * math.exp(-3 * math.pi**2 / 16) - math.pi * delta),
        "coarse_proved_rational_lower_bound": "3/256",
        "sign_preserving_variant_rational_lower_bound": "95/8192",
        "records": records,
    }, indent=2))


if __name__ == "__main__":
    main()
