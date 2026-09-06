"""Exact finite-product stress test of stable noise with a high-influence W.

Z is a normalized sum of disjoint degree-five Boolean monomials.
W is a fixed linear encoding of the first five spins, and is not Gaussian.
The characteristic function is evaluated exactly by 32 seed patterns and
an independent-block cosine product; only its floating evaluation is used.
"""

import itertools
import json

import numpy as np


def check(groups, t=0.8, u=8.0):
    local = np.array(list(itertools.product((-1.0, 1.0), repeat=5)))
    weights = 2.0 ** (-np.arange(1, 6))
    coherent = local @ weights
    monomial = np.prod(local, axis=1)
    coherent_characteristic = np.mean(np.exp(1j * u * coherent))
    actual = np.cos(t / np.sqrt(groups)) ** (groups - 1) * np.mean(
        np.exp(1j * (t * monomial / np.sqrt(groups) + u * coherent))
    )
    comparison = np.exp(-t**2 / 2) * coherent_characteristic
    mixed_gamma = np.sum(
        weights[None, :] * monomial[:, None] * local, axis=1
    ) / (5 * np.sqrt(groups))
    exact_gamma_l2 = np.linalg.norm(weights) / (5 * np.sqrt(groups))
    assert abs(np.sqrt(np.mean(mixed_gamma**2)) - exact_gamma_l2) < 1e-13
    return {"groups": groups, "seed_count": 5 * groups,
            "noise_degree": 5, "coherent_degree": 1,
            "coherent_maximum_influence": 0.25,
            "noise_proper_cut_norm": 1 / np.sqrt(groups),
            "mixed_Gamma_L2": exact_gamma_l2,
            "characteristic_error": float(abs(actual - comparison)),
            "error_times_sqrt_groups": float(abs(actual - comparison) * np.sqrt(groups))}


if __name__ == "__main__":
    print(json.dumps({"status": "exact finite-block formulas, floating evaluation; W remains a fixed 32-atom law",
                      "checks": [check(m) for m in (1, 4, 16, 64, 256, 1024, 4096)]}, indent=2))
