"""Replay exact algebra and finite feedback in the tau=0 counterexample.

The infinite fixed-response construction and its asymptotic proof are in
artifacts/continued_audit_zero_variance_threshold_counterexample_2026_09_06.md.
Only bounded finite checks are performed here. Grid membership uses
integer arithmetic, not floating-point equality of quadratic irrationals.
"""

from __future__ import annotations

import argparse
import json
import math

import numpy as np

from continued_director_steiner_signing_2026_09_06 import make_signing


def squarefree(number: int) -> bool:
    return all(number % (d * d) for d in range(2, math.isqrt(number) + 1))


def selected_row_size(order: int) -> int:
    row_size = order - 1
    while not squarefree(row_size):
        row_size -= 1
    assert order - row_size <= math.ceil(4 * math.sqrt(order))
    return row_size


def grid_matches(k: int, m: int, other: int) -> bool:
    if k == 0:
        return False
    numerator = k * k * other
    if numerator % m:
        return False
    other_k_squared = numerator // m
    other_k = math.isqrt(other_k_squared)
    return (
        other_k * other_k == other_k_squared
        and 0 < other_k <= other
        and other_k % 2 == other % 2
    )


def run(size: int, selected: list[int], samples: int, batch: int, seed: int):
    full, parameters = make_signing(size)
    full_order = len(full)
    m = selected_row_size(full_order)
    n = m + 1
    a = full[:n, :n]
    removed = full_order - n
    correction = full[:n, n:] @ full[n:, :n]
    coefficient = parameters["r"] + 3
    assert np.array_equal(
        a @ a,
        (full_order - 1) * np.eye(n, dtype=np.int64)
        + coefficient * a - correction,
    )
    for k in range(-m, m + 1, 2):
        matches = [other for other in selected if grid_matches(k, m, other)]
        assert matches == ([] if k == 0 else [m])
    b = a.astype(float) / np.sqrt(m)
    q = b @ b
    off_q = q - np.eye(n)
    assert np.max(np.abs(off_q)) <= (coefficient + removed) / m + 1e-12
    c = np.exp(1.5) / 2

    def u(g):
        return (np.sin(g) - c * np.sin(2 * g)) / (1 + c)

    def covariance(correlation):
        return (
            np.exp(-1) * np.sinh(correlation)
            + c * c * np.exp(-4) * np.sinh(4 * correlation)
            - 2 * c * np.exp(-2.5) * np.sinh(2 * correlation)
        ) / (1 + c) ** 2

    tau_squared = float(covariance(1))
    t = b @ covariance(q) @ b
    sigma = np.sqrt(np.diag(t))
    ideal_energy = np.sum(b * t / sigma[:, None] / sigma[None, :]) / (np.pi * n)
    epsilon = 1 / n
    rng = np.random.default_rng(seed)
    sums = np.zeros((3, 2))
    for start in range(0, samples, batch):
        count = min(batch, samples - start)
        s = rng.choice((-1.0, 1.0), (count, n))
        integer_g = s @ a
        assert np.all(np.remainder(integer_g, 2) == m % 2)
        g = integer_g / np.sqrt(m)
        response = u(g)
        # This is exactly the fixed Thomae function on this grid, including
        # g=0 since u(0)=0; no n-dependent function is substituted in proof.
        f = epsilon * response
        h = 1 - np.abs(f)
        returned_unscaled = response @ b
        hard = h * np.sign(returned_unscaled)
        smooth = h * np.tanh(epsilon * returned_unscaled)
        assert np.max(np.abs(f + hard)) <= 1 + 1e-15
        assert np.max(np.abs(-f + hard)) <= 1 + 1e-15
        values = np.stack((
            np.sum(hard * (hard @ b), axis=1) / (2 * n),
            np.sum(smooth * (smooth @ b), axis=1) / (2 * n),
            np.mean((epsilon * returned_unscaled) ** 2, axis=1),
        ))
        sums[:, 0] += values.sum(axis=1)
        sums[:, 1] += (values * values).sum(axis=1)
    means = sums[:, 0] / samples
    errors = np.sqrt(np.maximum(sums[:, 1] / samples - means**2, 0) / samples)
    return {
        "full_Steiner_order": full_order,
        "n": n,
        "squarefree_row_size": m,
        "deleted_vertices": removed,
        "exact_compression_and_disjoint_grid_checks": "PASS",
        "max_offdiagonal_Q": float(np.max(np.abs(off_q))),
        "normalized_cubic_trace": float(np.sum(b * q) / n),
        "genuine_u_gaussian_variance": tau_squared,
        "fixed_Thomae_f_gaussian_variance": 0,
        "T_minus_tau_squared_Q_frobenius": float(np.linalg.norm(t - tau_squared * q)),
        "positive_variance_u_ideal_energy": float(ideal_energy),
        "asymptotic_hard_energy_proved": 1 / (np.pi * np.sqrt(2)),
        "samples": samples,
        "finite_results": {
            name: {"mean": float(mean), "se": float(error)}
            for name, mean, error in zip(
                ("hard_energy", "smooth_energy", "raw_return_variance"), means, errors
            )
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", type=int, nargs="+", default=[3, 4, 5])
    parser.add_argument("--samples", type=int, default=6000)
    parser.add_argument("--batch", type=int, default=200)
    parser.add_argument("--seed", type=int, default=9061022)
    args = parser.parse_args()
    selected = [selected_row_size((2**size - 1) * 2 ** (size - 1)) for size in args.sizes]
    assert len(set(selected)) == len(selected)
    for size in args.sizes:
        print(json.dumps(run(size, selected, args.samples, args.batch, args.seed + size)), flush=True)


if __name__ == "__main__":
    main()
