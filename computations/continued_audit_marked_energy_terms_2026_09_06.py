"""Independent paired stress test of every first-marked energy remainder.

Actual hollow positive-twin Steiner signings, literal V=QD, f=sin(Y),
and C=sin(B f). Conditional Gaussian sine formulas require no quadrature.
Reported standard errors are conditional on an independent coefficient
pilot. This is a falsifier search, not a proof or a finite-n certificate.
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from continued_director_steiner_signing_2026_09_06 import make_signing


def run(size: int, samples: int, batch: int, seed: int) -> dict:
    signing, _ = make_signing(size)
    a = np.kron(np.ones((2, 2)), signing).astype(float)
    a[a == 0] = 1
    np.fill_diagonal(a, 0)
    n = len(a)
    assert np.array_equal(a, a.T)
    assert np.all(np.abs(a + np.eye(n)) == 1)
    b = a / np.sqrt(n - 1)
    q = b @ b
    first = np.exp(-0.5)
    t = b @ (np.exp(-1) * (np.sinh(q) - q)) @ b
    attenuation = np.exp(-np.diag(t) / 2)
    rng = np.random.default_rng(seed)

    def fields(count: int):
        s = rng.choice((-1.0, 1.0), (count, n))
        g = s @ b
        d = s * (g * g - 1) / np.sqrt(2)
        y = d @ b
        v = y @ b
        response = np.sin(y) @ b
        z = response - first * v
        c = np.sin(response)
        c0 = attenuation * np.sin(first * v)
        coefficient = attenuation * np.cos(first * v)
        return c, c0, coefficient, z

    coefficient_mean = np.zeros(n)
    for start in range(0, samples, batch):
        count = min(batch, samples - start)
        coefficient_mean += fields(count)[2].sum(axis=0)
    coefficient_mean /= samples
    ideal_noise_energy = np.sum(
        b * t * coefficient_mean[:, None] * coefficient_mean[None, :]
    ) / (2 * n)

    names = (
        "centered_coefficient_self",
        "centered_coefficient_coherent_cross",
        "centered_coefficient_deterministic_noise_cross",
        "higher_noise_self",
        "higher_noise_coherent_cross",
        "higher_noise_full_linear_cross",
        "raw_minus_ideal_noise_trace",
        "full_projection_error",
    )
    moments = np.zeros((len(names), 2))
    algebra_error = 0.0
    for start in range(0, samples, batch):
        count = min(batch, samples - start)
        c, c0, coefficient, z = fields(count)
        deterministic_noise = coefficient_mean * z
        centered_noise = (coefficient - coefficient_mean) * z
        linear = coefficient * z
        higher = c - c0 - linear
        bc, bc0, bcentered, bhigher, bdeterministic = (
            x @ b for x in (c, c0, centered_noise, higher, deterministic_noise)
        )

        def self_energy(x, bx):
            return np.sum(x * bx, axis=1) / (2 * n)

        def cross(x, by):
            return np.sum(x * by, axis=1) / n

        projection_error = (
            self_energy(c, bc)
            - self_energy(c0, bc0)
            - cross(deterministic_noise, bc0)
            - ideal_noise_energy
        )
        values = np.stack((
            self_energy(centered_noise, bcentered),
            cross(centered_noise, bc0),
            cross(centered_noise, bdeterministic),
            self_energy(higher, bhigher),
            cross(higher, bc0),
            cross(linear, bhigher),
            self_energy(deterministic_noise, bdeterministic) - ideal_noise_energy,
            projection_error,
        ))
        algebra_error = max(
            algebra_error,
            float(np.max(np.abs(values[:-1].sum(axis=0) - values[-1]))),
        )
        moments[:, 0] += values.sum(axis=1)
        moments[:, 1] += (values * values).sum(axis=1)
    means = moments[:, 0] / samples
    errors = np.sqrt(np.maximum(moments[:, 1] / samples - means**2, 0) / samples)
    assert algebra_error < 1e-11
    return {
        "family": "actual_positive_twin_Steiner",
        "n": n,
        "size_parameter": size,
        "samples_pilot_and_test_each": samples,
        "seed": seed,
        "exact_expansion_max_error": algebra_error,
        "ideal_noise_energy": float(ideal_noise_energy),
        "terms": {
            name: {"mean": float(mean), "conditional_pilot_se": float(error)}
            for name, mean, error in zip(names, means, errors)
        },
        "scope": "Monte Carlo falsifier only; no pilot-error or finite-n asymptotic bound",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", type=int, nargs="+", default=[3, 4, 5])
    parser.add_argument("--samples", type=int, default=6000)
    parser.add_argument("--batch", type=int, default=200)
    parser.add_argument("--seed", type=int, default=9061011)
    args = parser.parse_args()
    for size in args.sizes:
        print(json.dumps(run(size, args.samples, args.batch, args.seed + size)), flush=True)


if __name__ == "__main__":
    main()
