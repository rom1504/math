#!/usr/bin/env python3
"""Numerical falsification test for the coherent-return energy projection.

This is NOT a certificate. The seed and Monte Carlo sample size are
printed, and paired standard errors are retained. All signings are real
finite hollow signings. The nonlinear input is f(g)=sin(g).
"""

from __future__ import annotations

import argparse
import json

import numpy as np
from numpy.polynomial.hermite import hermgauss

from continued_director_steiner_signing_2026_09_06 import make_signing
from continued_feedback_variance_falsifier_2026_09_06 import SEED, sylvester


def signing(kind: str, size: int, rng: np.random.Generator) -> np.ndarray:
    if kind == "steiner":
        return make_signing(size)[0].astype(float)
    if kind == "random":
        a = np.triu(rng.choice((-1.0, 1.0), size=(size, size)), 1)
        return a + a.T
    if kind == "seed_lift":
        a = np.kron(SEED, sylvester(size)).astype(float)
        np.fill_diagonal(a, 0)
        return a
    if kind == "apex":
        a = np.ones((2 * size + 1, 2 * size + 1))
        a[1:, 1:] = np.kron([[1, -1], [-1, 1]], sylvester(size))
        np.fill_diagonal(a, 0)
        return a
    raise ValueError(kind)


def run(kind: str, size: int, samples: int, batch: int, seed: int, mask: bool) -> dict:
    rng = np.random.default_rng(seed)
    a = signing(kind, size, rng)
    n = len(a)
    b = a / np.sqrt(n - 1)
    q = b @ b
    first = np.exp(-0.5)
    residual_covariance = np.exp(-1) * (np.sinh(q) - q)
    t = b @ residual_covariance @ b
    sigma = np.sqrt(np.maximum(np.diag(t), 0))
    nodes, weights = hermgauss(24)
    nodes = nodes * np.sqrt(2)
    weights = weights / np.sqrt(np.pi)

    def fields(count: int):
        s = rng.choice((-1.0, 1.0), size=(count, n))
        g = s @ b
        v = s @ q
        f = np.sin(g)
        h = 1 - np.abs(f) if mask else np.ones_like(f)
        y = f @ b
        z = y - first * v
        c = h * np.tanh(y)
        c0 = np.zeros_like(c)
        derivative = np.zeros_like(c)
        for node, weight in zip(nodes, weights):
            response = np.tanh(first * v + node * sigma)
            c0 += weight * h * response
            derivative += weight * h * (1 - response**2)
        return c, c0, z, derivative

    # Independent pilot estimates the deterministic coefficient a_i.
    coefficient = np.zeros(n)
    done = 0
    while done < samples:
        count = min(batch, samples - done)
        _, _, _, derivative = fields(count)
        coefficient += derivative.sum(axis=0)
        done += count
    coefficient /= samples
    noise_energy = np.sum(b * (coefficient[:, None] * t * coefficient[None, :])) / (2 * n)
    moments = np.zeros((4, 2))
    done = 0
    while done < samples:
        count = min(batch, samples - done)
        c, c0, z, _ = fields(count)
        c_energy = np.sum((c @ b) * c, axis=1) / (2 * n)
        c0_energy = np.sum((c0 @ b) * c0, axis=1) / (2 * n)
        cross = np.sum((c0 @ b) * (z * coefficient), axis=1) / n
        projected = c0_energy + cross + noise_energy
        data = np.stack((c_energy, c0_energy, projected, c_energy - projected))
        moments[:, 0] += data.sum(axis=1)
        moments[:, 1] += (data**2).sum(axis=1)
        done += count
    means = moments[:, 0] / samples
    errors = np.sqrt(np.maximum(moments[:, 1] / samples - means**2, 0) / samples)
    return {
        "status": "Monte Carlo only; pilot and Gaussian quadrature errors excluded from paired SE",
        "family": kind,
        "size_parameter": size,
        "n": n,
        "seed": seed,
        "samples_pilot_and_test_each": samples,
        "mask": "1-|sin(G)|" if mask else "1",
        "operator_norm": float(np.linalg.norm(b, 2)),
        "minimum_noise_variance": float(np.min(sigma**2)),
        "actual_energy": float(means[0]),
        "coherent_energy": float(means[1]),
        "projected_energy": float(means[2]),
        "paired_difference": float(means[3]),
        "paired_difference_se_conditional_on_pilot": float(errors[3]),
        "noise_trace_energy": float(noise_energy),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("steiner", "random", "seed_lift", "apex"), default="steiner")
    parser.add_argument("--sizes", type=int, nargs="+", default=[3, 4, 5])
    parser.add_argument("--samples", type=int, default=20000)
    parser.add_argument("--batch", type=int, default=250)
    parser.add_argument("--seed", type=int, default=17092026)
    parser.add_argument("--mask", action="store_true")
    args = parser.parse_args()
    for size in args.sizes:
        print(json.dumps(run(args.kind, size, args.samples, args.batch, args.seed + size, args.mask)), flush=True)


if __name__ == "__main__":
    main()
