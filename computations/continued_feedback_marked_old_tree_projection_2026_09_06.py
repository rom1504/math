#!/usr/bin/env python3
"""Falsifier for a proposed energy projection involving genuine old history.

The actual degree-three marked input is D=S*h2(BS), its old field is
Y=B D, and its literal return is V=Q D. No theorem is asserted here.
"""

from __future__ import annotations

import argparse
import json

import numpy as np
from numpy.polynomial.hermite import hermgauss

from continued_feedback_coherent_energy_projection_2026_09_06 import signing
from continued_director_steiner_signing_2026_09_06 import make_signing
from continued_feedback_variance_falsifier_2026_09_06 import SEED


def run(kind: str, size: int, samples: int, batch: int, seed: int, mask: bool, angle: float) -> dict:
    rng = np.random.default_rng(seed)
    if kind in ("steiner_seed", "steiner_twins"):
        small = SEED if kind == "steiner_seed" else np.ones((2, 2))
        a = np.kron(small, make_signing(size)[0]).astype(float)
        # The base Steiner diagonal is zero: replace all remaining zero
        # off-diagonal entries by signs to obtain an actual dense signing.
        a[a == 0] = 1
        np.fill_diagonal(a, 0)
    else:
        a = signing(kind, size, rng)
    n = len(a)
    b = a / np.sqrt(n - 1)
    q = b @ b
    first = np.exp(-0.5)
    t = b @ (np.exp(-1) * (np.sinh(q) - q)) @ b
    sigma = np.sqrt(np.maximum(np.diag(t), 0))
    nodes, weights = hermgauss(32)
    nodes *= np.sqrt(2)
    weights /= np.sqrt(np.pi)

    def fields(count: int):
        s = rng.choice((-1.0, 1.0), size=(count, n))
        g = s @ b
        d = s * (g**2 - 1) / np.sqrt(2)
        old = d @ b
        v = np.cos(angle) * (s @ q) + np.sin(angle) * (old @ b)
        score = np.cos(angle) * g + np.sin(angle) * old
        f = np.sin(score)
        h = 1 - np.abs(f) if mask else np.ones_like(f)
        feedback = f @ b
        z = feedback - first * v
        c = h * np.tanh(feedback)
        c0 = np.zeros_like(c)
        derivative = np.zeros_like(c)
        for node, weight in zip(nodes, weights):
            response = np.tanh(first * v + node * sigma)
            c0 += weight * h * response
            derivative += weight * h * (1 - response**2)
        return c, c0, z, derivative, np.cos(v), np.sin(v)

    coefficient = np.zeros(n)
    cosine_mean = np.zeros(n)
    done = 0
    while done < samples:
        count = min(batch, samples - done)
        pilot = fields(count)
        coefficient += pilot[3].sum(axis=0)
        cosine_mean += pilot[4].sum(axis=0)
        done += count
    coefficient /= samples
    cosine_mean /= samples
    noise_energy = np.sum(b * (coefficient[:, None] * t * coefficient[None, :])) / (2 * n)
    moments = np.zeros((5, 2))
    done = 0
    while done < samples:
        count = min(batch, samples - done)
        c, c0, z, _, cosine, sine = fields(count)
        actual = np.sum((c @ b) * c, axis=1) / (2 * n)
        coherent = np.sum((c0 @ b) * c0, axis=1) / (2 * n)
        cross = np.sum((c0 @ b) * (z * coefficient), axis=1) / n
        projected = coherent + cross + noise_energy
        centered_cross = np.sum(((cosine - cosine_mean) * z) * (sine @ b), axis=1) / n
        data = np.stack((actual, coherent, projected, actual - projected, centered_cross))
        moments[:, 0] += data.sum(axis=1)
        moments[:, 1] += (data**2).sum(axis=1)
        done += count
    means = moments[:, 0] / samples
    errors = np.sqrt(np.maximum(moments[:, 1] / samples - means**2, 0) / samples)
    return {
        "status": "Conjecture falsification only; paired SE omits pilot and quadrature errors",
        "family": kind,
        "n": n,
        "size_parameter": size,
        "seed": seed,
        "samples_pilot_and_test_each": samples,
        "score_angle": angle,
        "score": "cos(angle)*G + sin(angle)*Y",
        "mask": "1-|sin(score)|" if mask else "1",
        "actual_energy": float(means[0]),
        "literal_coherent_energy": float(means[1]),
        "proposed_projected_energy": float(means[2]),
        "paired_difference": float(means[3]),
        "paired_difference_se_conditional_on_pilot": float(errors[3]),
        "noise_trace_energy": float(noise_energy),
        "centered_coherent_coefficient_cross": float(means[4]),
        "centered_cross_se_conditional_on_pilot": float(errors[4]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("steiner", "random", "seed_lift", "apex", "steiner_seed", "steiner_twins"), default="steiner")
    parser.add_argument("--sizes", type=int, nargs="+", default=[3, 4, 5])
    parser.add_argument("--samples", type=int, default=10000)
    parser.add_argument("--batch", type=int, default=250)
    parser.add_argument("--seed", type=int, default=6092026)
    parser.add_argument("--mask", action="store_true")
    parser.add_argument("--angle", type=float, default=float(np.pi / 2))
    args = parser.parse_args()
    for size in args.sizes:
        print(json.dumps(run(args.kind, size, args.samples, args.batch, args.seed + size, args.mask, args.angle)), flush=True)


if __name__ == "__main__":
    main()
