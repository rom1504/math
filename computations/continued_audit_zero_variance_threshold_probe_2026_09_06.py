"""Actual-signing probe of the positive-variance hard-threshold hypothesis.

Delete one vertex from a binary Steiner signing. At the resulting odd
order, G_i=0 can occur. The fixed regular function
f(g,y)=1_{g=0} sign(y) is Gaussian-a.e. zero, but not identically zero on
the actual discrete fields. H=1-|f| yields feasible cube endpoints.
No asymptotic theorem is claimed by this Monte Carlo test.
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from continued_director_steiner_signing_2026_09_06 import make_signing


def run(size: int, samples: int, batch: int, seed: int):
    complete, _ = make_signing(size)
    a = complete[:-1, :-1].astype(float)
    n = len(a)
    b = a / np.sqrt(n - 1)
    assert n % 2 == 1
    rng = np.random.default_rng(seed)
    sums = np.zeros((5, 2))
    for start in range(0, samples, batch):
        count = min(batch, samples - start)
        s = rng.choice((-1.0, 1.0), (count, n))
        integer_field = s @ a
        g = integer_field / np.sqrt(n - 1)
        d = s * (g**2 - 1) / np.sqrt(2)
        y = d @ b
        f = (integer_field == 0) * np.sign(y)
        h = 1 - np.abs(f)
        returned = f @ b
        c = h * np.sign(returned)
        smooth = h * np.tanh(returned)
        values = np.stack((
            np.mean(f**2, axis=1),
            np.sum((c @ b) * c, axis=1) / (2 * n),
            np.sum((smooth @ b) * smooth, axis=1) / (2 * n),
            np.mean(returned**2, axis=1),
            np.mean(returned == 0, axis=1),
        ))
        assert np.max(np.abs(f + c)) <= 1
        assert np.max(np.abs(-f + c)) <= 1
        sums[:, 0] += values.sum(axis=1)
        sums[:, 1] += (values**2).sum(axis=1)
    means = sums[:, 0] / samples
    errors = np.sqrt(np.maximum(sums[:, 1] / samples - means**2, 0) / samples)
    names = (
        "raw_f_squared_norm_per_root",
        "hard_feedback_energy",
        "smooth_feedback_energy",
        "raw_return_variance",
        "exact_zero_return_fraction",
    )
    return {
        "family": "one_vertex_deleted_Steiner",
        "n": n,
        "size_parameter": size,
        "samples": samples,
        "seed": seed,
        "gaussian_f_residual_variance_exact": 0,
        "values": {
            name: {"mean": float(mean), "se": float(error)}
            for name, mean, error in zip(names, means, errors)
        },
        "scope": "actual-signing finite stress test, not an asymptotic proof",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", type=int, nargs="+", default=[3, 4, 5])
    parser.add_argument("--samples", type=int, default=6000)
    parser.add_argument("--batch", type=int, default=200)
    parser.add_argument("--seed", type=int, default=9061013)
    args = parser.parse_args()
    for size in args.sizes:
        print(json.dumps(run(size, args.samples, args.batch, args.seed + size)), flush=True)


if __name__ == "__main__":
    main()
