"""Numerical replay of the stochastic sine-rounding approximation obstruction.

The proof is analytic; Monte Carlo outputs here are sanity checks only.
Run with the repository-local .venv/bin/python. No files are written.
"""

from __future__ import annotations

import argparse
import json
import math

import numpy as np


def normalized_walsh(a):
    out = np.array(a, dtype=float, copy=True)
    n = out.shape[-1]
    assert n >= 1 and n & (n - 1) == 0
    step = 1
    while step < n:
        blocks = out.reshape(*out.shape[:-1], -1, 2 * step)
        left = blocks[..., :step].copy()
        right = blocks[..., step:].copy()
        blocks[..., :step] = left + right
        blocks[..., step:] = left - right
        step *= 2
    return out / math.sqrt(n)


def sample(n, samples, batch, seed, k):
    rng = np.random.default_rng(seed + n)
    correlated = []
    uniform = []
    for start in range(0, samples, batch):
        size = min(batch, samples - start)
        g = rng.standard_normal((size, n))
        h = normalized_walsh(g)
        assert np.allclose(np.sum(g * g, axis=1), np.sum(h * h, axis=1))
        x = 2 * (rng.random(g.shape) < (1 + np.sin(g)) / 2) - 1
        y = 2 * (rng.random(g.shape) < (1 + np.sin(h)) / 2) - 1
        f = np.sum(x * normalized_walsh(y), axis=1) / n
        assert np.max(np.abs(f)) <= 1 + 1e-12
        correlated.extend(f.tolist())
        xu = rng.choice([-1.0, 1.0], size=g.shape)
        yu = rng.choice([-1.0, 1.0], size=g.shape)
        fu = np.sum(xu * normalized_walsh(yu), axis=1) / n
        uniform.extend(fu.tolist())
    correlated = np.asarray(correlated)
    uniform = np.asarray(uniform)
    p_corr = 1 - 2 * (1 - correlated**2) ** k
    p_unif = 1 - 2 * (1 - uniform**2) ** k
    return {
        "N": n, "samples": samples, "seed": seed + n, "k": k,
        "correlated_mean": float(correlated.mean()),
        "exact_correlated_mean": math.exp(-1) * math.sqrt(n) * math.sinh(1 / math.sqrt(n)),
        "N_times_correlated_variance": float(n * correlated.var()),
        "proved_upper_bound_N_times_correlated_variance": 7,
        "uniform_mean": float(uniform.mean()),
        "N_times_uniform_second_moment": float(n * np.mean(uniform**2)),
        "amplified_expectation_gap": float(p_corr.mean() - p_unif.mean()),
        "limiting_amplified_gap": 2 * (1 - (1 - math.exp(-2)) ** k),
        "status": "Monte Carlo sanity check only; analytic proof in artifact",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sizes", default="64,256,1024")
    parser.add_argument("--samples", type=int, default=8192)
    parser.add_argument("--batch", type=int, default=256)
    parser.add_argument("--seed", type=int, default=20260918)
    parser.add_argument("--k", type=int, default=8)
    args = parser.parse_args()
    for n in map(int, args.sizes.split(",")):
        print(json.dumps(sample(n, args.samples, args.batch, args.seed, args.k)), flush=True)


if __name__ == "__main__":
    main()
