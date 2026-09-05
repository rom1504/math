#!/usr/bin/env python3
"""Falsification laboratory for the rooted-channel Gaussian theorem.

Finite numerical evidence only. Each tested sequence has an analytic O(n^1.5)
cap bound (random families with overwhelming probability). No cap optimizer or
Gaussian limit is certified by this script. The planted clique deliberately
allows a growing normalized operator norm; the tensor family has persistent
local row correlations. Both test assumptions hidden in generic AMP claims.
"""

from __future__ import annotations

import argparse
import json
import math

import numpy as np


def regular_hadamard(power: int) -> np.ndarray:
    base = np.ones((4, 4), dtype=np.int64) - 2 * np.eye(4, dtype=np.int64)
    h = np.ones((1, 1), dtype=np.int64)
    for _ in range(power):
        h = np.kron(h, base)
    return h


def random_signing(n: int, rng: np.random.Generator) -> np.ndarray:
    a = np.triu(rng.choice((-1, 1), size=(n, n)), 1)
    return a + a.T


def models(rng: np.random.Generator, max_order: int):
    for power in (2, 3, 4, 5):
        n = 4**power
        if n > max_order:
            continue
        h = regular_hadamard(power)
        np.fill_diagonal(h, 0)
        yield "regular_hadamard", h
        a = random_signing(n, rng)
        yield "random", a
        planted = a.copy()
        clique = int(n**0.75)
        planted[:clique, :clique] = 1
        np.fill_diagonal(planted, 0)
        yield "planted_clique_n_three_quarters", planted
        gauge = a[0].copy()
        gauge[0] = 1
        yield "random_root_gauged", a * gauge[:, None] * gauge[None, :]
        if 3*n <= max_order:
            s = np.ones((3, 3), dtype=np.int64)-2*np.eye(3, dtype=np.int64)
            tensor = np.kron(regular_hadamard(power), s)
            np.fill_diagonal(tensor, 0)
            yield "hadamard_tensor_nonorthogonal_three", tensor


def test(label: str, a: np.ndarray, rng: np.random.Generator,
         samples: int, batch: int) -> dict:
    n = len(a)
    bmat = a / math.sqrt(n-1)
    t = 7/8
    b = math.erf(t/math.sqrt(2))
    aval = math.sqrt(2/math.pi)*math.exp(-t*t/2)
    selected = sorted(set((0, n//2, n-1)))
    acc = np.zeros(8, dtype=np.float64)
    row_acc = np.zeros((len(selected), 5), dtype=np.float64)
    for start in range(0, samples, batch):
        count = min(batch, samples-start)
        spins = rng.choice((-1., 1.), size=(count, n))
        fields = spins @ bmat
        mask = np.abs(fields) <= t
        v = spins * mask
        u = np.sign(fields) * (~mask)
        y = v @ bmat
        z = u @ bmat
        acc += np.array([
            mask.sum(), y.sum(), (y*y).sum(), (y**4).sum(),
            (y*fields).sum(), np.maximum(np.abs(y), np.abs(z)).sum(),
            (u*y).sum(), np.cos(y).sum(),
        ])
        yy = y[:, selected]
        gg = fields[:, selected]
        row_acc += np.stack((yy.sum(axis=0), (yy**2).sum(axis=0),
                             (yy**4).sum(axis=0), (yy*gg).sum(axis=0),
                             np.cos(yy).sum(axis=0)), axis=1)
    acc /= samples*n
    row_acc /= samples
    gram = bmat @ bmat
    return {
        "model": label, "order": n, "samples": samples,
        "threshold": "7/8", "gaussian_variance_target": b,
        "gaussian_fourth_target": 3*b*b,
        "gaussian_cosine_target": math.exp(-b/2),
        "one_probe_a": aval,
        "observed": dict(zip(("threshold_mass", "mean", "second", "fourth",
                              "covariance_with_first_field", "paired_response",
                              "paired_initial_half_energy", "cosine"), acc.tolist())),
        "gram_frobenius_squared_over_n_squared": float((gram*gram).sum()/n**2),
        "selected_rows": [{"row": row, "mean_second_fourth_cov_cos": vals.tolist()}
                          for row, vals in zip(selected, row_acc)],
        "classification": "finite numerical stress test, not asymptotic proof",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-order", type=int, default=1024)
    parser.add_argument("--samples", type=int, default=8192)
    parser.add_argument("--batch", type=int, default=256)
    parser.add_argument("--seed", type=int, default=20260905)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    for label, matrix in models(rng, args.max_order):
        print(json.dumps(test(label, matrix, rng, args.samples, args.batch)), flush=True)


if __name__ == "__main__":
    main()
