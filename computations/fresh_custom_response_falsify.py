#!/usr/bin/env python3
"""Counterexample-first Monte Carlo for the proposed custom response energy.

Tests e_y=E(yp^T B yp-ym^T B ym)/(4n), with
yp=sign(aS+Y), ym=sign(aS-Y), Y=B[S 1{|BS|<=7/8}].
Sample standard errors are computed over independent whole spin vectors,
not by treating correlated coordinates as independent observations.
No finite Monte Carlo output certifies an asymptotic signing family.
"""

from __future__ import annotations

import argparse
import json
import math

import numpy as np


def hollow(a):
    a = np.asarray(a, dtype=np.float64)
    np.fill_diagonal(a, 0)
    assert np.array_equal(a, a.T)
    assert np.all(np.abs(a + np.eye(len(a))) == 1)
    return a


def random_signing(n, rng):
    a = np.triu(rng.choice((-1., 1.), size=(n, n)), 1)
    return a + a.T


def hadamard(power):
    h = np.ones((1, 1))
    h4 = np.ones((4, 4)) - 2*np.eye(4)
    for _ in range(power):
        h = np.kron(h, h4)
    return h


def models(n, rng, quick=False):
    yield "random", random_signing(n, rng)
    ratios = (0.1, 0.25, 1., 4.) if quick else (0.03, 0.06, 0.1, 0.25, 0.5, 1., 2., 4.)
    for ratio in ratios:
        rank = max(2, round(ratio*n))
        x = rng.standard_normal((n, rank))
        gram = x @ x.T / math.sqrt(rank)
        yield f"sign_gram_rank_ratio_{ratio:g}", hollow(np.where(gram >= 0, 1., -1.))
        if not quick and ratio in (0.1, 0.5, 1.):
            noise = rng.standard_normal((n, n))
            noise = (noise + noise.T)/math.sqrt(2)
            for strength in (0.5, 2.):
                yield f"sign_gram_ratio_{ratio:g}_noise_{strength:g}", hollow(
                    np.where(gram + strength*noise >= 0, 1., -1.))
    a = random_signing(n, rng)
    k = int(n**0.75)
    a[:k, :k] = 1
    yield "planted_positive_clique_n_three_quarters", hollow(a)
    a = random_signing(n, rng)
    a[:k, :k] = -1
    yield "planted_negative_clique_n_three_quarters", hollow(a)
    for power in range(1, 6):
        if 4**power == n:
            yield "regular_hadamard", hollow(hadamard(power).copy())
        for seed_order in (2, 3, 5, 6, 8):
            if seed_order*4**power == n:
                seed = random_signing(seed_order, rng)
                np.fill_diagonal(seed, rng.choice((-1., 1.), seed_order))
                yield f"tensor_random_full_seed_{seed_order}", hollow(np.kron(hadamard(power), seed))


def test(label, a, rng, samples, batch):
    n = len(a)
    bmat = a/math.sqrt(n-1)
    eigs = np.linalg.eigvalsh(bmat)
    t = 7/8
    aval = math.sqrt(2/math.pi)*math.exp(-t*t/2)
    acc = np.zeros(5)
    acc2 = np.zeros(5)
    for start in range(0, samples, batch):
        count = min(batch, samples-start)
        s = rng.choice((-1., 1.), size=(count, n))
        g = s @ bmat
        mask = np.abs(g) <= t
        v = s*mask
        u = np.sign(g)*(~mask)
        y = v @ bmat
        bu = u @ bmat
        yp = np.where(aval*s+y >= 0, 1., -1.)
        ym = np.where(aval*s-y >= 0, 1., -1.)
        byp, bym = yp @ bmat, ym @ bmat
        vals = np.stack((
            ((yp*byp).sum(1)-(ym*bym).sum(1))/(4*n),
            (yp*(bu+y)+ym*(bu-y)).sum(1)/(2*n),
            (u*y).sum(1)/n,
            np.maximum(aval, np.abs(y)).sum(1)/n,
            (yp*byp+ym*bym).sum(1)/(4*n),
        ), axis=1)
        acc += vals.sum(0)
        acc2 += (vals*vals).sum(0)
    means = acc/samples
    se = np.sqrt(np.maximum(0, acc2/samples-means*means)/(samples-1))
    names = ("oriented_custom_energy", "custom_cross", "initial_energy",
             "rooted_max_field", "unoriented_custom_energy")
    return {
        "model": label, "order": n, "samples": samples,
        "threshold": "7/8", "a": aval,
        "spectral": {"op": float(np.max(np.abs(eigs))),
                     "third_moment": float(np.mean(eigs**3)),
                     "fourth_moment": float(np.mean(eigs**4)),
                     "positive_fraction": float(np.mean(eigs>0))},
        "estimates": {key: {"mean": float(value), "sample_se": float(error)}
                      for key, value, error in zip(names, means, se)},
        "classification": "finite Monte Carlo falsification diagnostic only",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orders", type=int, nargs="+", default=(64, 128, 256))
    parser.add_argument("--samples", type=int, default=8192)
    parser.add_argument("--batch", type=int, default=256)
    parser.add_argument("--seed", type=int, default=2026090520)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    for n in args.orders:
        for label, a in models(n, rng, args.quick):
            print(json.dumps(test(label, a, rng, args.samples, args.batch)), flush=True)


if __name__ == "__main__":
    main()
