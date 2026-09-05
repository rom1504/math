#!/usr/bin/env python3
"""Counterexample-first diagnostics for two successive rooted H2 layers.

G=B S, Y=B[S H2(G)]/sqrt(2), W=B[S H2(Y)]/sqrt(2).
Outputs raw moments, with whole-spin sample standard errors; these are
finite diagnostics, not a CLT certificate. Root zero is also tracked so
vanishing planted subsets are not averaged away.
"""
from __future__ import annotations

import argparse
import json
import math

import numpy as np
from scipy.integrate import quad

from fresh_custom_response_falsify import models, random_signing, hollow


def extra_models(n, rng):
    for count in (1, 2):
        a = random_signing(n, rng)
        k = max(2, round(math.sqrt(n)))
        for group in range(count):
            lo, hi = group*k, (group+1)*k
            pattern = rng.choice((-1., 1.), size=n-hi)
            a[lo:hi, hi:] = pattern
            a[hi:, lo:hi] = pattern[:, None]
            a[lo:hi, lo:hi] = 1
        yield f"planted_{count}_twin_groups_sqrt_n", hollow(a)


NAMES = ('G', 'Y', 'W', 'G2', 'Y2', 'W2', 'G4', 'Y4', 'W4', 'W6',
         'GY', 'GW', 'YW', 'SW', 'G2W2', 'Y2W2', 'G3W', 'Y3W',
         'G2Y2', 'GW3', 'YW3', 'Wtail1', 'Wtail2', 'Wtail3', 'Wcos1')


def test(label, a, rng, samples, batch, channel, scale):
    n = len(a)
    b = a/math.sqrt(n-1)
    eigs = np.linalg.eigvalsh(b)
    if channel == 'H2':
        h = lambda x: (x*x-1)/math.sqrt(2)
        moments = {'mean': 0., 'sd': math.sqrt(2)}
    else:
        def raw_h(x):
            return scale*np.tanh((x*x-1)/scale)
        density = lambda x: math.exp(-x*x/2)/math.sqrt(2*math.pi)
        mu = quad(lambda x: float(raw_h(x))*density(x), -12, 12, epsabs=1e-12)[0]
        var = quad(lambda x: (float(raw_h(x))-mu)**2*density(x), -12, 12, epsabs=1e-12)[0]
        h = lambda x: (raw_h(x)-mu)/math.sqrt(var)
        moments = {'mean': mu, 'sd': math.sqrt(var)}
    acc = np.zeros((2, len(NAMES)))
    sq = np.zeros_like(acc)
    for start in range(0, samples, batch):
        count = min(batch, samples-start)
        spins = rng.choice((-1., 1.), size=(count, n))
        g = spins@b
        y = (spins*h(g))@b
        w = (spins*h(y))@b
        raw = (g, y, w, g*g, y*y, w*w, g**4, y**4, w**4, w**6,
               g*y, g*w, y*w, spins*w, g*g*w*w, y*y*w*w, g**3*w,
               y**3*w, g*g*y*y, g*w**3, y*w**3,
               (np.abs(w)>1).astype(float), (np.abs(w)>2).astype(float),
               (np.abs(w)>3).astype(float), np.cos(w))
        for category in range(2):
            vals = np.stack([x.mean(1) if category == 0 else x[:, 0]
                             for x in raw], axis=1)
            acc[category] += vals.sum(0)
            sq[category] += (vals*vals).sum(0)
    means = acc/samples
    se = np.sqrt(np.maximum(0, sq/samples-means*means)/(samples-1))
    return {'model': label, 'order': n, 'samples': samples,
            'channel': channel, 'channel_scale': scale, 'normalizing_moments': moments,
            'op_B': float(np.max(np.abs(eigs))),
            'spectral_third': float(np.mean(eigs**3)),
            'spectral_fourth': float(np.mean(eigs**4)),
            'estimates': {category: {key: {'mean': float(value), 'sample_se': float(error)}
                          for key, value, error in zip(NAMES, mu, serr)}
                          for category, mu, serr in zip(('coordinate_average', 'root_zero'), means, se)},
            'classification': 'finite Monte Carlo falsifier, not proof'}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--orders', type=int, nargs='+', default=(128, 256, 512))
    p.add_argument('--samples', type=int, default=8192)
    p.add_argument('--batch', type=int, default=128)
    p.add_argument('--seed', type=int, default=2026090525)
    p.add_argument('--quick', action='store_true')
    p.add_argument('--only', nargs='+', default=None)
    p.add_argument('--channel', choices=('H2', 'tanh'), default='H2')
    p.add_argument('--scale', type=float, default=4.)
    p.add_argument('--output', default='computations/results/fresh_second_rooted_layer_falsify_2026_09_05.json')
    args = p.parse_args()
    rng = np.random.default_rng(args.seed)
    results = []
    for n in args.orders:
        for label, a in list(models(n, rng, args.quick)) + list(extra_models(n, rng)):
            if args.only and not any(pattern in label for pattern in args.only):
                continue
            result = test(label, a, rng, args.samples, args.batch, args.channel, args.scale)
            results.append(result)
            av = result['estimates']['coordinate_average']
            rz = result['estimates']['root_zero']
            print(json.dumps({'order': n, 'model': label,
                              'bulk': {k: av[k]['mean'] for k in ('Y2', 'Y4', 'W2', 'W4', 'GW', 'YW')},
                              'root': {k: rz[k]['mean'] for k in ('W2', 'W4', 'GW', 'YW')}}), flush=True)
    # This program's JSON output is a generated computation artifact.
    with open(args.output, 'w', encoding='utf-8') as handle:
        json.dump(results, handle, indent=2)
        handle.write('\n')


if __name__ == '__main__':
    main()
