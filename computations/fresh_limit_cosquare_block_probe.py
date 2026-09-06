"""Exact caps for a deterministic random cosquare block-flip probe.

The search is not exhaustive. Every reported witness is verified by integer
matrix arithmetic and exhaustive Boolean optimization.
"""
import argparse
import itertools
import json
import numpy as np


def beta(a, spins):
    return int(np.max(np.sum(np.abs(spins @ a), axis=1)))


def run(k, m, trials, seed, stop_first=True):
    rng = np.random.default_rng(seed)
    p = -np.ones((k, k), dtype=np.int16)
    # Symmetric circular row-balanced seed, including a positive diagonal.
    radius = (k // 2 - 1) // 2
    for i in range(k):
        for d in range(-radius, radius + 1):
            p[i, (i+d) % k] = 1
        if k % 4 == 0:
            p[i, (i+k//2) % k] = 1
    assert np.all(p.sum(axis=0) == 0)
    spins = np.array(list(itertools.product((-1, 1), repeat=k+m)), dtype=np.int16)
    witnesses = []
    for trial in range(trials):
        for _ in range(100):
            a,b,c,d = rng.choice(k, 4, replace=False)
            if p[a,b] == p[c,d] and p[a,c] == p[b,d] == -p[a,b]:
                for u,v in ((a,b),(c,d),(a,c),(b,d)):
                    p[u,v] = p[v,u] = -p[u,v]
        d = rng.choice(np.array([-1,1], dtype=np.int16), size=(m,m))
        d = np.triu(d) + np.triu(d,1).T
        j = np.ones((k,m), dtype=np.int16)
        c = np.block([[p,j],[j.T,d]])
        cp = np.block([[-p,j],[j.T,d]])
        assert np.array_equal(c @ c, cp @ cp)
        b0,b1 = beta(c,spins), beta(cp,spins)
        if b0 != b1:
            item = dict(trial=trial, beta_original=b0, beta_flipped=b1,
                        original=c.tolist(), flipped=cp.tolist())
            if not witnesses or abs(b0-b1) > abs(witnesses[0]['beta_original']-witnesses[0]['beta_flipped']):
                witnesses = [item]
            if stop_first:
                break
    print(json.dumps(dict(k=k,m=m,trials=trials,seed=seed,witnesses=witnesses)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--k', type=int, default=8)
    parser.add_argument('--m', type=int, default=4)
    parser.add_argument('--trials', type=int, default=500)
    parser.add_argument('--seed', type=int, default=20260905)
    parser.add_argument('--scan-all', action='store_true')
    args = parser.parse_args()
    run(args.k,args.m,args.trials,args.seed,not args.scan_all)
