"""Rank-two exact sign tiles: numerical lower bounds, NOT cap certificates."""
import argparse
import json
from pathlib import Path
import numpy as np


def sylvester(k):
    h = np.ones((1, 1), dtype=np.int64)
    while len(h) < k:
        h = np.kron(h, [[1, 1], [1, -1]])
    assert len(h) == k
    return h


def construct(m, rng, randomize=True):
    k = 2*m
    h = sylvester(k)
    src = Path(f'computations/results/exact_m{m}.json')
    s = np.array(json.loads(src.read_text())['matrix'], dtype=np.int64) if src.exists() else sylvester(m)
    np.fill_diagonal(s, 1)
    bases = []
    for i in range(m):
        order = rng.permutation(k) if randomize else np.arange(k)
        phases = rng.choice([-1, 1], k) if randomize else np.ones(k, dtype=int)
        bases.append(h[:, order]*phases)
    c = np.zeros((m*k, m*k), dtype=np.int64)
    h2 = np.array([[1, 1], [1, -1]], dtype=np.int64)
    for i in range(m):
        for j in range(m):
            c[i*k:(i+1)*k, j*k:(j+1)*k] = s[i, j]*(bases[i][:, 2*j:2*j+2] @ h2 @ bases[j][:, 2*i:2*i+2].T)//2
    assert np.array_equal(c, c.T) and np.all(np.abs(c) == 1)
    assert np.array_equal(c @ c, m*k*np.eye(m*k, dtype=np.int64))
    assert np.trace(c) == 0
    return c


def climb(c, rng, restarts=512):
    n = len(c)
    a = c.copy()
    np.fill_diagonal(a, 0)
    best = 0
    bestx = None
    for sigma in [-1, 1]:
        x = rng.choice([-1, 1], (restarts, n))
        field = x @ (sigma*a)
        for _ in range(20*n):
            gains = -2*x*field
            inds = np.argmax(gains, axis=1)
            active = np.flatnonzero(gains[np.arange(restarts), inds] > 0)
            if not len(active):
                break
            flips = inds[active]
            old = x[active, flips].copy()
            x[active, flips] *= -1
            field[active] -= 2*old[:, None]*(sigma*a[flips])
        vals = np.sum(x*field, axis=1)//2
        ind = int(np.argmax(vals))
        if int(vals[ind]) > best:
            best = int(vals[ind]); bestx = x[ind].tolist()
    assert abs(np.array(bestx) @ a @ np.array(bestx))//2 == best
    return best, bestx


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--m', type=int, default=4)
    p.add_argument('--instances', type=int, default=20)
    p.add_argument('--restarts', type=int, default=512)
    args = p.parse_args()
    rng = np.random.default_rng(20260907+args.m)
    records = []
    for index in range(args.instances):
        c = construct(args.m, rng, randomize=index != 0)
        lower, witness = climb(c, rng, args.restarts)
        row = dict(index=index, randomized=index != 0, cap_lower=lower,
                   normalized_lower=lower/len(c)**1.5, matrix=c.tolist(), witness=witness)
        records.append(row)
        print(json.dumps({k: v for k, v in row.items() if k not in ['matrix', 'witness']}), flush=True)
    Path(f'computations/results/flatify_construct_2026_09_07_rank_two_weave_m{args.m}.json').write_text(
        json.dumps(dict(status='heuristic lower bounds only; all sign and Hadamard identities exact', records=records), indent=2)+'\n')


if __name__ == '__main__':
    main()
