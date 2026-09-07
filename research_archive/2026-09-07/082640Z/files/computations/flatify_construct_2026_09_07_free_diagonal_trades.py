"""Exact whole-matrix closed-quad trades with free hollow diagonals."""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np


def options(c):
    n = len(c)
    out = []
    for q in itertools.combinations(range(n), 4):
        outside = [i for i in range(n) if i not in q]
        parity = np.prod(c[np.ix_(q, outside)], axis=0)
        if not np.all(parity == parity[0]):
            continue
        b = c[np.ix_(q, q)].copy()
        b[0] *= parity[0]
        b[:, 0] *= parity[0]
        np.fill_diagonal(b, 1)
        np.fill_diagonal(b, np.prod(b, axis=0))
        if int(b.sum()) % 8 == 0:
            out.append(q)
    return out


def trade(c, q):
    q = list(q)
    h = c + np.eye(len(c), dtype=np.int64)
    outside = [i for i in range(len(c)) if i not in q]
    parity = int(np.prod(c[q, outside[0]]))
    h[q[0]] *= parity
    h[:, q[0]] *= parity
    b = h[np.ix_(q, q)]
    h[q, q] = np.prod(b, axis=0)
    r = np.ones((4, 4), dtype=np.int64) - 2*np.eye(4, dtype=np.int64)
    assert np.all(np.prod(h[q], axis=0) == 1)
    h[q] = (r @ h[q]) // 2
    h[:, q] = (h[:, q] @ r) // 2
    assert np.all(np.abs(h) == 1) and np.array_equal(h, h.T)
    np.fill_diagonal(h, 0)
    return h


def score(c, x):
    e = np.einsum('bi,ij,bj->b', x, c, x, optimize=True)//2
    cap = int(np.max(np.abs(e)))
    return cap, int(np.count_nonzero(np.abs(e) == cap))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--start', choices=['cap30', 'sylvester'], default='cap30')
    p.add_argument('--steps', type=int, default=100)
    p.add_argument('--sample', type=int, default=32)
    args = p.parse_args()
    if args.start == 'cap30':
        raw = json.loads(Path('computations/results/flatify_adversary_2026_09_07_hadamard_bridge_phase_search.json').read_text())
        c = np.array(raw['records'][1]['best']['parent_matrix'], dtype=np.int64)
    else:
        c = np.array([[1]], dtype=np.int64)
        for _ in range(4):
            c = np.kron(c, [[1, 1], [1, -1]])
        np.fill_diagonal(c, 0)
    n = len(c)
    x = 1-2*((np.arange(1 << (n-1), dtype=np.int64)[:, None] >> np.arange(n)) & 1)
    rng = np.random.default_rng(20260907)
    best = score(c, x)
    records = [dict(step=0, score=best, matrix=c.tolist())]
    log = []
    for step in range(args.steps+1):
        opts = options(c)
        row = dict(step=step, score=score(c, x), options=len(opts))
        log.append(row)
        if step % 10 == 0:
            print(json.dumps(row), flush=True)
        if not opts or step == args.steps:
            break
        if len(opts) > args.sample:
            opts = [opts[i] for i in rng.choice(len(opts), args.sample, replace=False)]
        candidates = [(score(d := trade(c, q), x), q, d) for q in opts]
        candidates.sort(key=lambda a: a[0])
        if step == 0:
            print(json.dumps(dict(first_step_scores=[a[0] for a in candidates])), flush=True)
        selected = candidates[0] if rng.random() < .9 else candidates[int(rng.integers(len(candidates)))]
        current, q, c = selected
        if current < best:
            best = current
            records.append(dict(step=step+1, score=best, quad=q, matrix=c.tolist()))
            print(json.dumps(dict(improved=True, step=step+1, score=best)), flush=True)
    Path(f'computations/results/flatify_construct_2026_09_07_free_diagonal_trades_signed_{args.start}.json').write_text(
        json.dumps(dict(start=args.start, best=best, records=records, log=log), indent=2)+'\n')


if __name__ == '__main__':
    main()
