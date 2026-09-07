"""Exact n=16 cap screen under symmetric closed-quadruple switching."""
import itertools
import json
from pathlib import Path
import numpy as np


def sylvester(n):
    h = np.ones((1, 1), dtype=np.int64)
    while len(h) < n:
        h = np.block([[h, h], [h, -h]])
    return h


def switches(h):
    n = len(h)
    out = []
    for quad in itertools.combinations(range(n), 4):
        q = np.array(quad)
        if not np.all(np.prod(h[q], axis=0) == 1):
            continue
        block = h[np.ix_(q, q)]
        rr = np.ones((4, 4), dtype=np.int64) - 2*np.eye(4, dtype=np.int64)
        transformed = rr @ block @ rr
        if not np.all(np.abs(transformed) == 4):
            continue
        new = h.copy()
        new[q] = (rr @ new[q]) // 2
        new[:, q] = (new[:, q] @ rr) // 2
        assert np.all(np.abs(new) == 1)
        assert np.array_equal(new, new.T)
        assert np.array_equal(new @ new, n*np.eye(n, dtype=np.int64))
        if not np.array_equal(new, h):
            out.append((quad, new))
    return out


def main():
    n = 16
    h = sylvester(n)
    spins = 1 - 2*((np.arange(1 << (n-1), dtype=np.int64)[:, None]
                   >> np.arange(n)) & 1)
    rng = np.random.default_rng(20260907)
    records = []
    best = (n*n, 1 << n)
    for step in range(401):
        energies = np.einsum('bi,ij,bj->b', spins, h, spins, optimize=True)
        energies = (energies - np.trace(h)) // 2
        cap = int(np.max(np.abs(energies)))
        ground = int(np.count_nonzero(np.abs(energies) == cap))
        options = switches(h)
        if (cap, ground) < best:
            best = (cap, ground)
            ground_vectors = spins[np.abs(energies) == cap]
            ground_rank = int(np.linalg.matrix_rank(ground_vectors.astype(float)))
            print(json.dumps(dict(step=step, cap=cap, ground=ground,
                                  switches=len(options))), flush=True)
            records.append(dict(step=step, cap=cap, ground=ground,
                                ground_rank=ground_rank,
                                matrix=h.tolist()))
        if step % 20 == 0:
            print(json.dumps(dict(step=step, cap=cap, ground=ground,
                                  switches=len(options))), flush=True)
        if not options:
            break
        _, h = options[int(rng.integers(len(options)))]
    output = dict(n=n, steps=step, best=best, records=records)
    Path('computations/results/flatify_construct_2026_09_07_symmetric_hadamard_switch.json').write_text(
        json.dumps(output, indent=2) + '\n')


if __name__ == '__main__':
    main()
