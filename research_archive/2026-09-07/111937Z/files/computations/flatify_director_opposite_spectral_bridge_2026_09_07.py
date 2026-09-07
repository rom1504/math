"""Finite actual-child test, not an asymptotic certificate.

All sampled bridge caps are exact integer maxima; Gaussian sampling and its
spectral square roots are floating-point exploratory computations. Inputs are
tracked order-8 minimizer representatives. No ignored research input.
"""
import argparse
import json
from pathlib import Path
import numpy as np


def roots(a):
    val, vec = np.linalg.eigh(a)
    return [(vec * np.sqrt(np.maximum(s * val, 0))) @ vec.T for s in (1, -1)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--samples', type=int, default=1000)
    ap.add_argument('--seed', type=int, default=202609071031)
    ap.add_argument('--output', default='computations/results/flatify_director_opposite_spectral_bridge_2026_09_07.json')
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    inp = Path('computations/results/m8_minimizer_orbits.json')
    classes = json.loads(inp.read_text())['classes']
    n = 8
    spins = np.ones((1 << (n-1), n), dtype=np.int64)
    spins[:, 1:] = 1-2*((np.arange(1 << (n-1))[:, None] >> np.arange(n-1)) & 1)
    out = {'status': 'exploratory sampling; each reported cap exact by exhaustive integer spin enumeration',
           'seed': args.seed, 'samples_per_case': args.samples, 'input': str(inp), 'cases': []}
    for ai, ca in enumerate(classes):
        a = np.array(ca['representative_matrix'], dtype=np.int64)
        ha = np.einsum('bi,ij,bj->b', spins, a, spins)//2
        assert max(abs(ha)) == 10
        pa, na = roots(a)
        for di, cd in enumerate(classes):
            for polarity in (1, -1):
                d = polarity*np.array(cd['representative_matrix'], dtype=np.int64)
                hd = np.einsum('bi,ij,bj->b', spins, d, spins)//2
                pd, nd = roots(d)
                internal = np.abs(ha[:, None]+hd[None, :])
                for law in ('opposite_spectral', 'iid'):
                    hist = {}
                    best = 10**9
                    witness = None
                    for rep in range(args.samples):
                        if law == 'opposite_spectral':
                            g = pa @ rng.normal(size=(n,n)) @ nd + na @ rng.normal(size=(n,n)) @ pd
                        else:
                            g = rng.normal(size=(n,n))
                        c = np.where(g >= 0, 1, -1).astype(np.int64)
                        cap = int(np.max(internal + np.abs(spins @ c @ spins.T)))
                        hist[cap] = hist.get(cap, 0)+1
                        if cap < best:
                            best, witness = cap, c.tolist()
                    case = {'child_a_class': ai, 'child_d_class': di, 'polarity': polarity,
                            'law': law, 'best_exact_cap': best, 'cap_histogram': hist,
                            'best_bridge': witness}
                    out['cases'].append(case)
                    print(ai, di, polarity, law, best, flush=True)
                    Path(args.output).write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__':
    main()
