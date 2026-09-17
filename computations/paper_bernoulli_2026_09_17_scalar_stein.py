"""Numerically integrate the exact quantile W1 for finite weighted sign sums.

The normal antiderivative is analytic, but floating-point evaluation is a
diagnostic, not a certified bound.  The theorem is proved by Stein's identity.
"""

import itertools
import json
import math

import numpy as np
from scipy.special import ndtr, ndtri


def normal_density_at_quantile(probability):
    if probability <= 0 or probability >= 1:
        return 0.0
    return math.exp(-float(ndtri(probability))**2/2)/math.sqrt(2*math.pi)


def exact_atom_wasserstein(weights):
    weights = np.asarray(weights, dtype=np.int64)
    variance = int(weights @ weights)
    values = np.array([sum(s*w for s,w in zip(signs,weights))
                       for signs in itertools.product((-1,1), repeat=len(weights))])
    atoms, counts = np.unique(values, return_counts=True)
    left = 0.0
    answer = 0.0
    denominator = 2**len(weights)
    for atom, count in zip(atoms,counts):
        right = left + int(count)/denominator
        point = int(atom)/math.sqrt(variance)
        middle = min(right,max(left,float(ndtr(point))))
        pleft = normal_density_at_quantile(left)
        pmiddle = normal_density_at_quantile(middle)
        pright = normal_density_at_quantile(right)
        answer += point*(middle-left)+pmiddle-pleft
        answer += -pright+pmiddle-point*(right-middle)
        left = right
    assert abs(left-1) < 1e-12
    bound = float(np.sum(np.abs(weights)**3))/variance**1.5
    return answer,bound


if __name__ == '__main__':
    rng = np.random.default_rng(2026091728)
    cases = []
    for length in range(1,15):
        candidates = [np.ones(length,dtype=np.int64),
                      np.full(length,3,dtype=np.int64)]
        candidates += [rng.integers(1,8,size=length) for _ in range(5)]
        for weights in candidates:
            measured,bound = exact_atom_wasserstein(weights)
            assert measured <= bound+1e-12
            cases.append(dict(length=length, weights=list(map(int,weights)),
                              wasserstein=measured,stein_bound=bound))
    print(json.dumps(dict(status='PASS', cases=len(cases),
                          maximum_ratio=max(c['wasserstein']/c['stein_bound'] for c in cases),
                          equal_weight_examples=[c for c in cases
                              if all(x==1 for x in c['weights'])],
                          scope='Exact finite atoms; analytic normal integrals evaluated numerically'),indent=2))
