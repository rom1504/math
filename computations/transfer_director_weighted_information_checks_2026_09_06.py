"""Finite numerical falsification checks, not substitutes for the proof.

Exhausts every selector for each generated matrix. Includes sparse pairs,
uneven weights, and nonidentical retention probabilities. Seed is explicit.
"""
import itertools
import json
import math
from pathlib import Path

import numpy as np


def check_matrix(b, probabilities, ratio=0.35):
    n = len(b)
    op = float(np.linalg.norm(b, 2))
    if op == 0:
        return 0
    t = ratio / op
    c = t * t / (2 * (1 - ratio * ratio))
    w = b * b
    d = float(w.sum(axis=1).max())
    exponent = 2 * c * d
    weighted = sum(w[i, j] / d * math.log1p(probabilities[i] * probabilities[j] * math.expm1(exponent))
                   for i in range(n) for j in range(i + 1, n))
    determinant_mean = exponential_mean = 0.0
    count = 0
    for bits in itertools.product((0, 1), repeat=n):
        select = np.flatnonzero(bits)
        mass = math.prod(p if bit else 1 - p for p, bit in zip(probabilities, bits))
        small = b[np.ix_(select, select)]
        eigen = np.linalg.eigvalsh(small)
        logdet = -0.5 * float(np.log1p(-t * t * eigen * eigen).sum())
        upper = c * float((small * small).sum())
        assert logdet <= upper + 1e-12
        determinant_mean += mass * math.exp(logdet)
        exponential_mean += mass * math.exp(upper)
        count += 1
    assert math.log(determinant_mean) <= math.log(exponential_mean) + 1e-12
    assert math.log(exponential_mean) <= weighted + 1e-12
    return count


def main():
    rng = np.random.default_rng(20260906)
    instances = selectors = 0
    for n in range(2, 9):
        for kind in ('dense', 'sparse', 'pairs'):
            for _ in range(6):
                u = np.triu(rng.normal(size=(n, n)), 1)
                if kind == 'sparse':
                    u *= rng.random((n, n)) < 0.3
                if kind == 'pairs':
                    u[:] = 0
                    for i in range(0, n - 1, 2):
                        u[i, i + 1] = rng.uniform(0.1, 2)
                b = u + u.T
                probabilities = rng.uniform(0.001, 0.8, n)
                selectors += check_matrix(b, probabilities)
                instances += 1
    result = {'status': 'NUMERICAL FINITE CHECKS ONLY', 'seed': 20260906,
              'matrices': instances, 'selector_instances': selectors,
              'scope': 'weighted determinant/resampling/Finner inequalities; no Gaussian entropy quadrature'}
    target = Path(__file__).resolve().parent / 'results/transfer_director_weighted_information_checks_2026_09_06.json'
    target.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
