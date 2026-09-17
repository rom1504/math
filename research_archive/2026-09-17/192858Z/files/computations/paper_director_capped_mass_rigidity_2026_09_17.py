#!/usr/bin/env python3
"""Deterministic diagnostics for two proved inequalities; not a proof by sampling.

No ignored input. Output is explicitly chosen by --output and archived at
checkpoints. Uses NumPy/SciPy from the documented project environment.
"""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
from scipy.stats import binom


def logcosh(x):
    return np.logaddexp(x, -x) - math.log(2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    rng = np.random.default_rng(2026091711)
    worst_mgf = -math.inf
    max_roundoff = 0.0
    for case in range(3000):
        d = int(rng.integers(1, 65))
        coef = rng.normal(size=d)
        if case % 3 == 0:
            coef *= rng.random(d) ** 3
        cap = float(np.max(np.abs(coef))) * (1 + rng.random())
        variance = float(coef @ coef)
        fourth = float(np.sum(coef**4))
        count = int(math.floor(fourth / cap**4))
        remainder = max(0.0, fourth - count * cap**4)
        sea = variance - count * cap**2 - math.sqrt(remainder)
        max_roundoff = max(max_roundoff, -sea)
        assert sea > -1e-10
        t = float(rng.uniform(-3, 3))
        actual = float(np.sum(logcosh(t * coef)))
        packed = (t*t*sea/2 + count*float(logcosh(t*cap))
                  + float(logcosh(t*remainder**0.25)))
        chord = t*t*variance/2 - fourth/cap**4 * (
            t*t*cap**2/2 - float(logcosh(t*cap)))
        worst_mgf = max(worst_mgf, actual-packed, packed-chord)
        assert actual <= packed + 2e-9
        assert packed <= chord + 2e-9

    biased_cases = 0
    for n in range(2, 13):
        words = np.array(list(itertools.product((-1, 1), repeat=n)))
        for rho in (0.05, 0.2, 0.37, 0.5):
            negatives = np.sum(words < 0, axis=1)
            probs = rho**negatives * (1-rho)**(n-negatives)
            p = binom.pmf(np.arange(n//2+1), n//2, rho)
            lower = 2*float(np.sum(p[:, None]*p[None, :]*
                                    np.abs(np.arange(len(p))[:, None]
                                           - np.arange(len(p))[None, :])))
            for nplus in range(n+1):
                coeff = np.r_[np.ones(nplus), -np.ones(n-nplus)]
                mean = float(probs @ np.abs(words @ coeff))
                assert mean >= lower - 1e-10
                if n % 2 == 0 and nplus == n//2:
                    assert abs(mean-lower) < 1e-10
                biased_cases += 1

    variance_cases = 0
    finite_response_cases = 0
    for n in range(2, 9):
        words = np.array(list(itertools.product((-1, 1), repeat=n)))
        for rep in range(12):
            # General nonisotropic law, checking the actual covariance factor.
            support = words[rng.choice(len(words), size=min(20, len(words)),
                                       replace=False)]
            weights = rng.dirichlet(np.ones(len(support)))
            response = np.abs(words @ support.T) @ weights
            covariance = (support.T * weights) @ support
            covariance_norm = float(np.linalg.eigvalsh(covariance)[-1])
            upper = np.triu(rng.choice((-1, 1), size=(n, n)), 1)
            matrix = upper + upper.T
            energy = np.einsum('bi,ij,bj->b', words, matrix, words)/2
            cap = float(np.max(np.abs(energy)))
            ground = words[int(np.argmax(np.abs(energy)))]
            polarity = 1 if energy[np.argmax(np.abs(energy))] >= 0 else -1
            for rho in (0.08, 0.25, 0.4):
                distance = np.sum(words != ground, axis=1)
                probs = rho**distance * (1-rho)**(n-distance)
                mean = float(probs @ response)
                variance = float(probs @ (response-mean)**2)
                assert variance <= 4*covariance_norm + 1e-10
                deficit = cap-polarity*energy
                expected = 4*rho*(1-rho)*cap
                assert abs(float(probs @ deficit)-expected) < 1e-9
                variance_cases += 1
                tol = 1.4*expected + 0.01
                probability_floor = 1-expected/tol
                p = binom.pmf(np.arange(n//2+1), n//2, rho)
                bn = 2*float(np.sum(p[:, None]*p[None, :]*
                                      np.abs(np.arange(len(p))[:, None]
                                             - np.arange(len(p))[None, :])))
                bound = bn-2*math.sqrt(covariance_norm/probability_floor)
                actual = float(np.max(response[cap-np.abs(energy) <= tol]))
                assert actual >= bound - 1e-10
                finite_response_cases += 1

    result = dict(status='PASS: numerical diagnostics, analytic proofs separate',
                  seed=2026091711, capped_mgf_cases=3000,
                  largest_mgf_violation=worst_mgf,
                  gaussian_sea_negative_roundoff=max_roundoff,
                  exact_biased_pattern_cases=biased_cases,
                  convex_variance_cases=variance_cases,
                  finite_response_cases=finite_response_cases)
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
