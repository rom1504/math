#!/usr/bin/env python3
"""Replay the complex-evaluation identity used by the L2 compiler barrier.

These floating diagnostics supplement the analytic proof; no asymptotic
claim or exact rank is inferred from them.
"""
import json
import math

import numpy as np


def fwht(values):
    values = np.array(values, copy=True)
    width = 1
    while width < len(values):
        for start in range(0, len(values), 2*width):
            left = values[start:start+width].copy()
            right = values[start+width:start+2*width].copy()
            values[start:start+width] = left+right
            values[start+width:start+2*width] = left-right
        width *= 2
    return values


def main():
    rng = np.random.RandomState(17172026)
    cases = 0
    max_ratio = 0.0
    for n in range(1, 9):
        indices = np.arange(2**n)
        bits = ((indices[:, None] >> np.arange(n)) & 1)
        signs = 1-2*bits
        for trial in range(50):
            u = (math.pi/4)*2.0**(-rng.randint(0, 12))
            angles = rng.uniform(u/2, u, size=n)
            D = 1+trial % n
            coeff = rng.normal(size=2**n)+1j*rng.normal(size=2**n)
            coeff[bits.sum(axis=1) > D] = 0
            values = fwht(coeff)
            cap = np.max(np.abs(values))
            coeff /= cap
            values /= cap
            target = np.exp(1j*(signs @ angles))
            correlation = np.mean(target.conj()*values)
            point = -1j*np.tan(angles)
            monomials = np.prod(np.where(bits, point[None, :], 1), axis=1)
            evaluation = np.dot(coeff, monomials)
            exact_identity = np.prod(np.cos(angles))*evaluation
            assert abs(correlation-exact_identity) < 2e-12
            bound = np.prod(np.cos(angles))*math.exp(D*math.asinh(math.tan(u)))
            assert abs(correlation) <= bound+2e-12
            probabilities = np.sin(angles)**2
            entropy = sum(-p*math.log(p)-(1-p)*math.log1p(-p) for p in probabilities)
            entropy_bound = n*u*u*math.log(math.e/(u*u))
            assert entropy <= entropy_bound+2e-12
            assert math.asinh(math.tan(u)) <= math.sqrt(2)*u+1e-12
            assert sum(-np.log(np.cos(angles)))+1e-12 >= n*u*u/8
            max_ratio = max(max_ratio, abs(correlation)/bound)
            cases += 1
    dyadic_sum = sum((math.pi/4)*2.0**(-j)*
                     math.log(math.e/(((math.pi/4)*2.0**(-j))**2)) for j in range(100))
    exact_sum = (math.pi/2)*(1+2*math.log(8/math.pi))
    assert abs(dyadic_sum-exact_sum) < 1e-12
    assert 8*math.sqrt(2)*exact_sum+2*math.log(2) < 64
    print(json.dumps({"status": "PASS", "complex_evaluation_cases": cases,
                      "maximum_correlation_over_bound_diagnostic": max_ratio,
                      "dyadic_sum": exact_sum, "entropy_degree_constant": 64}, indent=2))


if __name__ == "__main__":
    main()
