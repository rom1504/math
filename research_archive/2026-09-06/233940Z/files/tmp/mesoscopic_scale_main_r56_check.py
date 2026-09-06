#!/usr/bin/env python3
"""Numerical checks for the Wave 56 entropy-matched core calculation."""

import math

import numpy as np
from scipy.special import gammaln
from scipy.stats import hypergeom


def logcomb(n: int, k):
    k = np.asarray(k)
    ans = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    return np.where((k >= 0) & (k <= n), ans, -np.inf)


p = 0.6
gamma = 0.75
kappa = (1 - p) ** 2 / (2 * p**2)
threshold_constant = p / (math.sqrt(2 * math.pi) * (1 - p))
print(f"p={p}, gamma={gamma}, predicted kappa={kappa:.9f}")
print(f"predicted threshold constant={threshold_constant:.9f}")
print(" n    ell   shift/ell  var/n    cost-ratio  max*ell/sqrt(n)  max/lambda2")

for n in (1000, 2000, 4000, 8000, 16000):
    m = int(round(p * n))
    ell = int(round(n**gamma))
    n_tail = n - ell
    d = m - ell

    mu0 = m * m / n
    mu = ell + d * d / n_tail
    shift = mu - mu0
    exact_shift = ell * (n - m) ** 2 / (n * (n - ell))
    assert abs(shift - exact_shift) < 1e-9

    exact_var = d * d * (n - m) ** 2 / (n_tail**2 * (n_tail - 1))
    scipy_var = hypergeom.var(n_tail, d, d)
    assert abs(exact_var - scipy_var) < 1e-8 * max(1, exact_var)

    j_lo = max(ell, 2 * m - n)
    js = np.arange(j_lo, m, dtype=int)  # omit the self-loop j=m
    log_a = logcomb(m, js) + logcomb(n - m, m - js)
    log_d = float(logcomb(m, ell) + logcomb(n - ell, m - ell))
    log_k = logcomb(js, ell) - log_d
    log_ball = np.logaddexp.accumulate(log_a[::-1])[::-1]
    log_cap = log_ball + log_k
    best = float(np.max(log_cap))
    max_load = math.exp(best)
    scaled = max_load * ell / math.sqrt(n)
    lambda2 = (
        ell
        * (ell - 1)
        * (n - m)
        * (n - m - 1)
        / (m * (m - 1) * (n - ell) * (n - ell - 1))
    )

    j_mu = int(round(mu))
    log_a_mu = float(logcomb(m, j_mu) + logcomb(n - m, m - j_mu))
    log_nslice = float(logcomb(n, m))
    moderate_cost = -(log_a_mu - log_nslice + 0.5 * math.log(n))
    cost_ratio = moderate_cost / (ell * ell / n)

    print(
        f"{n:5d} {ell:6d} {shift/ell:10.6f} {exact_var/n:8.5f} "
        f"{cost_ratio:11.7f} {scaled:17.8f} {max_load/lambda2:12.5f}"
    )

    assert abs(shift / ell - (1 - p) ** 2) < 0.08
    assert 0.02 < exact_var / n < 0.08
    assert 0.05 < cost_ratio < 0.5
    assert 0.1 < scaled < 20

print("PASS exact parameters, moderate-deviation scale, and threshold order")
