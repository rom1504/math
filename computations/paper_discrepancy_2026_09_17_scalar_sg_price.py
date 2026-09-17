"""Exact physical-slice checks and labeled scalar-MGF diagnostics.

The artifact contains the analytic uniform MGF and sharp asymptotic proof;
the finite grids and numerical maximizations below are not certificates.
"""
from fractions import Fraction as F
import itertools
import json
import math

import numpy as np
from scipy.optimize import minimize_scalar


def log_mgf(s, p):
    return float(np.logaddexp(math.log1p(-p),
                             math.log(p) + np.logaddexp(s, -s) - math.log(2)))


def proxy_upper(p):
    ell = math.log(2 / p)
    return min(1.0, max(2 * math.exp(-ell / 2),
                        1 / (2 * ell) + 8 * math.log(2) / ell**2)) / p


n = 12
words = np.asarray(list(itertools.product((-1, 1), repeat=n)), dtype=np.int64)
sums = words.sum(axis=1)
slice_checks = []
for magnitude in (4, 6, 8, 10, 12):
    p = F(n, magnitude**2)
    weights = [((1-p) / math.comb(n, n//2) if total == 0 else
                p / (2 * math.comb(n, (n+magnitude)//2))
                if abs(total) == magnitude else F(0)) for total in sums]
    denominator = math.lcm(*(weight.denominator for weight in weights))
    numerators = np.asarray([int(weight * denominator) for weight in weights],
                            dtype=object)
    assert sum(numerators) == denominator
    assert all(value == 0 for value in words.astype(object).T @ numerators)
    covariance = words.astype(object).T @ (numerators[:, None] * words)
    assert np.array_equal(covariance, denominator*np.eye(n, dtype=object))
    response = F(int(np.abs(sums).astype(object) @ numerators), denominator)
    assert response == F(n, magnitude)
    slice_checks.append({"magnetization": magnitude, "activation": str(p),
                         "mean_absolute_sum": str(response),
                         "normalized_response": float(response)/math.sqrt(n)})

diagnostics = []
for ell in (2, 4, 8, 16, 32, 64, 128):
    p = 2 * math.exp(-ell)
    # Optimize in s/ell near the asymptotic maximizer 2. This is diagnostic.
    optimum = minimize_scalar(lambda u: -2*log_mgf(u*ell, p)/(u*ell)**2,
                              bounds=(0.25, 4), method="bounded")
    unscaled_proxy = -optimum.fun
    lower = 1/(2*ell)
    upper = proxy_upper(p)*p
    assert lower <= unscaled_proxy*(1+1e-9) <= upper*(1+1e-9)
    for ratio in np.linspace(0.01, 8, 2000):
        s = ratio*ell
        assert 2*log_mgf(s,p)/s**2 <= upper*(1+1e-8)
    diagnostics.append({"log_2_over_p": ell,
                        "analytic_unscaled_lower": lower,
                        "analytic_unscaled_upper": upper,
                        "numerical_unscaled_proxy": unscaled_proxy,
                        "numerical_maximizer_s_over_ell": float(optimum.x),
                        "numerical_K_epsilon2_log_inverse_epsilon":
                            unscaled_proxy*math.log(1/p)/2})

print(json.dumps({"status": "PASS exact n12 covariance and response checks",
                  "physical_slice_checks": slice_checks,
                  "scalar_MGF_diagnostics_not_certificates": diagnostics}, indent=2))
