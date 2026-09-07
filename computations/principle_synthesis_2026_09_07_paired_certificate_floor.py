"""Outward interval check for a .497 floor of the paired Bellman certificate.

This bounds a specific two-witness certificate, NOT the actual anti-family.
The reduction to 700 fixed rational p intervals is proved in the artifact.
"""
from fractions import Fraction as F
import json

import mpmath as mp
from flatify_independent_2026_09_07_ternary_interval import endpoints, rational_iv


def entropy(p):
    if p in (F(0), F(1)):
        return mp.iv.mpf(0)
    x = rational_iv(p)
    return -x * mp.iv.ln(x) - (1-x) * mp.iv.ln(1-x)


def gaussian(t):
    lam = 2*t/(4*t+1+mp.iv.sqrt(16*t*t+1))
    return mp.iv.ln(lam*(2*t-lam)/(t*t))/4-lam


def verify(precision=50):
    mp.iv.dps = precision
    b, p0 = F(497, 1000), F(93, 100)
    bb, pp0 = rational_iv(b), rational_iv(p0)
    log2 = mp.iv.ln(2)
    first_endpoint = pp0*log2+mp.iv.ln(1-4*bb*bb*pp0)/4
    assert endpoints(first_endpoint)[0] > 0
    assert endpoints(entropy(p0))[1] < F(3, 10)
    assert endpoints(2*bb*mp.iv.sqrt(pp0))[0] > F(9, 10)
    # Thus h(p) < r/(1+r), uniformly for p>=p0.
    worst, worst_interval = None, None
    for j in range(9300, 10000):
        left, right = F(j, 10000), F(j+1, 10000)
        rright = 2*bb*mp.iv.sqrt(rational_iv(right))
        time_upper = entropy(left)/(2*(1-rright))
        lower = rational_iv(left)*log2+gaussian(time_upper)+entropy(right)/2
        lower_endpoint = endpoints(lower)[0]
        assert lower_endpoint > 0, (left, right, endpoints(lower))
        if worst is None or lower_endpoint < worst:
            worst, worst_interval = lower_endpoint, (left, right)
    return {"status": "DIRECTED_INTERVAL_CERTIFICATE",
            "scope": "PAIRED_BELLMAN_CERTIFICATE_FLOOR_NOT_ACTUAL_CAP",
            "floor": str(b), "precision": precision, "intervals": 700,
            "first_endpoint_lower": str(endpoints(first_endpoint)[0]),
            "boundary_uniform_lower": str(worst),
            "boundary_uniform_lower_display": float(worst),
            "worst_interval": list(map(str, worst_interval))}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
