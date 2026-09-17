"""Finite replay for the exact quartic isotropic laws (2026-09-17).

No optimizer geometry is inferred from these identities.  Outputs are kept
inside the campaign's workspace-owned temporary directory.
"""
import itertools
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tmp/paper_portfolio_2026_09_17/localization/quartic_isotropic_laws.json"
RNG = np.random.default_rng(20260917)


def cube(n):
    z = np.arange(1 << n, dtype=np.uint64)
    return (1 - 2 * ((z[:, None] >> np.arange(n, dtype=np.uint64)) & 1).astype(np.int8)).astype(float)


def hadamard(n):
    h = np.ones((1, 1))
    while len(h) < n:
        h = np.block([[h, h], [h, -h]])
    assert len(h) == n
    return h


def close(a, b, tol=2e-9):
    assert np.max(np.abs(np.asarray(a) - np.asarray(b))) <= tol, (a, b)


def features(h, u):
    n, r = u.shape
    pi = np.sum(u * u, axis=1)
    s = h @ u
    rr = np.sum(s * s, axis=1)
    vv = h @ (pi[:, None] * u)
    f = rr * rr - 2 * (r + 2) * rr + r * (r + 2) + 8 * np.sum(s * vv, axis=1) - 6 * np.sum(pi * pi)
    return f


def check_density(h, w):
    assert np.min(w) >= -2e-9
    close(w.mean(), 1)
    close(w @ h / len(h), np.zeros(h.shape[1]))
    close(h.T @ (w[:, None] * h) / len(h), np.eye(h.shape[1]))


feature_rows = []
query_checks = 0
for n in (8, 16):
    h = cube(n)
    for r in range(1, 5):
        u = hadamard(n)[:, :r] / math.sqrt(n)
        f = features(h, u)
        lam = 1 / (2 * (r + 2))
        w = 1 + lam * f
        check_density(h, w)
        close(f.mean(), 0)
        assert np.mean(f * f) <= 8 * r * (r + 2) + 1e-8
        assert np.mean(w * w) < 3
        assert np.min(f) >= -2 * (r + 2) - 1e-8
        m = n - 1
        mu = np.abs(h.sum(axis=1)).mean()
        indices = sorted(set([0, len(h) - 1] + RNG.integers(0, len(h), size=24).tolist()))
        for j in indices:
            response = np.mean(w * np.abs(h @ h[j]))
            predicted = mu * (1 - lam * f[j] / (m * (m - 2)))
            close(response, predicted)
            query_checks += 1
        feature_rows.append(dict(n=n, rank=r, lambda_=lam, min_density=float(w.min()), density_l2_squared=float(np.mean(w*w))))

# General nonconstant leverage, including a nearly localized direction.
for n, r in ((7, 2), (8, 3), (9, 1), (10, 4)):
    h = cube(n)
    u, _ = np.linalg.qr(RNG.normal(size=(n, r)))
    d = float(np.max(np.linalg.norm(u, axis=1)))
    er = 8*r*math.sqrt(r+2)*d + 3*(2*r*d)**(4/3) + 6*r*d*d
    lam = 1/(2*(r+2)+er)
    f = features(h, u)
    w = 1 + lam*f
    check_density(h, w)
    assert np.min(f) >= -2*(r+2)-er-1e-8
    assert np.mean(f*f) <= 8*r*(r+2)+1e-8
    m = n-1 if n % 2 == 0 else n
    mu = np.abs(h.sum(axis=1)).mean()
    for j in RNG.integers(0, len(h), size=20):
        close(np.mean(w*np.abs(h@h[j])), mu*(1-lam*f[j]/(m*(m-2))))
        query_checks += 1
    feature_rows.append(dict(n=n, rank=r, lambda_=lam, min_density=float(w.min()), density_l2_squared=float(np.mean(w*w))))

energy_rows = []
for n in (6, 8, 10, 12):
    h = cube(n)
    a = RNG.choice((-1.0, 1.0), size=(n, n))
    a = np.triu(a, 1)
    a += a.T
    hh = .5*np.sum((h@a)*h, axis=1)
    d_a = n*(n-1)/2
    p4 = hh*hh - np.sum((h@a)**2, axis=1) + d_a
    ll = max(1.0, float(np.max(np.abs(np.linalg.eigvalsh(a))))/math.sqrt(n))
    lam = 1/(2*ll*ll)
    w = 1 + lam*p4/(n*n)
    check_density(h, w)
    assert w.min() >= .5-1e-9
    pp = float(np.mean(p4*p4))
    assert 4*math.comb(n,4)-1e-8 <= pp <= 36*math.comb(n,4)+1e-8
    close(np.mean(w*hh*hh)/(n*n), d_a/(n*n)+lam*pp/n**4)
    mu = np.abs(h.sum(axis=1)).mean()
    m = n-1
    for j in RNG.integers(0,len(h),size=24):
        close(np.mean(w*np.abs(h@h[j])), mu*(1-lam*p4[j]/(n*n*m*(m-2))))
        query_checks += 1
    energy_rows.append(dict(n=n, op_normalized=ll, min_density=float(w.min()), energy_variance_gain=lam*pp/n**4))

# Inspect the limiting response in a Boolean feature direction at larger n.
limit_rows = []
kappa = math.sqrt(2/math.pi)
for r in (1, 2, 3, 4):
    for n in (64, 256, 1024, 4096):
        lam = 1/(2*(r+2))
        f = n*n-(2*(r+2)-8*r/n)*n+r*(r+2)-6*r*r/n
        mu = math.exp(math.log(n)+math.lgamma(n+1)-2*math.lgamma(n/2+1)-n*math.log(2))
        response = mu/math.sqrt(n)*(1-lam*f/((n-1)*(n-3)))
        limit_rows.append(dict(n=n, rank=r, response=response, limiting=kappa*(1-lam)))

result = dict(status="PASS", feature_cases=feature_rows, energy_cases=energy_rows,
              exact_query_checks=query_checks, limiting_responses=limit_rows)
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(result, indent=2)+"\n")
print(json.dumps(dict(status="PASS", feature_cases=len(feature_rows), energy_cases=len(energy_rows), exact_query_checks=query_checks)))
