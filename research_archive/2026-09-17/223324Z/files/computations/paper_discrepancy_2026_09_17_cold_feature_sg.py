"""Finite diagnostics for the all-direction cold-feature SG proof.

These are reproducible exact finite-sum calculations in floating point,
not certificates of the asymptotic or all-direction theorem.
"""
import itertools
import json
import math

import numpy as np
from scipy.special import gammaln, logsumexp


def sparse_binomial_sequence(n, v):
    half = n // 2
    j = np.arange(half + 1)
    sums = 2*j-half
    logs = (gammaln(half+1)-gammaln(j+1)-gammaln(half-j+1)
            - half*math.log(2))
    total = (sums[:, None]+sums[None, :])/math.sqrt(n)
    contrast = (sums[:, None]-sums[None, :])/math.sqrt(n)
    beta = 1/v-1
    log_weights = logs[:, None]+logs[None, :]-beta*total**2/2
    log_z = logsumexp(log_weights)
    law = np.exp(log_weights-log_z)
    results = {}
    for name, values in (("feature", total), ("perpendicular", contrast)):
        ratios = [2*(logsumexp(log_weights+amplitude*values)-log_z)/amplitude**2
                  for amplitude in np.geomspace(0.01, 8, 40)]
        results[name] = {"variance": float(np.sum(law*values**2)),
                         "largest_sampled_MGF_proxy": max(ratios)}
    return {"n": n, "v": v, "relative_partition": math.exp(log_z)/math.sqrt(v),
            "directions": results}


def finite_frame(u, v, rng):
    n, r = u.shape
    assert np.max(np.abs(u.T@u-np.eye(r))) < 1e-12
    words = np.asarray(list(itertools.product((-1.0, 1.0), repeat=n)))
    beta = 1/v-1
    feature = words@u
    log_density = -beta*np.sum(feature**2, axis=1)/2
    log_norm = logsumexp(log_density)
    probability = np.exp(log_density-log_norm)
    covariance = words.T@(probability[:, None]*words)
    assert np.max(np.abs(probability@words)) < 1e-12
    assert np.max(np.abs(np.diag(covariance)-1)) < 1e-12
    pmax = float(np.max(np.sum(u*u, axis=1)))
    max_ratio = 0.0
    determinant_checks = 0
    for _ in range(30):
        direction = rng.normal(size=n)
        direction /= np.linalg.norm(direction)
        values = words@direction
        for amplitude in np.geomspace(0.01, 8, 20):
            t = amplitude*direction
            d = 1-np.tanh(t)**2
            sigma = u.T@(d[:, None]*u)
            loss = (r*math.log(1+beta)-np.linalg.slogdet(np.eye(r)+beta*sigma)[1])/2
            bound = beta*pmax*float(t@t)/2
            assert loss <= bound+1e-11
            determinant_checks += 1
            ratio = 2*(logsumexp(log_density+amplitude*values)-log_norm)/amplitude**2
            max_ratio = max(max_ratio, float(ratio))
    return {"n": n, "rank": r, "v": v, "max_leverage": pmax,
            "covariance_operator_norm": float(np.linalg.eigvalsh(covariance)[-1]),
            "largest_sampled_MGF_proxy": max_ratio,
            "determinant_inequalities_checked": determinant_checks}


def quartic_extension_checks(v):
    beta = 1/v-1
    cutoff = min(0.5, 1/(2*math.sqrt(beta)))
    grid = np.linspace(-8, 8, 16001)
    magnitude = np.abs(grid)
    gap = np.maximum(0, magnitude-cutoff)
    extension = np.where(magnitude <= cutoff, grid**4,
                         cutoff**4+4*cutoff**3*gap+6*cutoff**2*gap**2)
    second_derivative = 12*np.minimum(grid*grid, cutoff*cutoff)
    assert np.min(extension) >= 0
    assert np.max(extension-grid**4) < 1e-10
    assert np.max(second_derivative) <= 12*cutoff**2+1e-12
    assert 2*beta*cutoff**2 <= 0.5+1e-12
    return {"v": v, "cutoff": cutoff,
            "interpolation_Hessian_upper": 2*beta*cutoff**2,
            "grid_points_checked": len(grid)}


rng = np.random.RandomState(170926)
frames = []
for n, r in ((6, 1), (8, 2), (10, 3)):
    u, _ = np.linalg.qr(rng.normal(size=(n, r)))
    frames.append(finite_frame(u, 0.25, rng))

# Direct scalar checks for the two elementary proof inequalities.
z = np.linspace(-1, 1, 10001)
assert np.all(np.sin(z)**2 >= z*z-z**4/3-1e-15)
for scale in np.linspace(0, 1, 101):
    assert np.all(np.cosh(scale*z)-1 <= scale**2*(np.cosh(z)-1)+1e-14)

sequences = [sparse_binomial_sequence(n, v)
             for v in (0.25, 0.75) for n in (16, 32, 64, 128, 256, 512)]
print(json.dumps({"status": "PASS finite floating replay; no numerical proof claim",
                  "frames": frames, "rank_one_sequences": sequences,
                  "convex_quartic_extension":
                      [quartic_extension_checks(v) for v in (0.1, 0.25, 0.75)]},
                 indent=2))
