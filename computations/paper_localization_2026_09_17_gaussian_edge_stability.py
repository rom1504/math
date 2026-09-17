"""Finite algebra/certificate replay; not a Monte Carlo proof of concentration."""

from itertools import product
from pathlib import Path
import json

import numpy as np


def words(n):
    return np.array(list(product((-1, 1), repeat=n)), dtype=np.int64)


def main():
    rng = np.random.default_rng(20260917)
    counts = {"edge_metric": 0, "tensor_comparison": 0,
              "shell_budget": 0, "continuum_certificate": 0,
              "tested_perturbations": 0}

    for n in range(2, 17):
        ii, jj = np.triu_indices(n, 1)
        for _ in range(300):
            x, y = rng.choice((-1, 1), size=(2, n))
            r = np.count_nonzero(x != y)
            vx, vy = x[ii] * x[jj], y[ii] * y[jj]
            inc = np.sum((vx - vy) ** 2)
            assert inc == 4 * r * (n-r)
            assert inc <= n * np.sum((x-y) ** 2)
            counts["edge_metric"] += 1

            anchor = rng.choice((-1, 1), size=n)
            va = anchor[ii] * anchor[jj]
            dx = rng.choice((-1, 1)) * vx - va
            dy = rng.choice((-1, 1)) * vy - va
            u, v = rng.normal(size=(2, 4))
            u, v = u / np.linalg.norm(u), v / np.linalg.norm(v)
            lhs = np.sum((np.outer(u, dx)-np.outer(v, dy)) ** 2)
            rhs = 2*np.sum((dx-dy)**2)+8*len(ii)*np.sum((u-v)**2)
            assert lhs <= rhs + 1e-8
            counts["tensor_comparison"] += 1

    for _ in range(1200):
        a, big_a, big_b = 10 ** rng.uniform(-4, 4, size=3)
        threshold = 16 * (a*big_a + a*a*big_b*big_b)
        for j in range(12):
            mean_bound = a*(big_a+big_b*np.sqrt(2**(j+1)*threshold))
            assert mean_bound <= 2**(j-1)*threshold*(1+1e-12)
            counts["shell_budget"] += 1

    for n in range(3, 10):
        x = words(n)
        ii, jj = np.triu_indices(n, 1)
        queries = x[:, ii] * x[:, jj]
        signed_queries = np.concatenate((queries, -queries), axis=0)
        for _ in range(8):
            signing = rng.choice((-1, 1), size=len(ii))
            energies = signed_queries @ signing
            cap = np.max(energies)
            ground = int(np.argmax(energies))
            deficits = cap - energies
            noises = rng.normal(size=(3, len(ii)))
            increments = (signed_queries-signed_queries[ground]) @ noises.T
            norms = np.linalg.norm(increments, axis=1)
            # Pick a nontrivial nearcode and the exact sufficient radius
            # from the anchored finite certificate, without invoking a
            # probabilistic width estimate for these arbitrary signings.
            threshold = float(np.quantile(deficits, 0.35))
            outside = deficits > threshold
            if not np.any(outside):
                continue
            radius = 0.9*np.min(deficits[outside]/norms[outside])
            assert np.all(radius*norms[outside] < deficits[outside])
            inside = ~outside
            for _ in range(100):
                direction = rng.normal(size=3)
                direction /= np.linalg.norm(direction)
                theta = direction*radius*rng.uniform(0, 1)
                shifted = energies+signed_queries @ (theta @ noises)
                assert abs(np.max(shifted)-np.max(shifted[inside])) < 1e-8
                # Exact anchor-residual representation of the cap.
                anchor_noise = signed_queries[ground] @ (theta @ noises)
                residual = np.max(increments @ theta-deficits)
                assert abs(np.max(shifted)-cap-anchor_noise-residual) < 1e-8
                counts["tested_perturbations"] += 1
            counts["continuum_certificate"] += 1

    result = {"status": "PASS", "checks": counts,
              "scope": "Finite algebra and sufficient witness certificates; no concentration inferred from sampling."}
    output = Path("tmp/paper_portfolio_2026_09_17/localization/gaussian_edge_stability_audit.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
