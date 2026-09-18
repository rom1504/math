"""Numerical diagnostics for the symmetric Boolean BH route.

Uses only the repository environment. These are floating-point searches,
not global certificates. The original degree constraint is represented in
the Chebyshev basis on the exact Hamming-weight grid. No cube enumeration
or interpolation between grid points is imposed.
"""

from __future__ import annotations

import argparse
import json
import math

import numpy as np
from numpy.polynomial.chebyshev import chebvander
from scipy.optimize import linprog
from scipy.special import gammaln
from scipy.stats import binom


def symmetric_data(n: int, m: int):
    grid = np.arange(n + 1)
    values = chebvander(2 * grid / n - 1, m)
    weights = binom.pmf(grid, n, 0.5)
    edge_weights = binom.pmf(np.arange(n), n - 1, 0.5)
    delta = np.diff(values, axis=0) / 2
    gram = n * delta.T @ (edge_weights[:, None] * delta)
    kraw = np.empty((n + 1, m + 1))
    kraw[:, 0] = 1
    if m:
        kraw[:, 1] = (n - 2 * grid) / math.sqrt(n)
    for r in range(1, m):
        kraw[:, r + 1] = (
            (n - 2 * grid) * kraw[:, r]
            - math.sqrt(r * (n - r + 1)) * kraw[:, r - 1]
        ) / math.sqrt((r + 1) * (n - r))
    transform = kraw.T @ (weights[:, None] * values)
    r = np.arange(m + 1)
    logcounts = gammaln(n + 1) - gammaln(r + 1) - gammaln(n - r + 1)
    return values, gram, transform, logcounts, delta


def diagnostic(n: int, m: int, starts: int, rounds: int, seed: int):
    values, gram, transform, logcounts, delta = symmetric_data(n, m)
    constraints = np.vstack([values, -values])
    rhs = np.ones(2 * (n + 1))
    bounds = [(None, None)] * (m + 1)
    rng = np.random.default_rng(seed + 101 * n + m)

    def optimize_linear(objective):
        sol = linprog(-objective, A_ub=constraints, b_ub=rhs, bounds=bounds,
                      method="highs")
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.x

    central = optimize_linear(delta[n // 2])
    center_ratio = abs(delta[n // 2] @ central) * n / m
    best_influence = 0.0
    best_bh = 0.0
    best_coefficients = None
    q = 2 * m / (m + 1)
    weights_q = np.exp(logcounts / (m + 1))
    max_violation = 0.0
    for run in range(starts):
        if run == 0:
            coeff = central.copy()
        elif run <= m:
            coeff = np.zeros(m + 1)
            coeff[run] = 1
        else:
            coeff = optimize_linear(rng.normal(size=m + 1))
        for _ in range(rounds):
            objective = gram @ coeff
            if not np.linalg.norm(objective):
                break
            nxt = optimize_linear(objective)
            if np.max(np.abs(nxt - coeff)) < 1e-9:
                coeff = nxt
                break
            coeff = nxt
        vals = values @ coeff
        max_violation = max(max_violation, float(np.max(np.abs(vals)) - 1))
        influence = float(coeff @ gram @ coeff)
        spectrum = transform @ coeff
        bh = float(np.sum(np.abs(spectrum) ** q * weights_q) ** (1 / q))
        best_bh = max(best_bh, bh)
        if influence > best_influence:
            best_influence = influence
            best_coefficients = coeff.tolist()
    return {
        "n": n, "m": m,
        "central_derivative_n_over_m": float(center_ratio),
        "max_found_total_influence": best_influence,
        "max_found_influence_n_over_m_squared": best_influence * n / m**2,
        "max_found_bh_ratio": best_bh,
        "max_constraint_violation": max_violation,
        "best_influence_chebyshev_coefficients": best_coefficients,
        "status": "numerical; local convex-maximization search, not a proof",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", default="4:2,8:2,16:2,16:4,32:4,64:4,64:8,128:8,256:8,256:16,512:16,1024:16")
    parser.add_argument("--starts", type=int, default=12)
    parser.add_argument("--rounds", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260918)
    args = parser.parse_args()
    for pair in args.pairs.split(","):
        n, m = map(int, pair.split(":"))
        print(json.dumps(diagnostic(n, m, args.starts, args.rounds, args.seed)), flush=True)


if __name__ == "__main__":
    main()
