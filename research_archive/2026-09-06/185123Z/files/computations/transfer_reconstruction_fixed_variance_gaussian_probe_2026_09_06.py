"""Numerical reconnaissance only: fixed-variance Gaussian versus exact Ising."""

import argparse
from itertools import product
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp, xlogy


def entropy(means):
    p = (1+means)/2
    return -np.sum(xlogy(p, p)+xlogy(1-p, 1-p))


def product_optimum(matrix, starts, rng):
    n = len(matrix)
    def objective(m):
        value = entropy(m)+m@matrix@m/2
        gradient = -np.arctanh(m)+matrix@m
        return -value, -gradient
    best = (n*math.log(2), np.zeros(n))
    for initial in starts:
        result = minimize(objective, initial, method="L-BFGS-B", jac=True,
                          bounds=[(-1+1e-13, 1-1e-13)]*n,
                          options={"ftol": 1e-13, "gtol": 1e-9, "maxiter": 400})
        if -result.fun > best[0]:
            best = (-result.fun, result.x)
    return best


def gaussian_fixed_variance(matrix):
    n = len(matrix)
    diagonal = np.full(n, 1+np.linalg.norm(matrix, 2))
    for _ in range(150):
        precision = np.diag(diagonal)-matrix
        covariance = np.linalg.inv(precision)
        gradient = (1-np.diag(covariance))/2
        value = (np.sum(diagonal)-n-np.linalg.slogdet(precision)[1])/2
        if np.max(np.abs(gradient)) < 1e-11:
            break
        hessian = covariance*covariance/2
        step = np.linalg.solve(hessian, -gradient)
        rate = 1.0
        for _ in range(60):
            candidate = diagonal+rate*step
            trial = np.diag(candidate)-matrix
            if np.linalg.eigvalsh(trial)[0] > 0:
                trial_value = (np.sum(candidate)-n-np.linalg.slogdet(trial)[1])/2
                if trial_value <= value+0.01*rate*gradient@step:
                    diagonal = candidate
                    break
            rate /= 2
        else:
            break
    return value, float(np.max(np.abs(gradient)))


def run(samples, seed):
    rng = np.random.default_rng(seed)
    records = []
    for n in range(3, 9):
        spins = np.array(list(product((-1., 1.), repeat=n)))
        edges = [(i, j) for i in range(n) for j in range(i+1, n)]
        patterns = [np.ones(len(edges)), -np.ones(len(edges))]
        patterns += [rng.choice((-1., 1.), len(edges)) for _ in range(samples)]
        for pattern_id, pattern in enumerate(patterns):
            a = np.zeros((n, n))
            for z, (i, j) in zip(pattern, edges):
                a[i, j] = a[j, i] = z
            bare_energies = np.einsum('bi,ij,bj->b', spins, a, spins)/2
            high_spins = spins[np.argsort(bare_energies)[-8:]]
            starts = [np.zeros(n)] + list(.97*high_spins)
            starts += list(rng.uniform(-.95, .95, (8, n)))
            for beta in (.2, .5, 1., 1.5, math.sqrt(3), 2., 4., 8.):
                interaction = beta*a/math.sqrt(n)
                product_value, means = product_optimum(interaction, starts, rng)
                variance = 1-means*means
                weighted = np.sqrt(variance)[:, None]*interaction*np.sqrt(variance)[None, :]
                spherical, residual = gaussian_fixed_variance(weighted)
                exact = logsumexp(beta*bare_energies/math.sqrt(n))
                difference = float(exact-product_value-spherical)
                weighted_norm = float(np.linalg.norm(weighted, 2))
                records.append({"n": n, "pattern_id": pattern_id,
                                "signs": [int(z) for z in pattern], "beta": beta,
                                "difference": difference, "weighted_norm": weighted_norm,
                                "max_variance": float(np.max(variance)),
                                "means": means.tolist(), "product_value": float(product_value),
                                "spherical_value": spherical, "logZ": float(exact),
                                "dual_residual": residual,
                                "product_globality": "numerical multistart only"})
    ordered = sorted(records, key=lambda r: r["difference"])
    return {"qualification": "numerical reconnaissance; no universal claim or global optimum certificate",
            "samples_per_order": samples, "seed": seed, "cases": len(records),
            "worst_overall": ordered[:12],
            "worst_weighted_norm_below_one": [r for r in ordered if r["weighted_norm"] < 1-1e-6][:8],
            "worst_small_variance": [r for r in ordered if r["max_variance"] < .1][:8]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--samples', type=int, default=3)
    parser.add_argument('--seed', type=int, default=20260906)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run(args.samples, args.seed)
    serialized = json.dumps(result, indent=2)
    if args.output:
        args.output.write_text(serialized+'\n')
    print(serialized)


if __name__ == '__main__':
    main()
