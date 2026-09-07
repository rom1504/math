#!/usr/bin/env python3
"""Noncertified finite optimization of actual child-sector subgaussian proxies."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.optimize import differential_evolution, minimize


ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / "extremal_information" / "experiments"
sys.path.insert(0, str(EXP))

import actual_child_bridge_law_exact as exact  # noqa: E402
import actual_child_escort_low_degree_falsifier as low_degree  # noqa: E402


BETA = 4.0
PLANS = ((6, 3, 3), (8, 4, 4), (10, 3, 7))
RANDOM_STARTS = 48
DE_MAXITER = 500
DE_POPSIZE = 20


def sector_law(matrix, total_order, sector):
    spins = exact.projective_spins(len(matrix)).astype(np.float64)
    energy = exact.energies_for_matrix(matrix, spins.astype(np.int16)).astype(
        np.float64
    )
    t = BETA / math.sqrt(total_order)
    log_weights = sector * t * energy
    log_weights -= float(np.max(log_weights))
    weights = np.exp(log_weights)
    weights /= float(np.sum(weights))
    # The full sector law is the even lift of the projective representatives.
    covariance = spins.T @ (weights[:, None] * spins)
    return spins, weights, covariance


def optimize_proxy(matrix, total_order, sector, seed):
    spins, weights, covariance = sector_law(matrix, total_order, sector)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    lambda_max = float(eigenvalues[-1])
    top_vector = eigenvectors[:, -1]
    dimension = len(matrix)

    # Since log E exp(<u,X>) <= ||u||_1 <= sqrt(k)||u|| and the tangent
    # lower bound is lambda_max, no optimizer lies outside this radius.
    radius_bound = 2.0 * math.sqrt(dimension) / lambda_max

    def proxy(u):
        u = np.asarray(u, dtype=np.float64)
        norm_square = float(np.dot(u, u))
        if norm_square < 1e-18:
            return lambda_max
        scores = spins @ u
        # Even lift: E exp(<u,X>)=E_projective cosh(<u,X>).
        # The log1p form preserves the quadratic behavior near zero.
        excess = float(np.dot(weights, 2.0 * np.sinh(scores / 2.0) ** 2))
        log_mgf = math.log1p(excess)
        return 2.0 * log_mgf / norm_square

    def objective(u):
        norm = float(np.linalg.norm(u))
        if norm > radius_bound:
            return 100.0 + (norm - radius_bound) ** 2
        return -proxy(u)

    bounds = [(-radius_bound, radius_bound)] * dimension
    de = differential_evolution(
        objective,
        bounds,
        seed=seed,
        maxiter=DE_MAXITER,
        popsize=DE_POPSIZE,
        tol=1e-11,
        atol=1e-12,
        polish=True,
        updating="immediate",
        workers=1,
    )

    rng = np.random.default_rng(seed + 1000)
    starts = [de.x]
    for fraction in (1e-4, 0.01, 0.1, 0.35, 0.7):
        starts.append(fraction * radius_bound * top_vector)
    for _ in range(RANDOM_STARTS):
        direction = rng.normal(size=dimension)
        direction /= np.linalg.norm(direction)
        # Log-uniform radii retain the covariance-tangent regime while also
        # probing the nonlinear boundary regime.
        fraction = math.exp(rng.uniform(math.log(1e-4), math.log(0.999)))
        starts.append(fraction * radius_bound * direction)

    local_rows = []
    for start in starts:
        local = minimize(
            objective,
            start,
            method="Nelder-Mead",
            options={
                "maxiter": 3000,
                "xatol": 1e-10,
                "fatol": 1e-12,
            },
        )
        candidate = np.asarray(local.x, dtype=np.float64)
        if np.linalg.norm(candidate) <= radius_bound * (1 + 1e-8):
            local_rows.append(
                {
                    "proxy": proxy(candidate),
                    "norm": float(np.linalg.norm(candidate)),
                    "u": candidate.tolist(),
                    "optimizer_success": bool(local.success),
                    "optimizer_message": str(local.message),
                }
            )
    local_rows.sort(key=lambda row: -row["proxy"])
    best = local_rows[0]
    estimate = max(lambda_max, float(-de.fun), best["proxy"])
    return {
        "child_order": dimension,
        "parent_order": total_order,
        "sector": sector,
        "covariance_tangent_lower_bound": lambda_max,
        "rigorous_radial_search_bound": radius_bound,
        "differential_evolution_proxy": float(-de.fun),
        "differential_evolution_success": bool(de.success),
        "differential_evolution_message": str(de.message),
        "best_multistart_candidate": best,
        "numerical_proxy_estimate": estimate,
        "top_distinct_local_proxy_values": sorted(
            {round(row["proxy"], 11) for row in local_rows}, reverse=True
        )[:10],
        "local_start_count": len(starts),
    }


def main():
    mp.mp.dps = 80
    child_cache = {}
    sector_cache = {}
    child_records = {}
    for total_order, left_order, right_order in PLANS:
        for child_order in sorted({left_order, right_order}):
            key = (total_order, child_order)
            matrix, certificate = low_degree.child_record(
                child_order, total_order, child_cache
            )
            child_records[f"N{total_order}_k{child_order}"] = certificate
            for sector in (-1, 1):
                print(
                    f"optimizing N={total_order} k={child_order} s={sector:+d}",
                    flush=True,
                )
                sector_cache[(total_order, child_order, sector)] = optimize_proxy(
                    matrix,
                    total_order,
                    sector,
                    seed=10000 * total_order + 100 * child_order + sector,
                )
                row = sector_cache[(total_order, child_order, sector)]
                print(
                    f"  tangent={row['covariance_tangent_lower_bound']:.9g} "
                    f"estimate={row['numerical_proxy_estimate']:.9g} "
                    f"norm={row['best_multistart_candidate']['norm']:.6g}",
                    flush=True,
                )

    products = []
    for total_order, left_order, right_order in PLANS:
        for epsilon in (-1, 1):
            sector_products = {}
            for sector in (-1, 1):
                left = sector_cache[(total_order, left_order, sector)][
                    "numerical_proxy_estimate"
                ]
                right = sector_cache[(total_order, right_order, epsilon * sector)][
                    "numerical_proxy_estimate"
                ]
                sector_products[str(sector)] = left * right
            products.append(
                {
                    "N": total_order,
                    "split": [left_order, right_order],
                    "orientation": epsilon,
                    "sector_products": sector_products,
                    "kappa_star_numerical_estimate": max(
                        sector_products.values()
                    ),
                }
            )

    result = {
        "schema": "scratch-actual-child-sector-subgaussian-wind-tunnel-v1",
        "classification": (
            "certified finite child input; noncertified numerical global and "
            "multistart optimization of exact finite sector MGFs"
        ),
        "definition": (
            "sigma^2(nu)=sup_{u!=0} 2 log E_nu exp(<u,X>)/||u||^2; "
            "kappa_*=max_s sigma^2(mu_A,s)sigma^2(mu_D,epsilon*s)"
        ),
        "parameters": {
            "beta": BETA,
            "plans": [list(plan) for plan in PLANS],
            "differential_evolution_maxiter": DE_MAXITER,
            "differential_evolution_popsize": DE_POPSIZE,
            "random_multistarts": RANDOM_STARTS,
        },
        "children": child_records,
        "sector_records": [
            sector_cache[key] for key in sorted(sector_cache)
        ],
        "orientation_products": products,
        "scope": (
            "The covariance tangent values are rigorous finite lower bounds. "
            "The optimized proxy and kappa values are reproducible numerical "
            "estimates, not certified upper bounds or asymptotic claims."
        ),
    }
    (ROOT / "tmp/actual_child_sector_subgaussian_wind_tunnel.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
