#!/usr/bin/env python3
"""Finite audit of the Wave 56 active-face minimax theorem.

The proof itself is analytic.  This checker exhausts a small signing problem,
forms the exact response game, and solves both its mixed-law LP and its pure
fractional-block side.  It also checks the conditioning inequality directly.
"""

from itertools import combinations, product
from math import log, sqrt

import numpy as np
from scipy.optimize import linprog


def edge_list(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def state_patterns(n, edges):
    """Distinct x_i x_j edge patterns; orientation is added later."""
    seen = set()
    rows = []
    for tail in product((-1, 1), repeat=n - 1):
        x = (1,) + tail
        row = tuple(x[i] * x[j] for i, j in edges)
        if row not in seen:
            seen.add(row)
            rows.append(row)
    return np.asarray(rows, dtype=int)


def exact_minimizer(n):
    edges = edge_list(n)
    patterns = state_patterns(n, edges)
    best_q = None
    best_a = None
    # Gauge-fix one edge sign.  Global negation has the same absolute cap.
    for tail in product((-1, 1), repeat=len(edges) - 1):
        a = np.asarray((1,) + tail, dtype=int)
        energies = 2 * (patterns @ a)
        q = int(np.max(np.abs(energies)))
        if best_q is None or q < best_q:
            best_q, best_a = q, a.copy()
    return edges, patterns, best_a, best_q


def all_blocks(edge_count, r):
    ans = [tuple()]
    for size in range(1, r + 1):
        ans.extend(combinations(range(edge_count), size))
    return ans


def main():
    n, r = 5, 3
    edges, patterns, a, q = exact_minimizer(n)

    energies0 = 2 * (patterns @ a)
    energies = np.concatenate((energies0, -energies0))
    tau = np.concatenate((np.ones(len(patterns), dtype=int),
                          -np.ones(len(patterns), dtype=int)))
    repeated_patterns = np.vstack((patterns, patterns))
    signed_edges = tau[:, None] * repeated_patterns * a[None, :]
    deficits = q - energies
    assert np.min(deficits) == 0 and np.all(deficits >= 0)

    blocks = all_blocks(len(edges), r)
    payoff = np.empty((len(blocks), len(deficits)))
    for j, block in enumerate(blocks):
        block_sum = (signed_edges[:, block].sum(axis=1)
                     if block else np.zeros(len(deficits)))
        payoff[j] = deficits + 4 * block_sum

    # min_mu max_block E_mu payoff.  Extreme blocks are exactly the vertices
    # needed for integer r in the uniform-matroid polytope.
    state_count = len(deficits)
    objective = np.r_[np.zeros(state_count), 1.0]
    aub = np.c_[payoff, -np.ones(len(blocks))]
    bub = np.zeros(len(blocks))
    aeq = np.r_[np.ones(state_count), 0.0][None, :]
    result = linprog(
        objective,
        A_ub=aub,
        b_ub=bub,
        A_eq=aeq,
        b_eq=np.array([1.0]),
        bounds=[(0.0, None)] * state_count + [(None, None)],
        method="highs",
    )
    assert result.success, result.message
    mu = result.x[:-1]
    game_value = result.x[-1]

    mean_delta = float(mu @ deficits)
    mean_edges = mu @ signed_edges
    hr = sum(sorted((max(0.0, z) for z in mean_edges), reverse=True)[:r])
    assert abs(game_value - (mean_delta + 4 * hr)) < 1e-7

    v = n
    eta = sqrt(32 * r * (v + 1) * log(2)) + (8 / 3) * (v + 1) * log(2)
    assert game_value <= eta + 1e-8

    # Solve max_p min_omega {Delta+4 p.s} directly.  A concave lower envelope
    # can attain its maximum in the interior, so checking only vertices would
    # not be sufficient here.
    edge_count = len(edges)
    pure_objective = np.r_[np.zeros(edge_count), -1.0]
    pure_aub = np.c_[np.vstack((-4 * signed_edges,
                               np.ones((1, edge_count)))),
                       np.r_[np.ones(state_count), 0.0]]
    pure_bub = np.r_[deficits, r]
    pure_result = linprog(
        pure_objective,
        A_ub=pure_aub,
        b_ub=pure_bub,
        bounds=[(0.0, 1.0)] * edge_count + [(None, None)],
        method="highs",
    )
    assert pure_result.success, pure_result.message
    pure_value = -pure_result.fun
    assert abs(pure_value - game_value) < 1e-7
    assert pure_value <= eta + 1e-8

    # Test additional interior fractional perturbations directly.
    rng = np.random.default_rng(5603)
    for _ in range(500):
        p = rng.random(len(edges))
        p *= min(1.0, r / p.sum())
        minimum_response = np.min(deficits + 4 * (signed_edges @ p))
        assert minimum_response <= eta + 1e-8

    # Directly audit the conditioning algebra.  Markov may make the finite
    # example's good set all states; the asserted inequality is still exact.
    epsilon = 0.2
    good = deficits <= eta / epsilon
    good_mass = float(mu[good].sum())
    assert good_mass >= 1 - epsilon - 1e-9
    conditional_edges = (mu[good] @ signed_edges[good]) / good_mass
    conditional_hr = sum(
        sorted((max(0.0, z) for z in conditional_edges), reverse=True)[:r]
    )
    conditioning_bound = (eta / 4 + epsilon * r) / (1 - epsilon)
    assert conditional_hr <= conditioning_bound + 1e-8

    print(f"n={n}, edges={len(edges)}, exact minimum cap q={q}")
    print(f"states={state_count}, blocks={len(blocks)}")
    print(f"mixed/pure game value={game_value:.9f}, eta_r={eta:.9f}")
    print(f"E Delta={mean_delta:.9f}, h_r(E s)={hr:.9f}")
    print("active-face minimax finite audit: PASS")


if __name__ == "__main__":
    main()
