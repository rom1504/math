#!/usr/bin/env python3
"""Finite-type adversary for a signed two-clique seed graphon.

The exact Boolean-type maximum is a LOWER bound on the full graphon cap;
continuous type magnetizations are not exhausted. Outer minimization is
heuristic and is never labeled as a minimum certificate.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution


def spin_table(n):
    return (2 * ((np.arange(1 << n)[:, None] >> np.arange(n)) & 1) - 1)


TYPES = spin_table(3)
SPINS = spin_table(8).astype(float)
SIGMA, TAU, SWITCH = TYPES.T
A = (SIGMA[:, None] + SIGMA[None, :]) / 2
B = (TAU[:, None] + TAU[None, :]) / 2 * SWITCH[:, None] * SWITCH[None, :]


def measure(params):
    corr, *biases = params
    group = ((SIGMA > 0).astype(int) + 2 * (TAU > 0).astype(int))
    return (1 + corr * SIGMA * TAU) * (1 + SWITCH * np.asarray(biases)[group]) / 8


def lower_cap(params, return_witness=False):
    mu = measure(params)
    wx = SPINS * mu
    qa = np.einsum("bi,ij,bj->b", wx, A, wx) / 2
    energy = qa[:, None] - qa[None, :] + wx @ B @ wx.T
    index = int(np.argmax(energy))
    cap = float(energy.flat[index])
    if not return_witness:
        return cap
    i, j = np.unravel_index(index, energy.shape)
    return {"params": np.asarray(params).tolist(), "mu": mu.tolist(),
            "boolean_type_lower_cap": cap, "old_graphon_cap": 1/8,
            "ratio_lower": 8 * cap, "x": SPINS[i].astype(int).tolist(),
            "y": SPINS[j].astype(int).tolist(),
            "types_sigma_tau_s": TYPES.tolist()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--maxiter", type=int, default=100)
    parser.add_argument("--popsize", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260918)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    baseline = lower_cap(np.zeros(5), True)
    print("Uniform independent types:", json.dumps(baseline), flush=True)
    result = differential_evolution(lower_cap, [(-1, 1)] * 5,
                                    seed=args.seed, maxiter=args.maxiter,
                                    popsize=args.popsize, tol=1e-9,
                                    polish=False, updating="immediate")
    data = {"seed": args.seed, "maxiter": args.maxiter,
            "popsize": args.popsize, "uniform_types": baseline,
            "best_outer_candidate": lower_cap(result.x, True),
            "outer_evaluations": int(result.nfev),
            "outer_certificate": False,
            "note": "Boolean-type caps are lower bounds; outer minimization is heuristic."}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(data, indent=2) + "\n")
    print(json.dumps(data), flush=True)


if __name__ == "__main__":
    main()
