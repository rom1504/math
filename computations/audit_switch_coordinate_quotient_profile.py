#!/usr/bin/env python3
"""Profile labelled switch quotients for one exact rare bridge.

This is a finite diagnostic, not asymptotic evidence.  It takes a bridge
gauge mask (normally the best full-switch bridge from
``audit_switch_fourier_actual_children.py``), keeps each possible subset of
the ``2*n-1`` independent vertex-switch coordinates, and computes the exact
conditional-expectation convolution certificate on that quotient.  With
``--all-subspaces`` it exhausts every binary linear character subspace; that
mode is capped at switch dimension seven.
"""

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

from audit_switch_fourier_actual_children import (
    energies,
    enumerate_child_minimizer,
    fwht_rows,
    logmeanexp,
    spins,
)


def rref_subspaces(dimension, ell):
    """Yield each ell-dimensional binary subspace once, by an RREF basis."""
    if ell == 0:
        yield (), ()
        return
    for pivots in itertools.combinations(range(dimension), ell):
        free_positions = []
        pivot_set = set(pivots)
        for row, pivot in enumerate(pivots):
            for column in range(pivot + 1, dimension):
                if column not in pivot_set:
                    free_positions.append((row, column))
        for free_mask in range(1 << len(free_positions)):
            basis = [1 << pivot for pivot in pivots]
            for bit, (row, column) in enumerate(free_positions):
                if (free_mask >> bit) & 1:
                    basis[row] |= 1 << column
            yield pivots, tuple(basis)


def span_indices(basis):
    indices = []
    for coefficient_mask in range(1 << len(basis)):
        value = 0
        for j, vector in enumerate(basis):
            if (coefficient_mask >> j) & 1:
                value ^= vector
        indices.append(value)
    return indices


def quotient_profile(n, beta, epsilon, bridge_mask, all_subspaces=False):
    X = spins(n)
    A, child_pressure = enumerate_child_minimizer(n, beta)
    t = beta / math.sqrt(2 * n)

    y_tail = spins(n - 1)
    Y = np.ones((1 << (n - 1), n), dtype=np.int8)
    Y[:, 1:] = y_tail
    Xg = np.repeat(X, len(Y), axis=0)
    Yg = np.tile(Y, (len(X), 1))
    dimension = 2 * n - 1
    group_size = 1 << dimension
    assert len(Xg) == group_size

    h_x = energies(A, Xg)
    h_y = energies(A, Yg)
    w = np.cosh(t * (h_x + epsilon * h_y))
    mean_w = float(np.mean(w))
    a_hat = fwht_rows((w / mean_w)[None, :], normalized=True)[0]

    variable_edges = [(i, j) for i in range(1, n) for j in range(1, n)]
    features = np.stack(
        [Xg[:, i] * Yg[:, j] for i, j in variable_edges], axis=1
    ).astype(np.float64)
    signs = np.array(
        [1.0 - 2.0 * ((bridge_mask >> bit) & 1)
         for bit in range(len(variable_edges))],
        dtype=np.float64,
    )
    bridge_energy = (
        Xg.sum(axis=1).astype(np.float64) * Yg.sum(axis=1)
        + features.dot(signs - 1.0)
    )
    scaled = t * bridge_energy
    psi = float(logmeanexp(scaled[None, :], axis=1)[0])
    shifted = scaled - np.max(scaled)
    b = np.exp(shifted)
    b /= np.mean(b)
    b_hat = fwht_rows(b[None, :], normalized=True)[0]
    products = a_hat * b_hat
    full_convolution = fwht_rows(products[None, :], normalized=False)[0]
    direct_base = math.log(mean_w) + psi - 2.0 * child_pressure

    records = []
    if all_subspaces and dimension > 7:
        raise ValueError("exhaustive linear-subspace mode is capped at dimension 7")

    for ell in range(dimension + 1):
        best = None
        if all_subspaces:
            candidates = rref_subspaces(dimension, ell)
        else:
            candidates = (
                (coordinates, tuple(1 << coordinate for coordinate in coordinates))
                for coordinates in itertools.combinations(range(dimension), ell)
            )
        candidate_count = 0
        for pivots, basis in candidates:
            candidate_count += 1
            mode_indices = span_indices(basis)
            quotient_fourier = products[np.asarray(mode_indices, dtype=np.int64)]
            quotient_values = fwht_rows(
                quotient_fourier[None, :], normalized=False
            )[0]
            quotient_min = float(np.min(quotient_values))
            candidate = {
                "rref_pivots": list(pivots),
                "basis_masks": list(basis),
                "minimum_quotient_convolution": quotient_min,
                "log_gain": math.log(max(quotient_min, 1e-300)),
                "certificate": direct_base
                + math.log(max(quotient_min, 1e-300)),
            }
            if best is None or candidate["certificate"] < best["certificate"]:
                best = candidate
        records.append(
            {
                "quotient_dimension": ell,
                "subspaces_checked": candidate_count,
                **best,
            }
        )

    return {
        "classification": (
            "exact finite quotient enumeration with floating-point "
            "transcendental evaluation"
        ),
        "child_order": n,
        "parent_order": 2 * n,
        "beta": beta,
        "epsilon": epsilon,
        "bridge_gauge_mask": bridge_mask,
        "switch_coordinate_dimension": dimension,
        "quotient_family": (
            "all binary linear subspaces" if all_subspaces
            else "coordinate-generated subspaces"
        ),
        "base_before_switch_gain": direct_base,
        "full_convolution_minimum": float(np.min(full_convolution)),
        "full_certificate": direct_base
        + math.log(max(float(np.min(full_convolution)), 1e-300)),
        "profiles": records,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--beta", type=float, required=True)
    parser.add_argument("--epsilon", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--bridge-mask", type=int, required=True)
    parser.add_argument("--all-subspaces", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.n > 5:
        raise SystemExit("The exact child audit is intentionally capped at n=5")
    result = quotient_profile(
        args.n, args.beta, args.epsilon, args.bridge_mask, args.all_subspaces
    )
    payload = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n")
    print(payload)


if __name__ == "__main__":
    main()
