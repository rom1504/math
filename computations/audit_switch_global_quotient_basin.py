#!/usr/bin/env python3
"""Exhaust all effective switch quotients for order-four actual children.

For a fixed quotient dimension, this program searches every bridge modulo
row/column switching, both relative child orientations, and every character
subspace contained in the effective even--even Fourier sector.  It validates
the winning quotient by decoding its affine fibre and averaging the true
fixed-child parent pressures there.

The output is a finite floating-point diagnostic, not an asymptotic claim or
an interval certificate for the transcendental comparisons.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from audit_switch_coordinate_quotient_profile import (
    rref_subspaces,
    span_indices,
)
from audit_switch_fourier_actual_children import (
    energies,
    enumerate_child_minimizer,
    fwht_rows,
    logmeanexp,
    spins,
)


def parity(value: int) -> int:
    return bin(int(value)).count("1") & 1


def lies_in_even_even(mask: int, n: int) -> bool:
    """Test membership in W in the script's Walsh-coordinate ordering.

    The low n-1 bits encode the free right switches and the next n bits the
    left switches.  An allowed character is in the even--even sector exactly
    when the parity of its left subset is even.
    """

    return parity(mask >> (n - 1)) == 0


def effective_subspaces(dimension: int, rank: int, n: int):
    """Yield every rank-dimensional character subspace U contained in W."""

    for pivots, basis in rref_subspaces(dimension, rank):
        if all(lies_in_even_even(vector, n) for vector in basis):
            yield pivots, basis, np.asarray(span_indices(basis), dtype=np.int64)


def quotient_sector(group_index: int, basis: tuple[int, ...]) -> int:
    sector = 0
    for j, vector in enumerate(basis):
        sector |= parity(group_index & vector) << j
    return sector


def search(n: int, beta: float, rank: int) -> dict[str, object]:
    if n != 4:
        raise ValueError("the global audit is intentionally restricted to n=4")

    x = spins(n)
    child, child_pressure = enumerate_child_minimizer(n, beta)
    t = beta / math.sqrt(2 * n)

    y_tail = spins(n - 1)
    y = np.ones((1 << (n - 1), n), dtype=np.int8)
    y[:, 1:] = y_tail
    xg = np.repeat(x, len(y), axis=0)
    yg = np.tile(y, (len(x), 1))
    dimension = 2 * n - 1
    group_size = 1 << dimension
    if len(xg) != group_size:
        raise AssertionError("switch group ordering mismatch")

    h_x = energies(child, xg)
    h_y = energies(child, yg)
    variable_edges = [(i, j) for i in range(1, n) for j in range(1, n)]
    features = np.stack(
        [xg[:, i] * yg[:, j] for i, j in variable_edges], axis=1
    ).astype(np.float64)
    base_bridge_energy = (
        xg.sum(axis=1).astype(np.float64)
        * yg.sum(axis=1).astype(np.float64)
    )

    bridge_count = 1 << len(variable_edges)
    bridge_masks = np.arange(bridge_count, dtype=np.uint64)[:, None]
    bit_positions = np.arange(len(variable_edges), dtype=np.uint64)
    bridge_bits = ((bridge_masks >> bit_positions) & 1).astype(np.float64)
    bridge_energy = base_bridge_energy[None, :] - 2.0 * bridge_bits.dot(features.T)
    scaled_bridge = t * bridge_energy
    psi = logmeanexp(scaled_bridge, axis=1)
    shifted = scaled_bridge - np.max(scaled_bridge, axis=1, keepdims=True)
    bridge_density = np.exp(shifted)
    bridge_density /= np.mean(bridge_density, axis=1, keepdims=True)
    bridge_hat = fwht_rows(bridge_density, normalized=True)

    subspaces = list(effective_subspaces(dimension, rank, n))
    best: dict[str, object] | None = None
    for epsilon in (1, -1):
        w = np.cosh(t * (h_x + epsilon * h_y))
        mean_w = float(np.mean(w))
        child_hat = fwht_rows((w / mean_w)[None, :], normalized=True)[0]
        products = bridge_hat * child_hat[None, :]
        base = math.log(mean_w) + psi - 2.0 * child_pressure

        for pivots, basis, modes in subspaces:
            quotient_values = fwht_rows(products[:, modes], normalized=False)
            minima = np.maximum(np.min(quotient_values, axis=1), 1e-300)
            certificates = base + np.log(minima)
            bridge_index = int(np.argmin(certificates))
            candidate_value = float(certificates[bridge_index])
            if best is not None and candidate_value >= float(best["certificate"]):
                continue
            sector = int(np.argmin(quotient_values[bridge_index]))
            best = {
                "certificate": candidate_value,
                "epsilon": epsilon,
                "bridge_gauge_mask": bridge_index,
                "rref_pivots": list(pivots),
                "basis_masks": list(basis),
                "base_before_switch_gain": float(base[bridge_index]),
                "minimizing_sector": sector,
                "quotient_density_values": quotient_values[bridge_index].tolist(),
                "quotient_certificates": (
                    base[bridge_index]
                    + np.log(np.maximum(quotient_values[bridge_index], 1e-300))
                ).tolist(),
            }

    if best is None:
        raise AssertionError("no quotient candidate was enumerated")

    # Rebuild the winner and validate (1.8) directly on its affine fibre.
    epsilon = int(best["epsilon"])
    bridge_index = int(best["bridge_gauge_mask"])
    basis = tuple(map(int, best["basis_masks"]))
    sector = int(best["minimizing_sector"])
    w = np.cosh(t * (h_x + epsilon * h_y))
    mean_w = float(np.mean(w))
    child_hat = fwht_rows((w / mean_w)[None, :], normalized=True)[0]
    full_convolution = fwht_rows(
        (child_hat * bridge_hat[bridge_index])[None, :], normalized=False
    )[0]
    base = float(best["base_before_switch_gain"])
    x_values = base + np.log(np.maximum(full_convolution, 1e-300))
    fibre = np.asarray(
        [
            group_index
            for group_index in range(group_size)
            if quotient_sector(group_index, basis) == sector
        ],
        dtype=np.int64,
    )
    direct_certificate = float(logmeanexp(x_values[fibre][None, :], axis=1)[0])
    if abs(direct_certificate - float(best["certificate"])) > 2e-9:
        raise AssertionError((direct_certificate, best["certificate"]))

    rounded_counts: dict[str, int] = {}
    for value in x_values[fibre]:
        key = f"{float(value):.12f}"
        rounded_counts[key] = rounded_counts.get(key, 0) + 1
    best["affine_fibre_size"] = int(len(fibre))
    best["affine_fibre_group_indices"] = fibre.tolist()
    best["affine_fibre_x_value_counts_rounded_12dp"] = rounded_counts
    best["direct_fibre_logmeanexp"] = direct_certificate

    return {
        "classification": (
            "exhaustive finite effective-quotient enumeration with "
            "floating-point transcendental evaluation"
        ),
        "child_order": n,
        "parent_order": 2 * n,
        "beta": beta,
        "effective_quotient_dimension": rank,
        "switch_coordinate_dimension": dimension,
        "effective_even_even_dimension": 2 * n - 2,
        "bridge_switching_classes_checked": bridge_count,
        "orientations_checked": 2,
        "effective_subspaces_checked_per_orientation": len(subspaces),
        "child_pressure": child_pressure,
        "child_matrix": child.tolist(),
        "best": best,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=4)
    parser.add_argument("--beta", type=float, required=True)
    parser.add_argument("--rank", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = search(args.n, args.beta, args.rank)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
