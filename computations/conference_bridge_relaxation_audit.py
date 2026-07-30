#!/usr/bin/env python3
"""Audit deterministic norm relaxations for a conference block bridge.

For a symmetric conference matrix split as C=[[A,R],[R.T,B]], exact block
multiplication gives

    A^2 + R R^T = (N-1) I,
    B^2 + R^T R = (N-1) I,
    A R + R B = 0.

This program exhaustively compares the actual Boolean parent cap with two
pointwise upper bounds derived only from these identities: a one-sided
Cauchy--Schwarz bound and a first-order anti-intertwining bound.  It is an
audit of proof loss, not an optimizer.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from exact_mn_milp import projective_spins, stable_matrix_hash


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.loads(args.source.read_text())
    parent = np.asarray(payload["parent_matrix"], dtype=np.int64)
    if len(parent) % 2:
        raise ValueError("the parent must have two equal blocks")
    r = len(parent) // 2
    a = parent[:r, :r]
    bridge = parent[:r, r:]
    b = parent[r:, r:]
    identity = np.eye(r, dtype=np.int64)
    conference_scale = len(parent) - 1
    identities = {
        "A2_plus_RRT": np.array_equal(
            a @ a + bridge @ bridge.T, conference_scale * identity
        ),
        "B2_plus_RTR": np.array_equal(
            b @ b + bridge.T @ bridge, conference_scale * identity
        ),
        "AR_plus_RB": np.array_equal(
            a @ bridge + bridge @ b, np.zeros((r, r), dtype=np.int64)
        ),
    }
    if not all(identities.values()):
        raise AssertionError(identities)

    spins = projective_spins(r).astype(np.float64)
    xa = spins @ a
    yb = spins @ b
    child_a = np.einsum("bi,bi->b", xa, spins) / 2
    child_b = np.einsum("bi,bi->b", yb, spins) / 2
    cross = spins @ bridge @ spins.T
    internal = child_a[:, None] + child_b[None, :]
    actual_absolute = np.abs(internal + cross)

    # ||R^T x||^2 = r(N-1)-||Ax||^2 and similarly on the right.
    rr_left = r * conference_scale - np.einsum("bi,bi->b", xa, xa)
    rr_right = r * conference_scale - np.einsum("bi,bi->b", yb, yb)
    rr_left = np.maximum(rr_left, 0)
    rr_right = np.maximum(rr_right, 0)
    one_left = np.sqrt(r * rr_left)[:, None]
    one_right = np.sqrt(r * rr_right)[None, :]
    cross_bound_one_sided = np.minimum(one_left, one_right)
    parent_bound_one_sided = np.abs(internal) + cross_bound_one_sided

    # Write Ax=alpha*x+u and By=beta*y+v.  AR+RB=0 implies
    # (alpha+beta) x^T R y = -u^T R y - x^T R v.
    alpha = 2 * child_a / r
    beta = 2 * child_b / r
    residual_a = xa - alpha[:, None] * spins
    residual_b = yb - beta[:, None] * spins
    residual_a_norm = np.linalg.norm(residual_a, axis=1)
    residual_b_norm = np.linalg.norm(residual_b, axis=1)
    anti_rhs = (
        residual_a_norm[:, None] * np.sqrt(rr_right)[None, :]
        + np.sqrt(rr_left)[:, None] * residual_b_norm[None, :]
    )
    anti_denominator = np.abs(alpha[:, None] + beta[None, :])
    cross_bound_anti = cross_bound_one_sided.copy()
    usable = anti_denominator > 1e-12
    cross_bound_anti[usable] = np.minimum(
        cross_bound_anti[usable], anti_rhs[usable] / anti_denominator[usable]
    )
    parent_bound_anti = np.abs(internal) + cross_bound_anti

    actual_cap = float(actual_absolute.max())
    one_cap = float(parent_bound_one_sided.max())
    anti_cap = float(parent_bound_anti.max())
    singular_values = np.linalg.svd(bridge, compute_uv=False)
    output = {
        "schema": "quadratic-signing-conference-bridge-relaxation-audit-v1",
        "classification": (
            "proved pointwise inequalities and exhaustive finite evaluation"
        ),
        "source": str(args.source),
        "parent_matrix_sha256": stable_matrix_hash(parent),
        "orders": {"child": r, "parent": 2 * r},
        "conference_scale": conference_scale,
        "verified_block_identities": identities,
        "bridge_singular_values": [float(value) for value in singular_values],
        "boolean_state_pair_count": int(actual_absolute.size),
        "actual_parent_cap": int(actual_cap),
        "bounds": {
            "one_sided_norm": {
                "parent_cap_bound": one_cap,
                "additive_loss_above_actual": one_cap - actual_cap,
            },
            "first_order_anti_intertwining": {
                "parent_cap_bound": anti_cap,
                "additive_loss_above_actual": anti_cap - actual_cap,
                "strictly_improved_pair_count": int(
                    np.count_nonzero(
                        cross_bound_anti + 1e-12 < cross_bound_one_sided
                    )
                ),
                "zero_denominator_pair_count": int(np.count_nonzero(~usable)),
            },
        },
        "conclusion": (
            "The conference identities alone, relaxed through norms and one "
            "anti-intertwining moment, do not recover the Boolean cap."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"verified identities={identities}")
    print(
        f"actual={actual_cap:.12f} one_sided={one_cap:.12f} "
        f"anti_intertwining={anti_cap:.12f}"
    )
    print(
        "losses="
        f"{one_cap - actual_cap:.12f},{anti_cap - actual_cap:.12f} "
        f"improved_pairs={output['bounds']['first_order_anti_intertwining']['strictly_improved_pair_count']}"
    )
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
