#!/usr/bin/env python3
"""Analyze how optimized bridges spend their state-dependent energy budget.

For every feasible bridge in a saved grid, this script groups all projective
state pairs by the aligned internal energy

    u = |H_A(x) + H_B(y)|

and records the realized cross-term magnitudes and slacks

    T - u - |x^T C y|.

It also computes, by exact integer arithmetic, the expected number of violated
constraints for an iid random sign bridge at the same cap.  The statistics are
exact for the saved finite witnesses; they do not constitute an asymptotic
bridge theorem.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np

from exact_mn_milp import projective_spins
from random_bridge_union_bound import strict_absolute_tail_numerator


def energy(matrix: np.ndarray, spins: np.ndarray) -> np.ndarray:
    return (
        np.einsum("bi,ij,bj->b", spins, matrix.astype(np.int64), spins) // 2
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "grid",
        nargs="?",
        type=Path,
        default=Path("computations/results/bridge_grid_through_12.json"),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    grid = json.loads(args.grid.read_text())
    records: list[dict[str, object]] = []
    for row in grid["rows"]:
        if row["status"] not in {"OPTIMAL", "FEASIBLE"}:
            continue
        payload = json.loads(Path(row["result"]).read_text())
        a_payload = json.loads(Path(payload["child_a"]).read_text())
        b_payload = json.loads(Path(payload["child_b"]).read_text())
        a = np.asarray(a_payload["matrix"], dtype=np.int64)
        b = payload["sign_b"] * np.asarray(b_payload["matrix"], dtype=np.int64)
        c = np.asarray(payload["bridge"], dtype=np.int64)
        x = projective_spins(len(a)).astype(np.int64)
        y = projective_spins(len(b)).astype(np.int64)
        ea = energy(a, x)
        eb = energy(b, y)
        internal = np.abs(ea[:, None] + eb[None, :])
        cross = np.abs(x @ c @ y.T)
        cap = int(payload["parent_profile"]["M"])
        slack = cap - internal - cross
        if int(slack.min()) < 0:
            raise AssertionError((row["result"], int(slack.min())))

        levels = []
        internal_hist = Counter(int(value) for value in internal.ravel())
        bridge_variables = len(a) * len(b)
        denominator = 1 << bridge_variables
        union_numerator = 0
        for value in sorted(internal_hist):
            mask = internal == value
            level_cross = cross[mask]
            level_slack = slack[mask]
            state_count = int(mask.sum())
            tail_numerator = strict_absolute_tail_numerator(
                bridge_variables, cap - value
            )
            union_numerator += state_count * tail_numerator
            levels.append(
                {
                    "internal_absolute_energy": value,
                    "state_pair_count": state_count,
                    "allowed_cross_magnitude": cap - value,
                    "realized_cross_maximum": int(level_cross.max()),
                    "realized_cross_mean": float(level_cross.mean()),
                    "active_constraint_count": int(np.count_nonzero(level_slack == 0)),
                    "minimum_slack": int(level_slack.min()),
                    "random_bridge_violation_probability": (
                        tail_numerator / denominator
                    ),
                }
            )

        flat_internal = internal.ravel().astype(float)
        flat_cross = cross.ravel().astype(float)
        correlation = (
            float(np.corrcoef(flat_internal, flat_cross)[0, 1])
            if np.std(flat_internal) and np.std(flat_cross)
            else None
        )
        records.append(
            {
                "m": len(a),
                "n": len(b),
                "sign_b": payload["sign_b"],
                "cap": cap,
                "state_pair_count": int(internal.size),
                "active_constraint_count": int(np.count_nonzero(slack == 0)),
                "minimum_slack": int(slack.min()),
                "internal_cross_pearson_correlation": correlation,
                "iid_random_expected_violation_count": union_numerator / denominator,
                "iid_random_union_numerator": str(union_numerator),
                "iid_random_denominator": str(denominator),
                "levels": levels,
                "source": row["result"],
            }
        )

    output = {
        "schema": "quadratic-signing-bridge-margin-profiles-v1",
        "classification": (
            "exact finite statistics for saved bridge witnesses; random comparison "
            "uses exact union-bound arithmetic; no asymptotic claim"
        ),
        "source": str(args.grid),
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"analyzed {len(records)} feasible bridge witnesses")
    for record in records:
        print(
            f"{record['m']}+{record['n']} sign={record['sign_b']:+d} "
            f"cap={record['cap']} active={record['active_constraint_count']}/"
            f"{record['state_pair_count']} random-Eviol="
            f"{record['iid_random_expected_violation_count']:.6g} "
            f"corr={record['internal_cross_pearson_correlation']}"
        )
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
