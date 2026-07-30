#!/usr/bin/env python3
"""Transport the conference 7+7 split to the saved exact order-7 signing.

The order-14 conference certificate contains two cap-9 principal blocks.  This
script exhaustively finds signed-permutation equivalences between those blocks
and the saved exact order-7 representative, transports the bridge, and checks
the resulting parent with diagonal blocks E and -E.  This gives an explicit
cap-21 bridge without relying on a bridge solver.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np

from exact_mn_milp import exact_profile, projective_spins, stable_matrix_hash
from random_bridge_union_bound import strict_absolute_tail_numerator


def equivalence(
    reference: np.ndarray, target: np.ndarray
) -> tuple[tuple[int, ...], np.ndarray, int]:
    n = len(reference)
    for sigma in (1, -1):
        for permutation in itertools.permutations(range(n)):
            permuted = reference[np.ix_(permutation, permutation)]
            switches = np.ones(n, dtype=np.int64)
            switches[1:] = sigma * target[0, 1:] * permuted[0, 1:]
            switched = switches[:, None] * permuted * switches[None, :]
            if np.array_equal(sigma * target, switched):
                return permutation, switches, sigma
    raise RuntimeError("the blocks are not signed-permutation equivalent")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("completion", type=Path)
    parser.add_argument("representative", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    completion = json.loads(args.completion.read_text())
    partition = completion["first_exact_7_7_partition"]
    representative_payload = json.loads(args.representative.read_text())
    e = np.asarray(representative_payload["matrix"], dtype=np.int64)
    left = np.asarray(partition["left_matrix"], dtype=np.int64)
    right = np.asarray(partition["right_matrix"], dtype=np.int64)
    original_bridge = np.asarray(partition["bridge"], dtype=np.int64)
    left_map = equivalence(e, left)
    right_map = equivalence(e, right)
    if left_map[2] != -1 or right_map[2] != 1:
        raise AssertionError((left_map[2], right_map[2]))

    # The saved instance has the right block already equal to E.  Switch the
    # left block, undo its permutation, and then swap the two parent blocks.
    left_permutation, left_switches, _ = left_map
    inverse = np.argsort(np.asarray(left_permutation, dtype=int))
    switched_bridge = left_switches[:, None] * original_bridge
    canonical_bridge = switched_bridge[inverse, :].T
    parent = np.block([[e, canonical_bridge], [canonical_bridge.T, -e]])
    profile = exact_profile(parent)
    if profile["M"] != 21:
        raise AssertionError(profile["M"])
    spins = projective_spins(7).astype(np.int64)
    energies = np.einsum("bi,ij,bj->b", spins, e, spins) // 2
    internal = np.abs(energies[:, None] - energies[None, :])
    cross = np.abs(spins @ canonical_bridge @ spins.T)
    slack = 21 - internal - cross
    if int(slack.min()) < 0:
        raise AssertionError(int(slack.min()))
    denominator = 1 << 49
    union_numerator = 0
    for value in np.unique(internal):
        union_numerator += int(np.count_nonzero(internal == value)) * (
            strict_absolute_tail_numerator(49, 21 - int(value))
        )
    output = {
        "schema": "quadratic-signing-transported-conference-7-7-v1",
        "classification": (
            "proved signed-permutation transport and exhaustive parent Boolean profile"
        ),
        "completion": str(args.completion),
        "representative": str(args.representative),
        "representative_matrix_sha256": stable_matrix_hash(e),
        "left_equivalence": {
            "permutation": list(left_map[0]),
            "switches": [int(value) for value in left_map[1]],
            "global_sign": left_map[2],
        },
        "right_equivalence": {
            "permutation": list(right_map[0]),
            "switches": [int(value) for value in right_map[1]],
            "global_sign": right_map[2],
        },
        "sign_b": -1,
        "bridge": [[int(value) for value in row] for row in canonical_bridge],
        "parent_matrix": [[int(value) for value in row] for row in parent],
        "parent_matrix_sha256": stable_matrix_hash(parent),
        "parent_profile": profile,
        "two_thirds_defect": 21 ** (2.0 / 3.0) - 2 * 9 ** (2.0 / 3.0),
        "margin_profile": {
            "state_pair_count": int(internal.size),
            "active_constraint_count": int(np.count_nonzero(slack == 0)),
            "minimum_slack": int(slack.min()),
            "internal_cross_pearson_correlation": float(
                np.corrcoef(internal.ravel(), cross.ravel())[0, 1]
            ),
            "iid_random_expected_violation_count": union_numerator / denominator,
            "iid_random_union_numerator": str(union_numerator),
            "iid_random_denominator": str(denominator),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"transported exact_m7 + (-exact_m7) to cap {profile['M']} "
        f"defect={output['two_thirds_defect']:+.12f}"
    )
    print(f"parent hash={output['parent_matrix_sha256']}")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
