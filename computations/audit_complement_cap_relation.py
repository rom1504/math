#!/usr/bin/env python3
"""Audit whether a small deleted-block cap determines its complement cap."""

from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np

from conference_prime_square import evaluate
from exact_mn_milp import exact_profile


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--matrix-key", default="conference_matrix")
    parser.add_argument("--deletions", type=int, required=True)
    parser.add_argument("--deleted-cap", type=int, required=True)
    parser.add_argument("--evaluator", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--progress-every", type=int, default=250)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text())
    matrix = np.asarray(payload[args.matrix_key], dtype=np.int8)
    n = len(matrix)
    selected = []
    deleted_cap_counts: Counter[int] = Counter()
    for deleted in itertools.combinations(range(n), args.deletions):
        block = matrix[np.ix_(deleted, deleted)]
        cap = int(exact_profile(block)["M"])
        deleted_cap_counts[cap] += 1
        if cap == args.deleted_cap:
            selected.append(deleted)
    pair_counts: Counter[tuple[int, int]] = Counter()
    examples = []
    all_vertices = set(range(n))
    for index, deleted in enumerate(selected, start=1):
        kept = tuple(sorted(all_vertices - set(deleted)))
        complement = matrix[np.ix_(kept, kept)]
        profile = evaluate(args.evaluator.resolve(), complement)
        pair_counts[(args.deleted_cap, int(profile["cap"]))] += 1
        if len(examples) < 10:
            examples.append(
                {
                    "deleted_vertices": list(deleted),
                    "kept_vertices": list(kept),
                    "complement_profile": profile,
                }
            )
        if index % args.progress_every == 0 or index == len(selected):
            print(
                f"checked={index}/{len(selected)} pair_counts={dict(pair_counts)}",
                flush=True,
            )
    output = {
        "schema": "quadratic-signing-complement-cap-audit-v1",
        "classification": "exhaustive exact finite conditional complement audit",
        "input": str(args.input),
        "parent_order": n,
        "deletions": args.deletions,
        "deleted_cap_counts": {
            str(key): value for key, value in sorted(deleted_cap_counts.items())
        },
        "conditioned_deleted_cap": args.deleted_cap,
        "conditioned_subsets_checked": len(selected),
        "pair_counts": {
            f"{key[0]},{key[1]}": value for key, value in sorted(pair_counts.items())
        },
        "examples": examples,
        "evaluator_source": "computations/exact_fixed_signing_gray.cpp",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
