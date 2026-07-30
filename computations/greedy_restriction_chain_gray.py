#!/usr/bin/env python3
"""Build a reproducible locally best principal-restriction chain."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np

from conference_prime_square import evaluate, stable_hash


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--matrix-key", default="conference_matrix")
    parser.add_argument("--initial-deleted", default="")
    parser.add_argument("--target-order", type=int, required=True)
    parser.add_argument("--evaluator", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text())
    parent = np.asarray(payload[args.matrix_key], dtype=np.int8)
    original_order = len(parent)
    deleted = {
        int(value) for value in args.initial_deleted.split(",") if value.strip()
    }
    kept = [vertex for vertex in range(original_order) if vertex not in deleted]
    records = []
    initial = parent[np.ix_(kept, kept)]
    initial_profile = evaluate(args.evaluator.resolve(), initial)
    records.append(
        {
            "order": len(kept),
            "kept_vertices": kept,
            "matrix_sha256": stable_hash(initial),
            "profile": initial_profile,
            "selection": "supplied initial restriction",
        }
    )
    print(f"order={len(kept)} cap={initial_profile['cap']} initial", flush=True)
    while len(kept) > args.target_order:
        candidates = []
        cap_counts: Counter[int] = Counter()
        for vertex in kept:
            child_kept = [item for item in kept if item != vertex]
            child = parent[np.ix_(child_kept, child_kept)]
            profile = evaluate(args.evaluator.resolve(), child)
            cap_counts[int(profile["cap"])] += 1
            candidates.append((int(profile["cap"]), vertex, child_kept, child, profile))
        cap, vertex, kept, child, profile = min(candidates, key=lambda row: (row[0], row[1]))
        records.append(
            {
                "order": len(kept),
                "deleted_vertex_at_step": vertex,
                "kept_vertices": kept,
                "matrix_sha256": stable_hash(child),
                "profile": profile,
                "candidate_cap_counts": {
                    str(key): value for key, value in sorted(cap_counts.items())
                },
                "selection": "minimum-cap one-vertex child; smallest original vertex breaks ties",
            }
        )
        print(
            f"order={len(kept)} cap={cap} deleted={vertex} "
            f"candidate_cap_counts={dict(sorted(cap_counts.items()))}",
            flush=True,
        )
    output = {
        "schema": "quadratic-signing-greedy-restriction-chain-v1",
        "classification": (
            "exact exhaustive choice at each local step; greedy chain not a globally exhaustive restriction search"
        ),
        "input": str(args.input),
        "original_order": original_order,
        "initial_deleted_vertices": sorted(deleted),
        "target_order": args.target_order,
        "records": records,
        "evaluator_source": "computations/exact_fixed_signing_gray.cpp",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
