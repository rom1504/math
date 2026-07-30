#!/usr/bin/env python3
"""Exhaustively search fixed-order principal restrictions with Gray evaluation."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import random
import time
from collections import Counter
from pathlib import Path

import numpy as np

from conference_prime_square import evaluate


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--matrix-key", default="conference_matrix")
    parser.add_argument("--deletions", type=int, required=True)
    parser.add_argument("--evaluator", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument("--sample-count", type=int)
    parser.add_argument("--seed", type=int, default=20260730)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text())
    matrix = np.asarray(payload[args.matrix_key], dtype=np.int8)
    n = len(matrix)
    choices = list(itertools.combinations(range(n), args.deletions))
    total_restrictions = len(choices)
    if args.sample_count is not None and args.sample_count < len(choices):
        choices = random.Random(args.seed).sample(choices, args.sample_count)
        choices.sort()
    cap_counts: Counter[int] = Counter()
    range_counts: Counter[str] = Counter()
    best_cap = None
    best_count = 0
    best_records = []
    ordered = []
    started = time.monotonic()
    for index, deleted in enumerate(choices, start=1):
        keep = [vertex for vertex in range(n) if vertex not in deleted]
        principal = matrix[np.ix_(keep, keep)]
        record = evaluate(args.evaluator.resolve(), principal)
        cap = int(record["cap"])
        cap_counts[cap] += 1
        range_key = f"min={record['min_energy']},max={record['max_energy']}"
        range_counts[range_key] += 1
        ordered.append([list(deleted), cap, range_key, record["matrix_sha256"]])
        if best_cap is None or cap < best_cap:
            best_cap = cap
            best_count = 0
            best_records = []
        if cap == best_cap:
            best_count += 1
            if len(best_records) < 10:
                best_records.append(
                    {
                        "deleted_vertices": list(deleted),
                        "kept_vertices": keep,
                        "profile": record,
                    }
                )
        if index % args.progress_every == 0 or index == len(choices):
            print(
                f"checked={index}/{len(choices)} best_cap={best_cap} "
                f"cap_counts={dict(sorted(cap_counts.items()))}",
                flush=True,
            )
    serialized = json.dumps(ordered, separators=(",", ":"))
    output = {
        "schema": "quadratic-signing-fixed-restriction-gray-v1",
        "classification": (
            "exhaustive exact finite fixed-signing computation"
            if len(choices) == total_restrictions
            else "exact profiles on a deterministic finite sample; restriction search not exhaustive"
        ),
        "input": str(args.input),
        "matrix_key": args.matrix_key,
        "parent_order": n,
        "deletions": args.deletions,
        "child_order": n - args.deletions,
        "restrictions_checked": len(choices),
        "total_restrictions": total_restrictions,
        "sample_seed": args.seed if len(choices) < total_restrictions else None,
        "best_cap": best_cap,
        "best_count": best_count,
        "best_record_examples": best_records,
        "cap_counts": {str(key): value for key, value in sorted(cap_counts.items())},
        "energy_range_counts": dict(sorted(range_counts.items())),
        "ordered_records_sha256": hashlib.sha256(serialized.encode()).hexdigest(),
        "elapsed_seconds": time.monotonic() - started,
        "evaluator_source": "computations/exact_fixed_signing_gray.cpp",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
