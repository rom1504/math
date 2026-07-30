#!/usr/bin/env python3
"""Sample equal principal splits and exactly profile both children."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
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
    parser.add_argument("--samples", type=int, required=True)
    parser.add_argument("--seed", type=int, default=20260730)
    parser.add_argument("--target-cap", type=int)
    parser.add_argument("--evaluator", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--progress-every", type=int, default=500)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text())
    matrix = np.asarray(payload[args.matrix_key], dtype=np.int8)
    n = len(matrix)
    if n % 2:
        raise ValueError("parent order must be even")
    half = n // 2
    total = math.comb(n - 1, half - 1)
    if args.samples > total:
        raise ValueError("sample count exceeds the number of rooted splits")
    rng = random.Random(args.seed)
    subsets = set()
    while len(subsets) < args.samples:
        subsets.add(tuple(sorted(rng.sample(range(1, n), half - 1))))
    subsets = sorted(subsets)
    pair_counts: Counter[tuple[int, int]] = Counter()
    target_count = 0
    target_examples = []
    best_max = None
    best_sum = None
    best_examples = []
    ordered = []
    started = time.monotonic()
    all_vertices = set(range(n))
    for index, tail in enumerate(subsets, start=1):
        left = (0,) + tail
        right = tuple(sorted(all_vertices - set(left)))
        left_profile = evaluate(args.evaluator.resolve(), matrix[np.ix_(left, left)])
        right_profile = evaluate(args.evaluator.resolve(), matrix[np.ix_(right, right)])
        caps = (int(left_profile["cap"]), int(right_profile["cap"]))
        pair_counts[tuple(sorted(caps))] += 1
        score = (max(caps), sum(caps))
        if best_max is None or score < (best_max, best_sum):
            best_max, best_sum = score
            best_examples = []
        record = {
            "left_vertices": list(left),
            "right_vertices": list(right),
            "left_profile": left_profile,
            "right_profile": right_profile,
        }
        if score == (best_max, best_sum) and len(best_examples) < 10:
            best_examples.append(record)
        if args.target_cap is not None and caps == (args.target_cap, args.target_cap):
            target_count += 1
            if len(target_examples) < 10:
                target_examples.append(record)
        ordered.append([list(left), caps, left_profile["matrix_sha256"], right_profile["matrix_sha256"]])
        if index % args.progress_every == 0 or index == len(subsets):
            print(
                f"checked={index}/{len(subsets)} best=(max {best_max},sum {best_sum}) "
                f"target_count={target_count}",
                flush=True,
            )
    serialized = json.dumps(ordered, separators=(",", ":"))
    output = {
        "schema": "quadratic-signing-equal-split-gray-sample-v1",
        "classification": "exact child profiles on a deterministic sample; split search not exhaustive",
        "input": str(args.input),
        "parent_order": n,
        "child_order": half,
        "total_rooted_splits": total,
        "samples": len(subsets),
        "seed": args.seed,
        "pair_cap_counts": {
            f"{key[0]},{key[1]}": value for key, value in sorted(pair_counts.items())
        },
        "best_max_cap": best_max,
        "best_sum_cap": best_sum,
        "best_examples": best_examples,
        "target_cap": args.target_cap,
        "both_target_count": target_count,
        "both_target_examples": target_examples,
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
