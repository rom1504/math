#!/usr/bin/env python3
"""Audit exact principal-submatrix landing gaps across structured parents.

This driver samples or exhausts subsets deterministically, evaluates every
selected child exactly with ``phase2_subset_caps_gray.cpp``, and reports both
cap and b=M^(2/3) landing gaps wherever M_k is certified.  It is designed to
compare conference and nonconference structured parents, not to certify that
an unexamined subset does not exist.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import random
import subprocess
import time
from collections import Counter
from pathlib import Path

import numpy as np


EXACT_M = {3: 3, 4: 4, 5: 4, 6: 5, 7: 9, 8: 10, 9: 12,
           10: 13, 11: 17, 12: 18, 13: 20, 14: 21}


def matrix_hash(matrix: np.ndarray) -> str:
    return hashlib.sha256(matrix.astype(np.int8).tobytes()).hexdigest()


def choose_subsets(n: int, k: int, samples: int, seed: int) -> tuple[list[tuple[int, ...]], bool]:
    total = math.comb(n, k)
    if samples >= total:
        return list(itertools.combinations(range(n), k)), True
    rng = random.Random(seed)
    selected: set[tuple[int, ...]] = set()
    while len(selected) < samples:
        selected.add(tuple(sorted(rng.sample(range(n), k))))
    return sorted(selected), False


def evaluate_batch(executable: Path, matrix: np.ndarray,
                   subsets: list[tuple[int, ...]]) -> list[tuple[int, int, int]]:
    n, k = len(matrix), len(subsets[0])
    chunks = [f"{n} {k} {len(subsets)}\n"]
    chunks.append(" ".join(map(str, matrix.reshape(-1))) + "\n")
    chunks.extend(" ".join(map(str, subset)) + "\n" for subset in subsets)
    completed = subprocess.run(
        [str(executable)], input="".join(chunks), text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,
    )
    rows = [tuple(map(int, line.split())) for line in completed.stdout.splitlines()]
    if len(rows) != len(subsets):
        raise AssertionError((len(rows), len(subsets), completed.stderr))
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--matrix-key", default="conference_matrix")
    parser.add_argument("--label", required=True)
    parser.add_argument("--child-orders", required=True,
                        help="comma-separated child orders")
    parser.add_argument("--samples", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--evaluator", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.source.read_text())
    matrix = np.asarray(payload[args.matrix_key], dtype=np.int8)
    if not np.array_equal(matrix, matrix.T) or np.any(np.diag(matrix)):
        raise ValueError("matrix is not a symmetric zero-diagonal signing")
    records = []
    started = time.monotonic()
    for k in map(int, args.child_orders.split(",")):
        subsets, exhaustive = choose_subsets(len(matrix), k, args.samples,
                                             args.seed + 1009 * k)
        rows = evaluate_batch(args.evaluator.resolve(), matrix, subsets)
        caps = [row[0] for row in rows]
        best = min(caps)
        examples = [list(subsets[i]) for i, cap in enumerate(caps) if cap == best][:10]
        exact = EXACT_M.get(k)
        record = {
            "child_order": k,
            "total_subsets": math.comb(len(matrix), k),
            "subsets_evaluated": len(subsets),
            "subset_search_exhaustive": exhaustive,
            "all_child_caps_exact": True,
            "best_cap": best,
            "best_examples": examples,
            "cap_histogram": {str(cap): count for cap, count in sorted(Counter(caps).items())},
            "ordered_subset_sha256": hashlib.sha256(
                json.dumps(subsets, separators=(",", ":")).encode()).hexdigest(),
        }
        if exact is not None:
            record.update({
                "certified_M_child": exact,
                "best_cap_gap": best - exact,
                "best_b_gap": best ** (2 / 3) - exact ** (2 / 3),
                "best_b_gap_over_child_order":
                    (best ** (2 / 3) - exact ** (2 / 3)) / k,
            })
        records.append(record)
        print(args.label, k, best, "exhaustive" if exhaustive else "sampled", flush=True)

    output = {
        "schema": "quadratic-signing-phase2-structured-landing-audit-v1",
        "classification": (
            "exact fixed-child caps on deterministic subset sets; minimum is "
            "certified only when subset_search_exhaustive is true"
        ),
        "source": str(args.source),
        "matrix_key": args.matrix_key,
        "label": args.label,
        "parent_order": len(matrix),
        "parent_matrix_sha256": matrix_hash(matrix),
        "seed": args.seed,
        "samples_per_order_cap": args.samples,
        "exact_M_values": EXACT_M,
        "records": records,
        "evaluator_source": "computations/phase2_subset_caps_gray.cpp",
        "elapsed_seconds": time.monotonic() - started,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
