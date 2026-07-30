#!/usr/bin/env python3
"""Construct PC(p^2+1) and exactly audit representative restrictions.

The field is GF(p)[t]/(t^2-d), where d is the first quadratic nonresidue.
The optional external evaluator must implement exact projective enumeration;
`exact_fixed_signing_gray.cpp` is the repository implementation.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import subprocess
import time
from collections import Counter
from pathlib import Path

import numpy as np


def first_nonresidue(p: int) -> int:
    squares = {a * a % p for a in range(p)}
    return next(a for a in range(2, p) if a not in squares)


class PrimeSquare:
    def __init__(self, p: int):
        self.p = p
        self.d = first_nonresidue(p)
        self.elements = tuple((a, b) for a in range(p) for b in range(p))

    def subtract(self, x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
        return ((x[0] - y[0]) % self.p, (x[1] - y[1]) % self.p)

    def multiply(self, x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
        return (
            (x[0] * y[0] + self.d * x[1] * y[1]) % self.p,
            (x[0] * y[1] + x[1] * y[0]) % self.p,
        )

    def power(self, x: tuple[int, int], exponent: int) -> tuple[int, int]:
        result = (1, 0)
        while exponent:
            if exponent & 1:
                result = self.multiply(result, x)
            x = self.multiply(x, x)
            exponent //= 2
        return result

    def character(self, x: tuple[int, int]) -> int:
        if x == (0, 0):
            return 0
        value = self.power(x, (self.p * self.p - 1) // 2)
        if value == (1, 0):
            return 1
        if value == ((self.p - 1) % self.p, 0):
            return -1
        raise AssertionError((x, value))

    def conference(self) -> np.ndarray:
        q = self.p * self.p
        matrix = np.zeros((q + 1, q + 1), dtype=np.int8)
        matrix[0, 1:] = matrix[1:, 0] = 1
        for i, x in enumerate(self.elements, start=1):
            for j, y in enumerate(self.elements, start=1):
                if i != j:
                    matrix[i, j] = self.character(self.subtract(x, y))
        product = matrix.astype(np.int64) @ matrix.astype(np.int64)
        if not np.array_equal(product, q * np.eye(q + 1, dtype=np.int64)):
            raise AssertionError("conference identity failed")
        return matrix


def stable_hash(matrix: np.ndarray) -> str:
    return hashlib.sha256(matrix.astype(np.int8).tobytes()).hexdigest()


def evaluate(
    executable: Path, matrix: np.ndarray, collect_extremizers: bool = False
) -> dict[str, object]:
    payload = str(len(matrix)) + "\n" + "\n".join(
        " ".join(str(int(value)) for value in row) for row in matrix
    ) + "\n"
    started = time.monotonic()
    command = [str(executable)]
    if collect_extremizers:
        command.append("--collect-extremizers")
    process = subprocess.run(
        command, input=payload, text=True, capture_output=True, check=True
    )
    record = json.loads(process.stdout)
    record["elapsed_seconds"] = time.monotonic() - started
    record["matrix_sha256"] = stable_hash(matrix)
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--evaluator", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--all-deletions",
        action="store_true",
        help="exhaustively profile all one- and two-vertex restrictions",
    )
    args = parser.parse_args()
    field = PrimeSquare(args.p)
    matrix = field.conference()
    restrictions = []
    for deleted in ((), (0,), (0, 1)):
        keep = [i for i in range(len(matrix)) if i not in deleted]
        principal = matrix[np.ix_(keep, keep)]
        record = evaluate(
            args.evaluator.resolve(), principal, collect_extremizers=True
        )
        record["deleted_vertices"] = list(deleted)
        restrictions.append(record)
        print(
            f"order={len(principal)} deleted={deleted} cap={record['cap']} "
            f"range=[{record['min_energy']},{record['max_energy']}] "
            f"seconds={record['elapsed_seconds']:.3f}",
            flush=True,
        )
    all_deletion_profiles = {}
    if args.all_deletions:
        for deletion_count in (1, 2):
            choices = list(itertools.combinations(range(len(matrix)), deletion_count))
            counts: Counter[str] = Counter()
            ordered_records = []
            started = time.monotonic()
            for index, deleted in enumerate(choices, start=1):
                keep = [i for i in range(len(matrix)) if i not in deleted]
                principal = matrix[np.ix_(keep, keep)]
                record = evaluate(args.evaluator.resolve(), principal)
                key = (
                    f"cap={record['cap']},min={record['min_energy']},"
                    f"max={record['max_energy']}"
                )
                counts[key] += 1
                ordered_records.append(
                    [list(deleted), key, record["matrix_sha256"]]
                )
                if index % 25 == 0 or index == len(choices):
                    print(
                        f"all-delete-{deletion_count}: {index}/{len(choices)}",
                        flush=True,
                    )
            serialized = json.dumps(ordered_records, separators=(",", ":"))
            all_deletion_profiles[str(deletion_count)] = {
                "deletion_count": deletion_count,
                "restrictions_checked": len(choices),
                "profile_counts": dict(sorted(counts.items())),
                "ordered_records_sha256": hashlib.sha256(
                    serialized.encode()
                ).hexdigest(),
                "elapsed_seconds": time.monotonic() - started,
            }
    def restrict_code(code: int, parent_order: int, deleted: tuple[int, ...]) -> int:
        spins = [1] + [(-1 if code & (1 << (i - 1)) else 1) for i in range(1, parent_order)]
        kept = [spins[i] for i in range(parent_order) if i not in deleted]
        root = kept[0]
        kept = [root * value for value in kept]
        result = 0
        for i, value in enumerate(kept[1:]):
            if value == -1:
                result |= 1 << i
        return result

    nesting = {}
    parent_order = len(matrix)
    for child in restrictions[1:]:
        deleted = tuple(child["deleted_vertices"])
        for kind in ("maximizer", "minimizer"):
            parent_codes = restrictions[0][f"{kind}_gray_codes"]
            child_codes = set(child[f"{kind}_gray_codes"])
            mapped = [restrict_code(code, parent_order, deleted) for code in parent_codes]
            nesting[f"delete_{'_'.join(map(str, deleted))}_{kind}"] = {
                "parent_extremizers": len(parent_codes),
                "distinct_restrictions": len(set(mapped)),
                "restrictions_that_are_child_extremizers": sum(
                    code in child_codes for code in set(mapped)
                ),
                "child_extremizers": len(child_codes),
            }
    output = {
        "schema": "quadratic-signing-prime-square-conference-v1",
        "classification": (
            "proved finite-field construction and conference identity; "
            "exact exhaustive fixed-signing Gray-code profiles"
        ),
        "field": f"GF({args.p ** 2})=GF({args.p})[t]/(t^2-{field.d})",
        "p": args.p,
        "q": args.p ** 2,
        "conference_order": args.p ** 2 + 1,
        "conference_identity": f"C^2={args.p ** 2}I",
        "conference_matrix_sha256": stable_hash(matrix),
        "conference_matrix": [[int(value) for value in row] for row in matrix],
        "profiles": restrictions,
        "extremizer_nesting": nesting,
        "all_deletion_profiles": all_deletion_profiles,
        "evaluator_source": "computations/exact_fixed_signing_gray.cpp",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
