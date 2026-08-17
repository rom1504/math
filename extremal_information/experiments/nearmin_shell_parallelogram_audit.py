#!/usr/bin/env python3
"""Audit the first cap-relative product-closed four-pole shell.

For projective Boolean spins, coordinatewise multiplication is XOR in
F_2^(n-1).  Four distinct spins form the three-port majority-selector orbit
exactly when their product is one, equivalently when two distinct unordered
pairs have the same XOR.  Sorting spins by cap deficit and inserting them
incrementally therefore finds the exact minimum shell width containing a
nondegenerate product parallelogram.

The repository blind-audit JSON supplies exact and one-step-near matrices.
Fresh random controls are deterministic under the recorded seed.  This is a
finite post-hypothesis audit, not an asymptotic theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "extremal_information/experiments/nearmin_blind_structural_results.json"
DEFAULT_OUTPUT = ROOT / "extremal_information/experiments/nearmin_shell_parallelogram_results.json"
DEFAULT_SEED = 20260817


def matrix_hash(a: np.ndarray) -> str:
    return hashlib.sha256(a.astype(np.int8).tobytes()).hexdigest()


def projective_spins(n: int) -> np.ndarray:
    masks = np.arange(1 << (n - 1), dtype=np.uint32)
    bits = ((masks[:, None] >> np.arange(n - 1, dtype=np.uint32)) & 1).astype(np.int8)
    spins = np.ones((len(masks), n), dtype=np.int8)
    spins[:, 1:] = 1 - 2 * bits
    return spins


def energies(a: np.ndarray, spins: np.ndarray) -> np.ndarray:
    return (np.einsum("bi,ij,bj->b", spins, a, spins, optimize=True) // 2).astype(np.int64)


def first_parallelogram(deficits: np.ndarray) -> dict[str, object]:
    order = np.argsort(deficits, kind="stable")
    pair_for_sum: dict[int, tuple[int, int]] = {}
    inserted: list[int] = []
    for rank, raw_u in enumerate(order):
        u = int(raw_u)
        for v in inserted:
            product = u ^ v
            old = pair_for_sum.get(product)
            if old is not None:
                a, b = old
                # Since u is new, equality of XORs forces disjoint pairs
                # unless the pairs are identical.  Keep the check explicit.
                if len({u, v, a, b}) == 4:
                    threshold = int(deficits[u])
                    shell_count = int(np.count_nonzero(deficits <= threshold))
                    return {
                        "threshold": threshold,
                        "quadruple_masks": [a, b, u, v],
                        "shell_count_at_threshold": shell_count,
                        "insertion_rank": rank + 1,
                    }
            else:
                pair_for_sum[product] = (u, v)
        inserted.append(u)
    raise RuntimeError("the complete projective group must contain a parallelogram for n >= 3")


def audit_matrix(a: np.ndarray, label: str, source: str) -> dict[str, object]:
    n = int(a.shape[0])
    spins = projective_spins(n)
    vals = energies(a, spins)
    cap = int(np.max(np.abs(vals)))
    deficits = cap - np.abs(vals)
    para = first_parallelogram(deficits)
    quad = [int(x) for x in para["quadruple_masks"]]
    assert quad[0] ^ quad[1] ^ quad[2] ^ quad[3] == 0
    assert len(set(quad)) == 4
    para["quadruple_deficits"] = [int(deficits[x]) for x in quad]
    para["quadruple_abs_energies"] = [int(abs(vals[x])) for x in quad]
    return {
        "n": n,
        "label": label,
        "source": source,
        "matrix_sha256": matrix_hash(a),
        "matrix": a.astype(int).tolist(),
        "cap": cap,
        "threshold_over_n32": float(para["threshold"] / (n ** 1.5)),
        "shell_counts": {
            "deficit_0": int(np.count_nonzero(deficits <= 0)),
            "deficit_2": int(np.count_nonzero(deficits <= 2)),
            "deficit_4": int(np.count_nonzero(deficits <= 4)),
        },
        "first_product_parallelogram": para,
    }


def unique_records(records: list[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    seen: set[str] = set()
    for rec in records:
        a = np.asarray(rec["matrix"], dtype=np.int8)
        key = matrix_hash(a)
        if key not in seen:
            seen.add(key)
            out.append(rec)
    return out


def summarize(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[tuple[int, str], list[int]] = {}
    for row in rows:
        key = (int(row["n"]), str(row["label"]))
        groups.setdefault(key, []).append(int(row["first_product_parallelogram"]["threshold"]))
    ans = []
    for (n, label), values in sorted(groups.items()):
        arr = np.asarray(values, dtype=float)
        ans.append({
            "n": n,
            "label": label,
            "count": len(values),
            "threshold_min": int(np.min(arr)),
            "threshold_median": float(np.median(arr)),
            "threshold_max": int(np.max(arr)),
            "normalized_median": float(np.median(arr) / (n ** 1.5)),
        })
    return ans


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--random-per-order", type=int, default=12)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    rows: list[dict[str, object]] = []

    exact = unique_records(payload["repository_exact_representatives"])
    near = unique_records(payload["repository_one_step_near_representatives"])
    for rec in exact:
        a = np.asarray(rec["matrix"], dtype=np.int8)
        rows.append(audit_matrix(a, "exact minimizer", str(rec.get("sources", "repository"))))
    for rec in near:
        a = np.asarray(rec["matrix"], dtype=np.int8)
        rows.append(audit_matrix(a, "one-step near", str(rec.get("sources", "repository"))))

    rng = np.random.default_rng(args.seed)
    for n in range(7, 15):
        for idx in range(args.random_per_order):
            upper = rng.choice(np.array([-1, 1], dtype=np.int8), size=(n, n))
            a = np.triu(upper, 1)
            a = a + a.T
            rows.append(audit_matrix(a, "uniform random", f"seed={args.seed};draw={idx}"))

    out = {
        "schema": "nearmin-shell-product-parallelogram-audit-v1",
        "status": "EXACT FINITE COMPUTATION; HYPOTHESIS AUDIT ONLY",
        "seed": args.seed,
        "definition": (
            "minimum absolute-energy deficit threshold whose projective spin shell "
            "contains four distinct masks with XOR zero"
        ),
        "records": rows,
        "summary": summarize(rows),
    }
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(args.output), "records": len(rows), "summary": out["summary"]}, indent=2))


if __name__ == "__main__":
    main()
