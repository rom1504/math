#!/usr/bin/env python3
"""Exhaustively test bridge Gram/conditional-moment composition states.

For small exact minimizer pairs, enumerate every sign bridge C and compute the
exact parent cap.  Compare increasingly rich nonlocal states:

* the full pair (C C^T, C^T C);
* energy-conditioned marginal variance distributions;
* energy-pair-conditioned cross second moments;
* conditioned second and fourth moments.

Any shared state with different exact parent caps is a rigorous finite
falsifier.  The shellwise maximum |x^T C y| is also checked as the exact
response-envelope fallback; it predicts the cap tautologically and is not
claimed to have a uniform construction lemma.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

from exact_mn_milp import projective_spins


def energy(matrix: np.ndarray, spins: np.ndarray) -> np.ndarray:
    return np.einsum("bi,ij,bj->b", spins, matrix, spins, optimize=True) // 2


def bridge_from_code(code: int, m: int, n: int) -> np.ndarray:
    return np.asarray(
        [1 if code & (1 << index) else -1 for index in range(m * n)],
        dtype=np.int64,
    ).reshape(m, n)


def hash_signature(signature: object) -> str:
    return hashlib.sha256(repr(signature).encode()).hexdigest()


def conditional_scalar_profile(
    energies: np.ndarray, values: np.ndarray
) -> tuple[tuple[int, tuple[tuple[int, int], ...]], ...]:
    records = []
    for level in sorted(set(map(int, energies))):
        histogram = Counter(map(int, values[energies == level]))
        records.append((level, tuple(sorted(histogram.items()))))
    return tuple(records)


def conditional_cross_states(
    ea: np.ndarray,
    eb: np.ndarray,
    cross: np.ndarray,
) -> tuple[object, object, object]:
    second = []
    second_fourth = []
    envelope = []
    for left in sorted(set(map(int, ea))):
        left_mask = ea == left
        for right in sorted(set(map(int, eb))):
            values = cross[np.ix_(left_mask, eb == right)].astype(np.int64)
            count = int(values.size)
            square_sum = int(np.sum(values * values))
            fourth_sum = int(np.sum(values**4))
            maximum = int(np.max(np.abs(values)))
            second.append((left, right, count, square_sum))
            second_fourth.append((left, right, count, square_sum, fourth_sum))
            envelope.append((left, right, maximum))
    return tuple(second), tuple(second_fourth), tuple(envelope)


def analyze_case(
    a: np.ndarray,
    b: np.ndarray,
    sign_b: int,
    max_bridges: int | None,
) -> dict[str, object]:
    m, n = len(a), len(b)
    x = projective_spins(m).astype(np.int64)
    y = projective_spins(n).astype(np.int64)
    ea = energy(a, x).astype(np.int64)
    eb = energy(sign_b * b, y).astype(np.int64)
    total_bridges = 1 << (m * n)
    limit = total_bridges if max_bridges is None else min(total_bridges, max_bridges)
    state_names = (
        "full_gram_pair",
        "conditional_marginal_variances",
        "conditional_cross_second_moments",
        "conditional_cross_second_fourth_moments",
        "conditional_cross_max_envelope",
    )
    tables: dict[str, dict[object, tuple[int, int]]] = {
        name: {} for name in state_names
    }
    collisions: dict[str, dict[str, object]] = {}
    cap_histogram: Counter[int] = Counter()

    for code in range(limit):
        c = bridge_from_code(code, m, n)
        cross = x @ c @ y.T
        parent_cap = int(
            np.max(np.abs(ea[:, None] + eb[None, :]) + np.abs(cross))
        )
        cap_histogram[parent_cap] += 1
        gram_left = c @ c.T
        gram_right = c.T @ c
        left_variance = np.einsum("bi,ij,bj->b", x, gram_left, x)
        right_variance = np.einsum("bi,ij,bj->b", y, gram_right, y)
        second, second_fourth, envelope = conditional_cross_states(
            ea, eb, cross
        )
        signatures = {
            "full_gram_pair": (
                tuple(map(int, gram_left.ravel())),
                tuple(map(int, gram_right.ravel())),
            ),
            "conditional_marginal_variances": (
                conditional_scalar_profile(ea, left_variance),
                conditional_scalar_profile(eb, right_variance),
            ),
            "conditional_cross_second_moments": second,
            "conditional_cross_second_fourth_moments": second_fourth,
            "conditional_cross_max_envelope": envelope,
        }
        predicted = max(
            abs(left + right) + maximum for left, right, maximum in envelope
        )
        if predicted != parent_cap:
            raise AssertionError((predicted, parent_cap))
        for name, signature in signatures.items():
            previous = tables[name].get(signature)
            if previous is None:
                tables[name][signature] = (parent_cap, code)
            elif previous[0] != parent_cap and name not in collisions:
                previous_cap, previous_code = previous
                collisions[name] = {
                    "state_sha256": hash_signature(signature),
                    "first_bridge_code": previous_code,
                    "first_parent_cap": previous_cap,
                    "first_bridge": bridge_from_code(previous_code, m, n).tolist(),
                    "second_bridge_code": code,
                    "second_parent_cap": parent_cap,
                    "second_bridge": c.tolist(),
                }

    states = {}
    for name in state_names:
        cap_sets: dict[object, set[int]] = {}
        for signature, (cap, _code) in tables[name].items():
            cap_sets.setdefault(signature, set()).add(cap)
        states[name] = {
            "distinct_state_count": len(tables[name]),
            "different_cap_collision_found": name in collisions,
            "collision": collisions.get(name),
        }
    if states["conditional_cross_max_envelope"]["different_cap_collision_found"]:
        raise AssertionError("exact response envelope cannot have a cap collision")
    return {
        "orders": [m, n],
        "sign_b": sign_b,
        "bridges_enumerated": limit,
        "all_bridges_enumerated": limit == total_bridges,
        "exact_parent_cap_histogram": {
            str(cap): count for cap, count in sorted(cap_histogram.items())
        },
        "states": states,
    }


def find_fourth_moment_collision(
    a: np.ndarray,
    b: np.ndarray,
    sign_b: int,
) -> dict[str, object]:
    """Stop at the first same-(second,fourth)-state/different-cap bridge."""

    m, n = len(a), len(b)
    x = projective_spins(m).astype(np.int64)
    y = projective_spins(n).astype(np.int64)
    ea = energy(a, x).astype(np.int64)
    eb = energy(sign_b * b, y).astype(np.int64)
    seen: dict[object, tuple[int, int]] = {}
    total = 1 << (m * n)
    for code in range(total):
        c = bridge_from_code(code, m, n)
        cross = x @ c @ y.T
        parent_cap = int(
            np.max(np.abs(ea[:, None] + eb[None, :]) + np.abs(cross))
        )
        _second, signature, _envelope = conditional_cross_states(ea, eb, cross)
        previous = seen.get(signature)
        if previous is None:
            seen[signature] = (parent_cap, code)
        elif previous[0] != parent_cap:
            previous_cap, previous_code = previous
            return {
                "orders": [m, n],
                "sign_b": sign_b,
                "bridges_examined": code + 1,
                "all_bridges_enumerated": False,
                "collision_found": True,
                "state_sha256": hash_signature(signature),
                "first_bridge_code": previous_code,
                "first_parent_cap": previous_cap,
                "first_bridge": bridge_from_code(previous_code, m, n).tolist(),
                "second_bridge_code": code,
                "second_parent_cap": parent_cap,
                "second_bridge": c.tolist(),
            }
    return {
        "orders": [m, n],
        "sign_b": sign_b,
        "bridges_examined": total,
        "all_bridges_enumerated": True,
        "collision_found": False,
        "distinct_state_count": len(seen),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--m3-classes",
        type=Path,
        default=Path("computations/results/m3_minimizer_orbits.json"),
    )
    parser.add_argument(
        "--m4-classes",
        type=Path,
        default=Path("computations/results/m4_minimizer_orbits.json"),
    )
    parser.add_argument(
        "--m5-classes",
        type=Path,
        default=Path("computations/results/m5_minimizer_orbits.json"),
    )
    parser.add_argument(
        "--m6-classes",
        type=Path,
        default=Path("computations/results/m6_minimizer_orbits.json"),
    )
    parser.add_argument("--max-bridges", type=int)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    m3 = json.loads(args.m3_classes.read_text())
    m4 = json.loads(args.m4_classes.read_text())
    m5 = json.loads(args.m5_classes.read_text())
    m6 = json.loads(args.m6_classes.read_text())
    a3 = np.asarray(m3["classes"][0]["representative_matrix"], dtype=np.int64)
    a4 = np.asarray(m4["classes"][0]["representative_matrix"], dtype=np.int64)
    a5 = np.asarray(m5["classes"][0]["representative_matrix"], dtype=np.int64)
    a6 = np.asarray(m6["classes"][0]["representative_matrix"], dtype=np.int64)
    cases = []
    for sign in (1, -1):
        cases.append(analyze_case(a3, a3, sign, args.max_bridges))
        cases.append(analyze_case(a3, a4, sign, args.max_bridges))
        cases.append(analyze_case(a3, a5, sign, args.max_bridges))
        cases.append(analyze_case(a4, a4, sign, args.max_bridges))
    payload = {
        "schema": "quadratic-signing-bridge-gram-response-collision-v1",
        "classification": (
            "exhaustive exact finite collision search when all_bridges_enumerated; "
            "no uniform composition lemma"
        ),
        "cases": cases,
        "heldout_fourth_moment_test": find_fourth_moment_collision(
            a3, a6, 1
        ),
        "conclusion": (
            "Gram and conditional moment states are falsified wherever a recorded "
            "same-state/different-cap collision appears. The conditional maximum "
            "response envelope predicts finite cap exactly but remains an open "
            "uniform construction obligation."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    for case in cases:
        print(
            f"case={case['orders']} sign={case['sign_b']:+d} "
            f"caps={case['exact_parent_cap_histogram']}"
        )
        for name, state in case["states"].items():
            print(
                f"  {name}: states={state['distinct_state_count']} "
                f"collision={state['different_cap_collision_found']}"
            )
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
