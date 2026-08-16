#!/usr/bin/env python3
"""Finite falsification checks for mixed-circuit holonomy.

The proof in ``drafts/phase3_mixed_circuit_hierarchy.md`` is analytic.  This
script independently checks its identities on every small lifted family in
a fixed exhaustive range and on deterministic larger samples.  It also
checks the mixed-cycle dimension formula, the all-arity construction, the
sharp nullity amplification examples, and a finite response packing.
"""

from __future__ import annotations

import itertools
import json
import math
import random
from pathlib import Path


def weight(x: int) -> int:
    return bin(x).count("1")


def xor_sum(mask: int, values: tuple[int, ...]) -> int:
    answer = 0
    for i, value in enumerate(values):
        if mask >> i & 1:
            answer ^= value
    return answer


def vector_rank(vectors: tuple[int, ...], dimension: int) -> int:
    rows = list(vectors)
    rank = 0
    for coordinate in reversed(range(dimension)):
        pivot = next(
            (i for i in range(rank, len(rows)) if rows[i] >> coordinate & 1),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(len(rows)):
            if i != rank and (rows[i] >> coordinate & 1):
                rows[i] ^= rows[rank]
        rank += 1
    return rank


def independent_basis(vectors: tuple[int, ...], dimension: int) -> tuple[int, ...]:
    basis: list[int] = []
    current_rank = 0
    for vector in vectors:
        candidate = tuple(basis + [vector])
        new_rank = vector_rank(candidate, dimension)
        if new_rank > current_rank:
            basis.append(vector)
            current_rank = new_rank
    return tuple(basis)


def cycles(q_values: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        mask
        for mask in range(1 << len(q_values))
        if xor_sum(mask, q_values) == 0
    )


def circuits(q_values: tuple[int, ...]) -> tuple[int, ...]:
    answer = []
    for mask in cycles(q_values):
        if mask == 0:
            continue
        submask = (mask - 1) & mask
        minimal = True
        while submask:
            if xor_sum(submask, q_values) == 0:
                minimal = False
                break
            submask = (submask - 1) & mask
        if minimal:
            answer.append(mask)
    return tuple(answer)


def cycle_excess(a_values: tuple[int, ...], q_values: tuple[int, ...]) -> int:
    return max(
        weight(xor_sum(mask, a_values)) - weight(mask)
        for mask in cycles(q_values)
    )


def elimination_length(
    d: int,
    a_values: tuple[int, ...],
    q_values: tuple[int, ...],
    u: int,
    q_target: int,
) -> int | None:
    candidates = [
        weight(mask) + weight(u ^ xor_sum(mask, a_values))
        for mask in range(1 << len(a_values))
        if xor_sum(mask, q_values) == q_target
    ]
    return min(candidates) if candidates else None


def brute_length(
    d: int,
    k: int,
    a_values: tuple[int, ...],
    q_values: tuple[int, ...],
    u: int,
    q_target: int,
) -> int | None:
    generators = tuple(1 << i for i in range(d)) + tuple(
        a | (q << d) for a, q in zip(a_values, q_values)
    )
    target = u | (q_target << d)
    best: int | None = None
    for mask in range(1 << len(generators)):
        if best is not None and weight(mask) >= best:
            continue
        if xor_sum(mask, generators) == target:
            best = weight(mask)
    return best


def exhaustive_identity_audit() -> dict[str, int]:
    families = 0
    endpoint_checks = 0
    circuit_checks = 0
    approximate_checks = 0
    for d in range(1, 4):
        for k in range(0, 3):
            column_count = 1 << (d + k)
            for m in range(0, 5):
                if (d + k) * m > 12:
                    continue
                for columns in itertools.product(range(column_count), repeat=m):
                    a_values = tuple(column & ((1 << d) - 1) for column in columns)
                    q_values = tuple(column >> d for column in columns)
                    delta = cycle_excess(a_values, q_values)
                    observed = brute_length(
                        d, k, a_values, q_values, (1 << d) - 1, 0
                    )
                    if observed != d - delta:
                        raise AssertionError(
                            f"MC.1 failed: d={d}, k={k}, columns={columns}"
                        )

                    # A second endpoint is chosen deterministically from the family.
                    signature = sum((i + 1) * column for i, column in enumerate(columns))
                    u = signature & ((1 << d) - 1)
                    q_target = (signature >> d) & ((1 << k) - 1)
                    direct = brute_length(d, k, a_values, q_values, u, q_target)
                    eliminated = elimination_length(
                        d, a_values, q_values, u, q_target
                    )
                    if direct != eliminated:
                        raise AssertionError("general elimination identity failed")

                    cycle_contracts = all(
                        weight(xor_sum(mask, a_values)) <= weight(mask)
                        for mask in cycles(q_values)
                    )
                    circuit_list = circuits(q_values)
                    circuit_contracts = all(
                        weight(xor_sum(mask, a_values)) <= weight(mask)
                        for mask in circuit_list
                    )
                    if cycle_contracts != circuit_contracts:
                        raise AssertionError("circuit criterion failed")

                    nullity = m - vector_rank(q_values, k)
                    max_circuit_defect = max(
                        (
                            max(
                                0,
                                weight(xor_sum(mask, a_values)) - weight(mask),
                            )
                            for mask in circuit_list
                        ),
                        default=0,
                    )
                    if delta > nullity * max_circuit_defect:
                        raise AssertionError("nullity defect bound failed")
                    families += 1
                    endpoint_checks += 2
                    circuit_checks += 1
                    approximate_checks += 1
    return {
        "lifted_families": families,
        "endpoint_identity_checks": endpoint_checks,
        "circuit_criterion_checks": circuit_checks,
        "approximate_defect_checks": approximate_checks,
    }


def random_identity_audit(seed: int = 20260816, trials: int = 1_000) -> dict[str, int]:
    rng = random.Random(seed)
    endpoint_checks = 0
    for _ in range(trials):
        d = rng.randint(1, 6)
        k = rng.randint(0, 4)
        m = rng.randint(0, max(0, 11 - d))
        a_values = tuple(rng.randrange(1 << d) for _ in range(m))
        q_values = tuple(rng.randrange(1 << k) for _ in range(m))
        for _ in range(3):
            u = rng.randrange(1 << d)
            q_target = rng.randrange(1 << k)
            direct = brute_length(d, k, a_values, q_values, u, q_target)
            eliminated = elimination_length(d, a_values, q_values, u, q_target)
            if direct != eliminated:
                raise AssertionError("random endpoint identity failed")
            endpoint_checks += 1
    return {"seed": seed, "families": trials, "endpoint_checks": endpoint_checks}


def mixed_rank_audit(seed: int = 1701, trials: int = 2_000) -> dict[str, int]:
    rng = random.Random(seed)
    gluing_pattern_checks = 0
    gauge_classification_checks = 0
    for _ in range(trials):
        k = rng.randint(0, 5)
        m = rng.randint(0, 8)
        source_count = rng.randint(1, 4)
        q_values = tuple(rng.randrange(1 << k) for _ in range(m))
        colors = tuple(rng.randrange(source_count) for _ in range(m))
        z = cycles(q_values)
        z_basis = independent_basis(z, m)
        local_cycles = tuple(
            mask
            for mask in z
            if mask == 0
            or len({colors[i] for i in range(m) if mask >> i & 1}) <= 1
        )
        z_local_basis = independent_basis(local_cycles, m)
        rank_all = vector_rank(q_values, k)
        rank_sum = sum(
            vector_rank(
                tuple(q_values[i] for i in range(m) if colors[i] == color), k
            )
            for color in range(source_count)
        )
        kappa = rank_sum - rank_all
        if len(z_basis) - len(z_local_basis) != kappa:
            raise AssertionError("mixed-cycle dimension identity failed")

        zero_holonomy_offsets = {
            d_mask
            for d_mask in range(1 << m)
            if not any(weight(d_mask & cycle) & 1 for cycle in z_basis)
        }
        global_shear_offsets = {
            sum(
                (weight(linear_mask & q_values[i]) & 1) << i
                for i in range(m)
            )
            for linear_mask in range(1 << k)
        }
        if zero_holonomy_offsets != global_shear_offsets:
            raise AssertionError("complete gauge classification failed")
        gauge_classification_checks += len(zero_holonomy_offsets)

        # For W=F_2, enumerate offset differences which vanish locally and
        # count their distinct restrictions to the global cycle space.
        patterns = set()
        for d_mask in range(1 << m):
            if any(weight(d_mask & cycle) & 1 for cycle in z_local_basis):
                continue
            patterns.add(tuple(weight(d_mask & cycle) & 1 for cycle in z_basis))
        if len(patterns) != 1 << kappa:
            raise AssertionError(
                f"gluing class count failed: got {len(patterns)}, kappa={kappa}"
            )
        gluing_pattern_checks += len(patterns)
    return {
        "seed": seed,
        "partitions": trials,
        "scalar_gluing_patterns_enumerated": gluing_pattern_checks,
        "scalar_gauge_patterns_checked": gauge_classification_checks,
    }


def arity_hierarchy_audit(max_r: int = 7) -> dict[str, int]:
    proper_unions = 0
    shear_endpoint_checks = 0
    for r in range(2, max_r + 1):
        d = r + 2
        q_values = tuple(1 << i for i in range(r - 1)) + ((1 << (r - 1)) - 1,)
        a_values = tuple(
            sum(1 << coordinate for coordinate in range(d) if coordinate % r == i)
            for i in range(r)
        )
        if xor_sum((1 << r) - 1, a_values) != (1 << d) - 1:
            raise AssertionError("offset partition does not sum to t")

        for source_mask in range((1 << r) - 1):
            chosen = tuple(i for i in range(r) if source_mask >> i & 1)
            qs = tuple(q_values[i] for i in chosen)
            offsets = tuple(a_values[i] for i in chosen)
            if cycle_excess(offsets, qs) != 0:
                raise AssertionError("proper union was not cycle-contracting")

            # The quotient columns are independent.  The shear sending q_i
            # to a_i identifies every endpoint with the zero-offset model.
            for selected_mask in range(1 << len(chosen)):
                q_target = xor_sum(selected_mask, qs)
                shift = xor_sum(selected_mask, offsets)
                for u in range(1 << d):
                    bad = elimination_length(d, offsets, qs, u ^ shift, q_target)
                    control = elimination_length(
                        d, tuple(0 for _ in offsets), qs, u, q_target
                    )
                    if bad != control:
                        raise AssertionError("proper-union shear isometry failed")
                    shear_endpoint_checks += 1
            proper_unions += 1

        full_length = d - cycle_excess(a_values, q_values)
        control_length = d - cycle_excess(tuple(0 for _ in a_values), q_values)
        if full_length != r or control_length != d:
            raise AssertionError("full arity jump failed")
    return {
        "maximum_arity": max_r,
        "proper_unions": proper_unions,
        "shear_endpoint_checks": shear_endpoint_checks,
    }


def sharp_nullity_audit() -> dict[str, int]:
    checks = 0
    for nullity in range(1, 5):
        for defect in range(0, 4):
            d = nullity * (defect + 2)
            q_values = tuple(
                1 << j for j in range(nullity) for _ in range(2)
            )
            offsets = []
            for j in range(nullity):
                block = ((1 << (defect + 2)) - 1) << (j * (defect + 2))
                offsets.extend((0, block))
            a_values = tuple(offsets)
            delta = cycle_excess(a_values, q_values)
            circuit_defects = tuple(
                max(
                    0,
                    weight(xor_sum(mask, a_values)) - weight(mask),
                )
                for mask in circuits(q_values)
            )
            if delta != nullity * defect:
                raise AssertionError("sharp nullity example failed")
            if max(circuit_defects, default=0) != defect:
                raise AssertionError("sharp circuit defect failed")
            if d - delta != 2 * nullity:
                raise AssertionError("sharp target length failed")
            checks += 1
    return {"parameter_pairs": checks}


def response_packing_audit(d: int = 16, target_size: int = 256) -> dict[str, int]:
    minimum_distance = d // 4
    candidates = (v for v in range(1 << d) if weight(v) >= math.ceil(d / 2))
    packed: list[int] = []
    for candidate in candidates:
        if all(weight(candidate ^ previous) >= minimum_distance for previous in packed):
            packed.append(candidate)
            if len(packed) == target_size:
                break
    if len(packed) < target_size:
        raise AssertionError("finite greedy packing was unexpectedly small")

    witness_checks = 0
    observed_separation = d
    for i, v in enumerate(packed):
        own = min(weight(v), 2)
        for v_prime in packed[i + 1 :]:
            cross = min(weight(v), 2 + weight(v ^ v_prime))
            gap = cross - own
            if gap < min(minimum_distance, math.ceil(d / 2) - 2):
                raise AssertionError("response packing witness failed")
            observed_separation = min(observed_separation, gap)
            witness_checks += 1
    return {
        "kernel_dimension": d,
        "packing_size": len(packed),
        "minimum_hamming_distance": minimum_distance,
        "minimum_witnessed_response_gap": observed_separation,
        "pair_witness_checks": witness_checks,
    }


def main() -> None:
    result = {
        "exhaustive_identity": exhaustive_identity_audit(),
        "random_identity": random_identity_audit(),
        "mixed_rank": mixed_rank_audit(),
        "arity_hierarchy": arity_hierarchy_audit(),
        "sharp_nullity": sharp_nullity_audit(),
        "response_packing": response_packing_audit(),
    }
    destination = Path(__file__).with_name(
        "phase3_mixed_circuit_hierarchy_results.json"
    )
    destination.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
