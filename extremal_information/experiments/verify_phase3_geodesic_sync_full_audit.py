#!/usr/bin/env python3
"""Independent finite checks for the full geodesic-synchronization audit.

This verifier deliberately targets claims not settled by the scalar BLR
calculation itself:

* the selector-cube word-length formula in GS.8, exhaustively through h=4;
* the three-coupled-bent-pairs family in GS.11 at both k=4 and k=8;
* the exact distance-to-linear formula and the diametral assertion for those
  two bent instances.

The all-support, all-future checks through ambient dimension three live in
``verify_phase3_vector_blr_response_audit.py`` and are not duplicated here.
"""

from __future__ import annotations

import json
from collections import deque
from pathlib import Path


BLOCKS = (
    (0b1000, 0b0100, 0b0010, 0b0001),
    (0b1100, 0b1000, 0b0001, 0b0011),
    (0b0100, 0b1100, 0b0011, 0b0010),
)


def weight(x: int) -> int:
    return bin(x).count("1")


def rank(vectors: tuple[int, ...], dimension: int) -> int:
    rows = list(vectors)
    result = 0
    for coordinate in reversed(range(dimension)):
        pivot = next(
            (i for i in range(result, len(rows)) if rows[i] >> coordinate & 1),
            None,
        )
        if pivot is None:
            continue
        rows[result], rows[pivot] = rows[pivot], rows[result]
        for i in range(len(rows)):
            if i != result and (rows[i] >> coordinate & 1):
                rows[i] ^= rows[result]
        result += 1
    return result


def distances(support: tuple[int, ...], dimension: int) -> list[int]:
    answer = [-1] * (1 << dimension)
    answer[0] = 0
    queue: deque[int] = deque([0])
    while queue:
        x = queue.popleft()
        for generator in support:
            y = x ^ generator
            if answer[y] == -1:
                answer[y] = answer[x] + 1
                queue.append(y)
    if -1 in answer:
        raise AssertionError("support does not span")
    return answer


def selector_audit(max_h: int = 4) -> dict[str, int]:
    source_query_checks = 0
    pair_metric_checks = 0
    for h in range(1, max_h + 1):
        kernel_dimension = 2 * h
        basis = tuple(1 << i for i in range(kernel_dimension))
        selected = tuple(1 << (kernel_dimension + i) for i in range(h))
        alternatives = tuple(
            selected[i] ^ (1 << (2 * i)) ^ (1 << (2 * i + 1))
            for i in range(h)
        )
        profiles: dict[int, tuple[int, ...]] = {}
        for j_mask in range(1 << h):
            support = basis + selected + tuple(
                alternatives[i] for i in range(h) if j_mask >> i & 1
            )
            metric = distances(support, 3 * h)
            values = []
            for p_mask in range(1 << h):
                target = 0
                for i in range(h):
                    if p_mask >> i & 1:
                        target ^= alternatives[i]
                expected = weight(p_mask) + 2 * weight(p_mask & ~j_mask)
                if metric[target] != expected:
                    raise AssertionError(
                        f"GS.8 formula failed at h={h}, J={j_mask:b}, P={p_mask:b}"
                    )
                values.append(metric[target])
                source_query_checks += 1
            profiles[j_mask] = tuple(values)
        for j_mask, profile_j in profiles.items():
            for k_mask, profile_k in profiles.items():
                observed = max(abs(a - b) for a, b in zip(profile_j, profile_k))
                expected = 2 * max(
                    weight(j_mask & ~k_mask), weight(k_mask & ~j_mask)
                )
                if observed != expected:
                    raise AssertionError(
                        f"GS.8 selector metric failed at h={h}, J={j_mask:b}, K={k_mask:b}"
                    )
                pair_metric_checks += 1
    return {
        "maximum_h": max_h,
        "rooted_query_checks": source_query_checks,
        "selector_pair_metric_checks": pair_metric_checks,
    }


def block_rows(rows4: tuple[int, ...], blocks: int) -> tuple[int, ...]:
    return tuple(
        row << (4 * block)
        for block in range(blocks)
        for row in rows4
    )


def quadratic(q: int, rows: tuple[int, ...]) -> int:
    answer = 0
    for i in range(len(rows)):
        if q >> i & 1:
            for j in range(i + 1, len(rows)):
                answer ^= ((rows[i] >> j) & 1) & ((q >> j) & 1)
    return answer


def bent_audit(k: int, d: int = 6) -> dict[str, object]:
    if k % 4:
        raise ValueError("k must be a positive multiple of four")
    forms = tuple(block_rows(block, k // 4) for block in BLOCKS)
    if any(rank(rows, k) != k for rows in forms):
        raise AssertionError("singular polar form")
    if any(forms[2][i] != (forms[0][i] ^ forms[1][i]) for i in range(k)):
        raise AssertionError("third polar form is not the sum of the first two")

    def f(q: int) -> int:
        if q == 0:
            return 0
        answer = 0
        for i, form in enumerate(forms):
            value = quadratic(q, form)
            answer |= value << (2 * i)
            answer |= (value ^ 1) << (2 * i + 1)
        return answer

    values = tuple(f(q) for q in range(1 << k))
    triple_checks = 0
    for x in range(1, 1 << k):
        for y in range(1, 1 << k):
            if x == y:
                continue
            if weight(values[x] ^ values[y] ^ values[x ^ y]) != 3:
                raise AssertionError(f"triple defect failed at k={k}, x={x}, y={y}")
            triple_checks += 1

    coordinate_minima = []
    best_masks = []
    for coordinate in range(6):
        best = 1 << k
        best_mask = 0
        for mask in range(1 << k):
            mismatch = sum(
                ((values[q] >> coordinate) & 1) != (weight(mask & q) & 1)
                for q in range(1, 1 << k)
            )
            if mismatch < best:
                best = mismatch
                best_mask = mask
        coordinate_minima.append(best)
        best_masks.append(best_mask)
    expected = 3 * ((1 << k) - (1 << (k // 2)) - 1)
    if sum(coordinate_minima) != expected:
        raise AssertionError(
            f"bent distance failed at k={k}: {sum(coordinate_minima)} != {expected}"
        )

    support = tuple(1 << i for i in range(d)) + tuple(
        values[q] | (q << d) for q in range(1, 1 << k)
    )
    metric = distances(support, d + k)
    target = (1 << d) - 1
    if metric[target] != d or max(metric) != d:
        raise AssertionError(f"diametral claim failed at (D,k)=({d},{k})")

    synchronized_errors = [
        weight(
            values[q]
            ^ sum(
                (weight(mask & q) & 1) << coordinate
                for coordinate, mask in enumerate(best_masks)
            )
        )
        for q in range(1 << k)
    ]
    if max(synchronized_errors) > 9:
        raise AssertionError(f"uniform radius-nine conclusion failed at k={k}")

    return {
        "D": d,
        "k": k,
        "triple_checks": triple_checks,
        "coordinate_distances": coordinate_minima,
        "distance_to_linear": sum(coordinate_minima),
        "theorem_formula": expected,
        "diameter": max(metric),
        "diametral_target_distance": metric[target],
        "maximum_selected_linear_error": max(synchronized_errors),
    }


def main() -> None:
    output = {
        "claim": (
            "The GS.8 selector formula is exact, and the GS.11 bent family "
            "has the claimed distance and diameter beyond its smallest seed."
        ),
        "selector_cube": selector_audit(),
        "bent_instances": [bent_audit(4), bent_audit(8)],
    }
    path = Path(__file__).with_name(
        "phase3_geodesic_sync_full_audit_results.json"
    )
    path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
