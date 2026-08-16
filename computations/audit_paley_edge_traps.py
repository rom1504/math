#!/usr/bin/env python3
"""Verify that saved Paley extremizers oppose every edge.

If an extremizer of cap Q opposes edge e, flipping e raises its energy to
Q+2.  The one-edge Lipschitz bound then certifies the flipped cap exactly.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def energy(matrix: list[list[int]], spin: list[int]) -> int:
    return sum(
        matrix[i][j] * spin[i] * spin[j]
        for i in range(len(matrix))
        for j in range(i + 1, len(matrix))
    )


def audit(
    path: Path,
    matrix: list[list[int]],
    cap: int,
    signed_spins: list[tuple[int, list[int]]],
) -> dict[str, object]:
    n = len(matrix)
    all_edges = {(i, j) for i in range(n) for j in range(i + 1, n)}
    covered: set[tuple[int, int]] = set()
    disagreement_counts = []
    for sign, spin in signed_spins:
        assert energy(matrix, spin) == sign * cap
        opposed = {
            (i, j)
            for i, j in all_edges
            if sign * matrix[i][j] * spin[i] * spin[j] == -1
        }
        covered.update(opposed)
        disagreement_counts.append(len(opposed))
    assert covered == all_edges
    return {
        "input": str(path.relative_to(ROOT)),
        "input_sha256": digest(path),
        "n": n,
        "cap": cap,
        "projective_extremizers": len(signed_spins),
        "opposed_edges_per_extremizer": {
            "minimum": min(disagreement_counts),
            "maximum": max(disagreement_counts),
        },
        "edges_covered": len(covered),
        "total_edges": len(all_edges),
        "all_edges_covered": True,
        "certified_cap_after_each_single_edge_flip": cap + 2,
    }


def main() -> int:
    path10 = ROOT / "computations/results/conference_order10_gf9.json"
    data10 = json.loads(path10.read_text())
    profile10 = data10["conference_profile"]
    spins10 = [(1, x) for x in profile10["top_spins"]]
    spins10 += [(-1, x) for x in profile10["bottom_spins"]]

    path26 = ROOT / "computations/results/conference_order26_gf25.json"
    data26 = json.loads(path26.read_text())
    profile26 = data26["profiles"][0]
    n26 = profile26["n"]

    def decode(code: int) -> list[int]:
        return [1] + [
            -1 if code & (1 << (i - 1)) else 1 for i in range(1, n26)
        ]

    spins26 = [
        (1, decode(code)) for code in profile26["maximizer_gray_codes"]
    ]
    spins26 += [
        (-1, decode(code)) for code in profile26["minimizer_gray_codes"]
    ]

    output = {
        "schema": "paley-conference-edge-trap-audit-v1",
        "classification": (
            "exact checks from exhaustive fixed-signing extremizer records"
        ),
        "cases": [
            audit(
                path10,
                data10["conference_matrix"],
                profile10["Q"],
                spins10,
            ),
            audit(
                path26,
                data26["conference_matrix"],
                profile26["cap"],
                spins26,
            ),
        ],
    }
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
