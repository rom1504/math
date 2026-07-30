#!/usr/bin/env python3
"""Consolidate rigorous cap information for conference-doubled signings.

Each input is a construction JSON from ``conference_double_construction.py``.
The script combines exhaustive profiles, explicit Hadamard eigenvector
certificates, exact integer-shell upper bounds, and explicit CP-SAT witnesses.
Every interval endpoint is independently recomputed or read from a result with
the same matrix hash.  Missing decisions remain intervals.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def load_if_present(path: Path) -> dict[str, object] | None:
    return json.loads(path.read_text()) if path.exists() else None


def parity_floor(bound: float, parity: int) -> int:
    value = math.floor(bound + 1e-10)
    if value % 2 != parity:
        value -= 1
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("constructions", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    records = []
    for path in args.constructions:
        construction = json.loads(path.read_text())
        prime = int(construction["paley_prime"])
        q = int(construction["conference_order"])
        n = int(construction["parent_order"])
        matrix_hash = construction["parent_matrix_sha256"]
        stem = path.with_suffix("")
        cap_run = load_if_present(stem.parent / f"{stem.name}_cap.json")
        shell = load_if_present(stem.parent / f"{stem.name}_integrality.json")
        eigen = load_if_present(stem.parent / f"{stem.name}_eigen.json")
        parent_heuristic = load_if_present(
            stem.parent / f"{stem.name}_parent_heuristic.json"
        )
        child_heuristic = load_if_present(
            stem.parent / f"{stem.name}_child_heuristic.json"
        )

        lower = 0
        upper = parity_floor(
            n * math.sqrt(n) / 2,
            (n * (n - 1) // 2) % 2,
        )
        evidence: list[str] = ["symmetric-Hadamard spectral upper bound"]
        if "parent_profile" in construction:
            exact = int(construction["parent_profile"]["M"])
            lower = upper = exact
            evidence.append("exhaustive parent Boolean profile")
        if cap_run is not None:
            if cap_run["matrix_sha256"] != matrix_hash:
                raise AssertionError(f"cap hash mismatch for {path}")
            if cap_run.get("cap_lower_bound") is not None:
                lower = max(lower, int(cap_run["cap_lower_bound"]))
                evidence.append("explicit CP-SAT energy witness")
            if cap_run.get("certified_cap") is not None:
                lower = upper = int(cap_run["certified_cap"])
                evidence.append("CP-SAT exact fixed-cap optimization")
        if shell is not None:
            if shell["matrix_sha256"] != matrix_hash:
                raise AssertionError(f"integer-shell hash mismatch for {path}")
            upper = min(upper, int(shell["certified_cap_upper_bound"]))
            evidence.append("exact integer-shell upper bound")
        if eigen is not None:
            if eigen["matrix_sha256"] != matrix_hash:
                raise AssertionError(f"eigenvector hash mismatch for {path}")
            if eigen.get("certified_cap") is not None:
                exact = int(eigen["certified_cap"])
                lower = upper = exact
                evidence.append(
                    "explicit Boolean Hadamard eigenvector and spectral bound"
                )
        if parent_heuristic is not None:
            if parent_heuristic["matrix_sha256"] != matrix_hash:
                raise AssertionError(f"parent heuristic hash mismatch for {path}")
            lower = max(lower, int(parent_heuristic["cap_lower_bound"]))
            evidence.append("heuristic parent search with explicit energy witness")
        if lower > upper:
            raise AssertionError((path, lower, upper))

        conference_lower = 0
        conference_upper = parity_floor(
            q * math.sqrt(q - 1) / 2,
            (q * (q - 1) // 2) % 2,
        )
        if "conference_profile" in construction:
            conference_lower = conference_upper = int(
                construction["conference_profile"]["M"]
            )
        if child_heuristic is not None:
            if child_heuristic["matrix_sha256"] != construction[
                "conference_matrix_sha256"
            ]:
                raise AssertionError(f"child heuristic hash mismatch for {path}")
            conference_lower = max(
                conference_lower, int(child_heuristic["cap_lower_bound"])
            )
            evidence.append("heuristic child search with explicit energy witness")
        defect_interval = [
            lower ** (2.0 / 3.0) - 2 * conference_upper ** (2.0 / 3.0),
            upper ** (2.0 / 3.0) - 2 * conference_lower ** (2.0 / 3.0),
        ]
        records.append(
            {
                "paley_prime": prime,
                "conference_order": q,
                "conference_cap_interval": [conference_lower, conference_upper],
                "parent_order": n,
                "parent_cap_interval": [lower, upper],
                "parent_normalized_interval": [
                    lower / n**1.5,
                    upper / n**1.5,
                ],
                "two_thirds_defect_interval": defect_interval,
                "matrix_sha256": matrix_hash,
                "evidence": evidence,
                "source": str(path),
            }
        )

    output = {
        "schema": "quadratic-signing-conference-double-family-audit-v1",
        "classification": (
            "proved/certified fixed-family intervals assembled from explicit "
            "certificates; absent endpoints remain open"
        ),
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    for record in records:
        print(
            f"p={record['paley_prime']} q={record['conference_order']} "
            f"child={record['conference_cap_interval']} parent={record['parent_cap_interval']} "
            f"defect={record['two_thirds_defect_interval']}"
        )
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
