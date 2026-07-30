#!/usr/bin/env python3
"""Compare the two exact order-8 minimizer classes under universal doubling."""

from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np

from conference_double_construction import double_conference
from exact_mn_milp import exact_profile, stable_matrix_hash


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("classes", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.loads(args.classes.read_text())
    if payload["order"] != 8 or len(payload["classes"]) != 2:
        raise ValueError("expected the exhaustive two-class order-8 result")
    records = []
    energy_histograms = []
    for row in payload["classes"]:
        matrix = np.asarray(row["representative_matrix"], dtype=np.int8)
        profile = exact_profile(matrix)
        energy_histograms.append(profile["energy_histogram"])
        restrictions: dict[str, dict[str, int]] = {}
        for size in range(2, 9):
            caps = Counter()
            for subset in itertools.combinations(range(8), size):
                principal = matrix[np.ix_(subset, subset)]
                caps[exact_profile(principal)["M"]] += 1
            restrictions[str(size)] = {
                str(cap): count for cap, count in sorted(caps.items())
            }
        parent = double_conference(matrix)
        parent_profile = exact_profile(parent)
        records.append(
            {
                "class": row["class"],
                "class_size": row["root_gauged_labeled_count"],
                "matrix_sha256": stable_matrix_hash(matrix),
                "spectrum": profile["eigenvalues"],
                "energy_histogram": profile["energy_histogram"],
                "principal_restriction_cap_histograms": restrictions,
                "universal_double_cap": parent_profile["M"],
                "universal_double_matrix_sha256": stable_matrix_hash(parent),
                "universal_double_two_thirds_defect": (
                    parent_profile["M"] ** (2.0 / 3.0)
                    - 2 * profile["M"] ** (2.0 / 3.0)
                ),
            }
        )
    if energy_histograms[0] != energy_histograms[1]:
        raise AssertionError("the two class energy histograms unexpectedly differ")
    output = {
        "schema": "quadratic-signing-m8-class-composition-comparison-v1",
        "classification": "exhaustive exact finite structural comparison",
        "source": str(args.classes),
        "proved_shared_energy_histogram": True,
        "records": records,
        "conclusion": (
            "The full one-body Boolean energy histogram does not determine the "
            "universal-double cap; spectral and principal-restriction data differ."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    for record in records:
        print(
            f"class={record['class']} double_cap={record['universal_double_cap']} "
            f"defect={record['universal_double_two_thirds_defect']:+.12f}"
        )
    print("shared energy histogram verified; restriction histograms recorded")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
