#!/usr/bin/env python3
"""Certify common-parent coverage of all exact minimizer classes for n=3..8."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


CASES = {
    3: (
        "computations/results/m3_minimizer_orbits.json",
        "computations/results/exact_parent_3_3_partition_analysis.json",
        True,
    ),
    4: (
        "computations/results/m4_minimizer_orbits.json",
        "computations/results/exact_parent_4_4_partition_analysis.json",
        True,
    ),
    5: (
        "computations/results/m5_minimizer_orbits.json",
        "computations/results/exact_parent_5_5_partition_analysis.json",
        True,
    ),
    6: (
        "computations/results/m6_minimizer_orbits.json",
        "computations/results/conference_double_6_6_partition_analysis.json",
        True,
    ),
    7: (
        "computations/results/m7_minimizer_orbits.json",
        "computations/results/conference_7_7_partition_analysis.json",
        True,
    ),
    8: (
        "computations/results/m8_minimizer_orbits.json",
        "computations/results/universal_double_8_8_class1_partition_analysis.json",
        False,
    ),
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("computations/results/small_equal_composition_closure.json"),
    )
    args = parser.parse_args()
    records = []
    for n, (classes_text, partitions_text, parent_known_optimal) in CASES.items():
        classes = json.loads(Path(classes_text).read_text())
        partitions = json.loads(Path(partitions_text).read_text())
        class_ids = sorted(
            row["canonical_orbit_sha256"] for row in classes["classes"]
        )
        parent_ids = sorted(
            row["canonical_orbit_sha256"]
            for row in partitions["class_orbit_sizes"]
        )
        if class_ids != parent_ids:
            raise AssertionError((n, class_ids, parent_ids))
        child_cap = int(partitions["child_cap"])
        parent_cap = int(partitions["parent_cap"])
        records.append(
            {
                "child_order": n,
                "child_cap": child_cap,
                "minimizer_class_count": len(class_ids),
                "root_gauged_minimizer_count": classes[
                    "minimizing_signing_count"
                ],
                "common_parent_order": 2 * n,
                "common_parent_cap": parent_cap,
                "common_parent_known_globally_optimal": parent_known_optimal,
                "exact_child_partitions_in_parent": partitions[
                    "exact_child_partition_count"
                ],
                "all_minimizer_classes_covered": True,
                "two_thirds_defect": (
                    parent_cap ** (2.0 / 3.0)
                    - 2 * child_cap ** (2.0 / 3.0)
                ),
                "class_source": classes_text,
                "partition_source": partitions_text,
            }
        )
    output = {
        "schema": "quadratic-signing-small-equal-composition-closure-v1",
        "classification": (
            "exhaustive exact finite class enumeration and parent-partition verification"
        ),
        "statement": (
            "For every n=3,...,8, one saved order-2n parent contains, up to "
            "switching, permutation, and global sign, every exact order-n "
            "minimizer class on one half and an exact minimizer on the other."
        ),
        "scope": (
            "The parents are globally optimal through n=7; at n=8 only the "
            "explicit cap-32 upper bound is claimed."
        ),
        "records": records,
    }
    output_path = args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    for record in records:
        print(
            f"n={record['child_order']} classes={record['minimizer_class_count']} "
            f"common_parent_cap={record['common_parent_cap']} "
            f"defect={record['two_thirds_defect']:+.12f}"
        )
    print(f"wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
