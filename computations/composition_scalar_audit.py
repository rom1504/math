#!/usr/bin/env python3
"""Audit scalar approximate-subadditivity on certified exact values.

For b_n = M_n^(2/3), the principal composition target has scalar defect

    d(m,n) = b_(m+n) - b_m - b_n.

This script does not construct a parent signing.  It records the unavoidable
defect in the optimal values themselves, along with the energy-scale slack
above the ideal target (b_m+b_n)^(3/2).  Inputs are explicit exact values and
the output is finite certified evidence only.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


DEFAULT_VALUES = {
    3: 3,
    4: 4,
    5: 4,
    6: 5,
    7: 9,
    8: 10,
    9: 12,
    10: 13,
    11: 17,
    12: 18,
}


def audit(values: dict[int, int]) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for m in sorted(values):
        for n in sorted(values):
            if n < m or m + n not in values:
                continue
            child_sum = values[m] ** (2.0 / 3.0) + values[n] ** (2.0 / 3.0)
            parent_power = values[m + n] ** (2.0 / 3.0)
            ideal_energy = child_sum ** 1.5
            rows.append(
                {
                    "m": m,
                    "n": n,
                    "M_m": values[m],
                    "M_n": values[n],
                    "M_parent": values[m + n],
                    "parent_order": m + n,
                    "two_thirds_child_sum": child_sum,
                    "two_thirds_parent": parent_power,
                    "two_thirds_defect": parent_power - child_sum,
                    "ideal_energy_target": ideal_energy,
                    "energy_slack": values[m + n] - ideal_energy,
                }
            )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--value",
        action="append",
        default=[],
        metavar="N=M",
        help="add or override a certified exact value",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    values = dict(DEFAULT_VALUES)
    for item in args.value:
        n_text, value_text = item.split("=", 1)
        values[int(n_text)] = int(value_text)
    rows = audit(values)
    payload = {
        "schema": "quadratic-signing-composition-scalar-audit-v1",
        "classification": "proved arithmetic consequences of supplied exact M_n values",
        "exact_values": {str(n): values[n] for n in sorted(values)},
        "rows": rows,
    }
    for row in rows:
        print(
            f"{row['m']}+{row['n']}={row['parent_order']}: "
            f"d_2/3={row['two_thirds_defect']:+.12f} "
            f"energy_slack={row['energy_slack']:+.12f}"
        )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
