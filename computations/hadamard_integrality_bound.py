#!/usr/bin/env python3
"""Rigorous integrality upper bound for a symmetric-Hadamard signing cap.

If H is a symmetric Hadamard matrix of even order N and x is a sign vector,
put z_i=x_i(Hx)_i.  Then every z_i is even,

    sum_i z_i^2=N^2,    x^T H x=sum_i z_i.

Dynamic programming maximizes the latter sum over this relaxed integer shell.
For a zero-diagonal signing A with x^T A x=x^T H x, half of that maximum is a
rigorous cap upper bound.  The relaxation ignores compatibility with H and x,
so it need not be sharp.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from exact_mn_milp import stable_matrix_hash


def integer_shell_maximum(order: int) -> tuple[int, list[int]]:
    values = list(range(order % 2, order + 1, 2))
    target = order * order
    states: dict[int, tuple[int, list[int]]] = {0: (0, [])}
    for _ in range(order):
        new: dict[int, tuple[int, list[int]]] = {}
        for square_sum, (linear_sum, chosen) in states.items():
            for value in values:
                next_square = square_sum + value * value
                if next_square > target:
                    continue
                next_linear = linear_sum + value
                old = new.get(next_square)
                if old is None or next_linear > old[0]:
                    new[next_square] = (next_linear, chosen + [value])
        states = new
    if target not in states:
        raise AssertionError("integer shell is unexpectedly empty")
    return states[target]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("construction", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.loads(args.construction.read_text())
    a = np.asarray(payload["parent_matrix"], dtype=np.int64)
    h = np.asarray(payload["symmetric_hadamard_matrix"], dtype=np.int64)
    n = len(a)
    if not np.array_equal(h @ h.T, n * np.eye(n, dtype=np.int64)):
        raise AssertionError("H is not Hadamard")
    if np.trace(h - a) != 0 or not np.array_equal(h - a, np.diag(np.diag(h - a))):
        raise AssertionError("H-A is not a trace-zero diagonal correction")
    maximum, values = integer_shell_maximum(n)
    if maximum % 2:
        raise AssertionError(maximum)
    cap_bound = maximum // 2
    edge_parity = (n * (n - 1) // 2) % 2
    if cap_bound % 2 != edge_parity:
        cap_bound -= 1
    output = {
        "schema": "quadratic-signing-hadamard-integrality-bound-v1",
        "classification": "proved fixed-signing cap upper bound by exact integer dynamic programming",
        "source": str(args.construction),
        "n": n,
        "matrix_sha256": stable_matrix_hash(a),
        "integer_shell_maximum_xTHx": maximum,
        "shell_absolute_values": values,
        "energy_parity": edge_parity,
        "certified_cap_upper_bound": cap_bound,
        "proof": [
            "z_i=x_i(Hx)_i has the parity of N",
            "sum z_i^2=N^2 because H H^T=NI",
            "x^T A x=x^T H x=sum z_i",
            "dynamic programming exactly maximizes sum z_i on the relaxed integer shell",
        ],
    }
    print(
        f"n={n}: integer-shell xTHx<={maximum}, signing cap<={cap_bound}, "
        f"shell={values}"
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
