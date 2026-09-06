#!/usr/bin/env python3
"""Search a precise exact-minimizer orientation-gap falsifier.

MILP discovery only: any returned signing is independently enumerated with
integer arithmetic. An infeasible/timeout solver status is not a standalone
proof certificate. Global minimum caps above order eight are imported from
the repository's existing certificates, not reproved by this program.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import csr_matrix, vstack


CAPS = {3: 3, 4: 4, 5: 4, 6: 5, 7: 9, 8: 10, 9: 12,
        10: 13, 11: 17, 12: 18, 13: 20, 14: 21}


def run(n, gap, seconds):
    cap = CAPS[n]
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    d = len(edges)
    states = np.ones((1 << (n - 1), n), dtype=np.int16)
    states[:, 1:] = 1 - 2 * ((np.arange(len(states))[:, None]
                              >> np.arange(n - 1)[None, :]) & 1)
    chars = states[:, [i for i, _ in edges]] * states[:, [j for _, j in edges]]
    shift = chars.sum(axis=1)
    lower = (-cap + gap + shift) / 2
    upper = (cap + shift) / 2
    lower[0] = upper[0]  # switch an actual positive absolute ground to 1.
    rows = [csr_matrix(chars, dtype=float)]
    lows = [lower]
    highs = [upper]

    # Choose vertex zero with smallest signed row sum, then permute the other
    # vertices to sort its incident signs. These two reductions are compatible.
    degree = np.zeros((n, d), dtype=np.int16)
    for e, (i, j) in enumerate(edges):
        degree[i, e] = degree[j, e] = 1
    sym = [degree[j] - degree[0] for j in range(1, n)]
    for j in range(1, n - 1):
        row = np.zeros(d, dtype=np.int16)
        row[edges.index((0, j))] = 1
        row[edges.index((0, j + 1))] = -1
        sym.append(row)
    rows.append(csr_matrix(np.asarray(sym), dtype=float))
    lows.append(np.zeros(len(sym)))
    highs.append(np.full(len(sym), np.inf))
    constraint = LinearConstraint(vstack(rows, format="csr"),
                                  np.concatenate(lows), np.concatenate(highs))
    start = time.monotonic()
    result = milp(np.zeros(d), integrality=np.ones(d),
                  bounds=Bounds(np.zeros(d), np.ones(d)),
                  constraints=constraint,
                  options={"time_limit": seconds, "mip_rel_gap": 0.0,
                           "presolve": True, "disp": False})
    record = {"n": n, "imported_minimum_cap": cap, "required_gap": gap,
              "solver_status": int(result.status), "solver_message": result.message,
              "seconds": time.monotonic() - start,
              "milp_node_count": getattr(result, "mip_node_count", None),
              "scope": "solver discovery; infeasibility is NOT independently certified"}
    if result.x is not None:
        z = np.rint(result.x).astype(np.int16)
        a = 2 * z - 1
        energy = chars @ a
        assert np.all((z == 0) | (z == 1))
        assert int(energy.max()) == cap
        assert int(energy.min()) >= -cap + gap
        matrix = np.zeros((n, n), dtype=np.int16)
        for value, (i, j) in zip(a, edges):
            matrix[i, j] = matrix[j, i] = value
        values, counts = np.unique(energy, return_counts=True)
        record.update(verified_integer_witness=True, matrix=matrix.tolist(),
                      maximum=int(energy.max()), minimum=int(energy.min()),
                      orientation_gap=int(energy.max() + energy.min()),
                      projective_histogram={str(int(v)): int(c) for v, c in zip(values, counts)},
                      matrix_sha256=hashlib.sha256(matrix.astype("i1").tobytes()).hexdigest())
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orders", type=int, nargs="+", default=[9, 10, 11, 12, 13, 14])
    parser.add_argument("--gap", type=int, default=4)
    parser.add_argument("--seconds", type=float, default=120)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    output = {"question": "Can an exact minimum-cap signing have oriented gap >=4?",
              "method": "exact full spin inequalities; SciPy/HiGHS MILP discovery",
              "cases": []}
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    for n in args.orders:
        record = run(n, args.gap, args.seconds)
        output["cases"].append(record)
        path.write_text(json.dumps(output, indent=2) + "\n")
        print(json.dumps({k: v for k, v in record.items()
                          if k not in {"matrix", "projective_histogram"}}), flush=True)


if __name__ == "__main__":
    main()
