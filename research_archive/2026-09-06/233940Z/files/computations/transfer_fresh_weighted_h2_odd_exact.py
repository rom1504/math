#!/usr/bin/env python3
"""Exact linear Boolean feasibility for weighted-H2 odd Walsh saturation.

This is not a floating catalyst optimum. A feasible output is independently
recomputed in integers. INFEASIBLE is a solver-dependent finite claim;
UNKNOWN never constitutes a lower bound. No files are written by this script.
"""

import argparse
import json

from ortools.sat.python import cp_model


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--r", type=int, default=3)
    parser.add_argument("--seconds", type=float, default=180)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=260906)
    parser.add_argument("--log", action="store_true")
    args = parser.parse_args()
    if args.r < 1:
        raise ValueError("r must be positive")
    n = 2 ** (2 * args.r + 1)
    scale = 2 ** args.r
    w = [[1 if bin(i & j).count("1") % 2 == 0 else -1
          for j in range(n)] for i in range(n)]
    model = cp_model.CpModel()
    # A Boolean variable t represents the sign 1-2t.
    f = [model.NewBoolVar(f"f_{j}") for j in range(n)]
    h = [model.NewBoolVar(f"h_{j}") for j in range(n)]
    g1 = [model.NewBoolVar(f"g1_{j}") for j in range(n)]
    g2 = [model.NewBoolVar(f"g2_{j}") for j in range(n)]
    for i in range(n):
        # Keep the un-divided equations so r=1 is represented exactly too.
        model.Add(sum(w[i][j] * (1 - 2 * f[j]) for j in range(n))
                  == scale * (-(1 - 2 * g1[i]) + 2 * (1 - 2 * g2[i])))
        model.Add(sum(w[i][j] * (1 - 2 * h[j]) for j in range(n))
                  == (scale // 2) * ((1 - 2 * g1[i]) + 2 * (1 - 2 * g2[i])))
    for left, right, tag in [(f, h, "in"), (g1, g2, "out")]:
        diff = [model.NewBoolVar(f"d_{tag}_{j}") for j in range(n)]
        for j in range(n):
            model.AddBoolXOr([left[j], right[j], diff[j].Not()])
        model.Add(sum(diff) == n // 8)
    # Positive input correlation supplies an agreeing coordinate. Translate
    # it to zero and use common reversal to arrange f(0)=h(0)=+1.
    model.Add(f[0] == 0)
    model.Add(h[0] == 0)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = args.seconds
    solver.parameters.num_search_workers = args.workers
    solver.parameters.random_seed = args.seed
    solver.parameters.log_search_progress = args.log
    status = solver.Solve(model)
    result = {"r": args.r, "walsh_dimension": 2 * args.r + 1,
              "order": n, "status": solver.StatusName(status),
              "wall_time": solver.WallTime(), "branches": solver.NumBranches(),
              "conflicts": solver.NumConflicts(), "seed": args.seed}
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        arrays = [[1 - 2 * solver.Value(t) for t in family]
                  for family in [f, h, g1, g2]]
        ff, hh, gg1, gg2 = arrays
        wf = [sum(w[i][j] * ff[j] for j in range(n)) for i in range(n)]
        wh = [sum(w[i][j] * hh[j] for j in range(n)) for i in range(n)]
        assert all(wf[i] == scale * (-gg1[i] + 2 * gg2[i]) for i in range(n))
        assert all(wh[i] == scale // 2 * (gg1[i] + 2 * gg2[i]) for i in range(n))
        assert sum(ff[i] * hh[i] for i in range(n)) == 3 * n // 4
        assert sum(gg1[i] * gg2[i] for i in range(n)) == 3 * n // 4
        result.update({"integer_verified": True, "arrays": arrays,
                       "walsh_f": wf, "walsh_h": wh})
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
