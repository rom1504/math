"""Find a symmetric representative in the 2026 nonregular Hadamard classes.

Read-only source download; an exact candidate is independently checked before
printing. UNKNOWN is not evidence that a symmetric representative is absent.
"""
import argparse
import json
import urllib.request

import numpy as np
from ortools.sat.python import cp_model

SOURCE = "https://www.cs.uleth.ca/~hadi/research/28%20matrices%20and%20excess.txt"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=int, default=0)
    parser.add_argument("--seconds", type=float, default=60)
    args = parser.parse_args()
    raw = urllib.request.urlopen(SOURCE).read().decode()
    lines = [line.split() for line in raw.splitlines() if line.strip()]
    mats = []
    for k in range(0, len(lines), 37):
        cap = int(lines[k][0])
        H = np.array([[1 if x == "1" else -1 for x in row]
                      for row in lines[k + 1:k + 37]], dtype=int)
        assert H.shape == (36, 36)
        assert np.array_equal(H @ H.T, 36 * np.eye(36, dtype=int))
        mats.append((cap, H))
    cap, H = mats[args.index]
    # Column signs normalize row zero; independent equivalence preserves beta.
    H = H * H[0, :]
    n = len(H)
    model = cp_model.CpModel()
    p = [model.new_int_var(0, n - 1, f"p{j}") for j in range(n)]
    model.add_all_different(p)
    # d_j = H[j,p_0], since H[0,p_j]=1 and d_0=1.
    d = [model.new_bool_var(f"d{j}") for j in range(n)]
    for j in range(n):
        model.add_element(p[0], ((1 - H[j]) // 2).tolist(), d[j])
    for i in range(n):
        for j in range(i + 1, n):
            hij = model.new_bool_var(f"h{i}_{j}")
            hji = model.new_bool_var(f"h{j}_{i}")
            model.add_element(p[j], ((1 - H[i]) // 2).tolist(), hij)
            model.add_element(p[i], ((1 - H[j]) // 2).tolist(), hji)
            # Four bit XOR must be zero.
            model.add_bool_xor([hij, d[j], hji, d[i].Not()])
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = args.seconds
    solver.parameters.num_search_workers = 2
    status = solver.solve(model)
    result = {"source": SOURCE, "index": args.index, "published_cap": cap,
              "status": solver.status_name(status), "seconds": solver.wall_time}
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        perm = np.array([solver.value(x) for x in p])
        signs = np.array([1 - 2 * solver.value(x) for x in d])
        K = H[:, perm] * signs
        assert np.array_equal(K, K.T)
        assert np.array_equal(K @ K, 36 * np.eye(36, dtype=int))
        result.update(permutation=perm.tolist(), signs=signs.tolist(),
                      symmetric_hadamard=K.tolist())
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
