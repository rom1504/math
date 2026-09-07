"""Adaptive exact Hadamard trades; SAT tests Boolean eigenvectors."""
import argparse
import itertools
import json
import time
from pathlib import Path
import numpy as np
from ortools.sat.python import cp_model
from computations.flatify_construct_2026_09_07_symmetric_hadamard_switch import sylvester


def candidates(h):
    n = len(h)
    groups = {}
    for i, j in itertools.combinations(range(n), 2):
        key = np.packbits(h[i] != h[j]).tobytes()
        groups.setdefault(key, []).append((i, j))
    out = set()
    for pairs in groups.values():
        for ab, cd in itertools.combinations(pairs, 2):
            quad = tuple(sorted(ab+cd))
            if len(set(quad)) != 4:
                continue
            block = h[np.ix_(quad, quad)]
            if int(block.sum()) % 8 == 0:
                out.add(quad)
    return sorted(out)


def trade(h, quad):
    q = list(quad)
    rr = np.ones((4, 4), dtype=np.int64)-2*np.eye(4, dtype=np.int64)
    new = h.copy()
    new[q] = (rr @ new[q])//2
    new[:, q] = (new[:, q] @ rr)//2
    assert np.all(np.abs(new) == 1)
    assert np.array_equal(new, new.T)
    assert np.array_equal(new @ new, len(h)*np.eye(len(h), dtype=np.int64))
    return new


def eigen(h, sigma, seconds):
    n = len(h)
    root = int(np.sqrt(n))
    assert root*root == n
    model = cp_model.CpModel()
    b = [model.new_bool_var(f'b{i}') for i in range(n)]
    model.add(b[0] == 1)
    lam = sigma*root
    for i in range(n):
        model.add(sum(int(h[i,j]-(lam if i == j else 0))*b[j] for j in range(n))
                  == int((h[i].sum()-lam)//2))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = 1
    status = solver.solve(model)
    x = None
    if status in [cp_model.OPTIMAL, cp_model.FEASIBLE]:
        x = np.array([2*solver.value(v)-1 for v in b], dtype=np.int64)
        assert np.array_equal(h @ x, lam*x)
    return solver.status_name(status), x


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, default=64)
    parser.add_argument('--steps', type=int, default=80)
    parser.add_argument('--seconds', type=float, default=2)
    args = parser.parse_args()
    h = sylvester(args.n)
    rng = np.random.default_rng(20260907)
    log = []
    witnesses = []
    start = time.monotonic()
    for step in range(args.steps+1):
        statuses = []
        fresh = []
        for sigma in [1, -1]:
            status, x = eigen(h, sigma, args.seconds)
            statuses.append(status)
            if x is not None:
                fresh.append(x)
        witnesses.extend(fresh)
        options = candidates(h)
        if options and witnesses:
            score = [sum(np.prod(x[list(q)]) == -1 for x in witnesses) for q in options]
            maximum = int(max(score))
            choices = [i for i,s in enumerate(score) if s == maximum]
            selected = options[int(rng.choice(choices))]
        else:
            maximum = 0
            selected = options[int(rng.integers(len(options)))] if options else None
        row = dict(step=step, status=statuses, valid_trades=len(options),
                   archived_witnesses=len(witnesses), killed=maximum,
                   selected=selected, elapsed=time.monotonic()-start)
        log.append(row)
        print(json.dumps(row), flush=True)
        if statuses == ['INFEASIBLE', 'INFEASIBLE'] or selected is None or step == args.steps:
            break
        q = list(selected)
        rr = np.ones((4,4), dtype=np.int64)-2*np.eye(4,dtype=np.int64)
        survivors = []
        for x in witnesses:
            if np.prod(x[q]) == 1:
                y = x.copy()
                y[q] = (rr @ x[q])//2
                survivors.append(y)
        witnesses = survivors
        h = trade(h, selected)
    dest = Path(f'computations/results/flatify_construct_2026_09_07_hadamard_eigen_escape_n{args.n}.json')
    dest.write_text(json.dumps(dict(n=args.n, matrix=h.tolist(), log=log), indent=2)+'\n')


if __name__ == '__main__':
    main()
