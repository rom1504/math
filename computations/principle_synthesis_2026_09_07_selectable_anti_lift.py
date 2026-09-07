"""Exact finite selectable-optimal-child anti-lift diagnostic.

Both the original signing A and skew bridge C are optimization variables.
The child bound is an archived certified M_n; all Boolean child tests and
all directed cut tests are imposed. Solver statuses are computational
certificates, not formal standalone proofs. No asymptotic claim follows.
"""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model


def solve(n, seconds, workers, atlas_index=None):
    archive = Path(__file__).parent / "results" / ("exact_m%d.json" % n)
    payload = json.loads(archive.read_text())
    matrix = np.asarray(payload["matrix"], dtype=np.int64)
    spins = np.asarray([(1,) + s for s in itertools.product((-1, 1), repeat=n-1)],
                       dtype=np.int64)
    seed_cap = int(np.max(np.abs(np.sum(spins * (spins @ matrix), axis=1)))) // 2
    edges = list(itertools.combinations(range(n), 2))
    model = cp_model.CpModel()
    av = [model.new_bool_var("a_%d_%d" % edge) for edge in edges]
    cv = [model.new_bool_var("c_%d_%d" % edge) for edge in edges]
    if atlas_index is not None:
        import networkx as nx
        graph = nx.graph_atlas_g()[atlas_index]
        assert graph.number_of_nodes() == n-1
        if graph.number_of_edges() > (n-1)*(n-2)//4:
            graph = nx.complement(graph)
        root = min(graph.nodes(), key=lambda i: (graph.degree(i), i))
        others = sorted((i for i in graph.nodes() if i != root),
                        key=lambda i: (graph.has_edge(root, i), graph.degree(i), i))
        order = [root] + others
        for e, (i, j) in enumerate(edges):
            value = 0 if i == 0 else int(graph.has_edge(order[i-1], order[j-1]))
            model.add(av[e] == value)
    # Simultaneous switching of both copies gauges A's first row positive.
    for e, (i, j) in enumerate(edges):
        if i == 0:
            model.add(av[e] == 0)
    if n >= 3:
        model.add(av[edges.index((1, 2))] == 0)
    internal = [e for e, (i, j) in enumerate(edges) if i > 0]
    model.add(sum(av[e] for e in internal) <= len(internal) // 2)
    degrees = {i: sum(av[e] for e, edge in enumerate(edges)
                      if i in edge and 0 not in edge) for i in range(1, n)}
    for i in range(2, n):
        model.add(degrees[1] <= degrees[i])
    model.add(degrees[1] <= (n-2)//2)
    for j in range(2, n-1):
        left = av[edges.index((1, j))]
        right = av[edges.index((1, j+1))]
        model.add(left <= right)
        model.add(degrees[j] <= degrees[j+1]).only_enforce_if([left.Not(), right.Not()])
        model.add(degrees[j] <= degrees[j+1]).only_enforce_if([left, right])
    for spin in spins:
        expression = sum(int(spin[i] * spin[j]) * (1 - 2 * av[e])
                         for e, (i, j) in enumerate(edges))
        model.add(expression <= seed_cap)
        model.add(expression >= -seed_cap)

    bound = model.new_int_var(0, len(edges), "directed_bound")
    cut_tests = []
    for subset in range(1, (1 << n) - 1):
        shore = [i for i in range(n) if subset & (1 << i)]
        other = [i for i in range(n) if not subset & (1 << i)]
        free = [i for i in range(n) if i not in (shore[0], other[0])]
        for choices in itertools.product((-1, 1), repeat=len(free)):
            spin = np.ones(n, dtype=np.int64)
            spin[free] = choices
            constant = 0
            ca = np.zeros(len(edges), dtype=np.int64)
            cc = np.zeros(len(edges), dtype=np.int64)
            for e, (i, j) in enumerate(edges):
                i_in = bool(subset & (1 << i))
                j_in = bool(subset & (1 << j))
                if i_in == j_in:
                    continue
                s = int(spin[i] * spin[j])
                ca[e] = -s
                if i_in:
                    cc[e] = s
                else:
                    constant += s
                    cc[e] = -s
            expression = constant + sum(int(ca[e]) * av[e] + int(cc[e]) * cv[e]
                                        for e in range(len(edges)) if ca[e])
            model.add(expression <= bound)
            model.add(expression >= -bound)
            cut_tests.append((constant, ca, cc))
    model.minimize(bound)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = 20260907
    status = solver.solve(model)
    result = dict(n=n, child_cap=seed_cap, status=solver.status_name(status),
                  directed_lower_bound=solver.best_objective_bound,
                  cut_tests=len(cut_tests), wall_time=solver.wall_time,
                  selectable_target=seed_cap / math.sqrt(2),
                  qualification="FINITE SOLVER DIAGNOSTIC; NO ASYMPTOTIC CLAIM")
    if atlas_index is not None:
        result["fixed_atlas_index"] = atlas_index
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        aa = np.asarray([solver.value(v) for v in av], dtype=np.int64)
        cc = np.asarray([solver.value(v) for v in cv], dtype=np.int64)
        observed = max(abs(k + int(u @ aa) + int(v @ cc)) for k, u, v in cut_tests)
        assert observed == solver.value(bound)
        a = np.zeros((n, n), dtype=np.int64)
        c = np.zeros((n, n), dtype=np.int64)
        for e, (i, j) in enumerate(edges):
            a[i, j] = a[j, i] = 1 - 2 * aa[e]
            c[i, j] = 1 - 2 * cc[e]
            c[j, i] = -c[i, j]
        child_actual = int(np.max(np.abs(np.sum(spins * (spins @ a), axis=1)))) // 2
        assert child_actual == seed_cap
        lift = np.block([[a, c], [-c, -a]])
        maximum = 0
        minimum = 0
        # Batch all projective Boolean tests without a large persistent array.
        for first in range(0, 1 << (2*n-1), 4096):
            words = np.arange(first, min(first+4096, 1 << (2*n-1)), dtype=np.int64)
            tails = 1 - 2 * ((words[:, None] >> np.arange(2*n-1)) & 1)
            test = np.concatenate((np.ones((len(words), 1), dtype=np.int64), tails), axis=1)
            energy2 = np.sum(test * (test @ lift), axis=1)
            maximum = max(maximum, int(energy2.max()))
            minimum = min(minimum, int(energy2.min()))
        assert maximum == -minimum == 8 * observed
        result.update(directed_cap=observed, lift_cap=4*observed,
                      child_normalized=seed_cap/n**1.5,
                      lift_normalized=4*observed/(2*n)**1.5,
                      exact_child_a=a.tolist(), skew_c=c.tolist(),
                      independent_full_spin_check=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--orders", type=int, nargs="+", default=[5, 6, 7, 8])
    parser.add_argument("--seconds", type=float, default=30)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--atlas-index", type=int)
    args = parser.parse_args()
    for order in args.orders:
        print(json.dumps(solve(order, args.seconds, args.workers, args.atlas_index)), flush=True)
