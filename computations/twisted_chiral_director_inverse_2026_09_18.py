#!/usr/bin/env python3
"""Exact inverse audit of stored parent matrices; searches are explicitly bounded.

Read tracked signing witnesses, find signed order-four antiautomorphisms via
rooted triangle graphs, then test their transversals for the proposed twist.
No solver/global-optimality assertion is made. Integer Walsh transform checks
each input cap independently. Run using the repository .venv.
"""
import argparse
import hashlib
import itertools
import json
import time
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def walsh_values(a):
    n = len(a)
    out = np.zeros(1 << n, dtype=np.int64)
    for i in range(n):
        for j in range(i):
            out[(1 << i) | (1 << j)] = int(a[i, j])
    step = 1
    while step < len(out):
        blocks = out.reshape(-1, 2 * step)
        left, right = blocks[:, :step].copy(), blocks[:, step:].copy()
        blocks[:, :step] = left + right
        blocks[:, step:] = left - right
        step *= 2
    return out


def cap(a):
    return int(np.abs(walsh_values(a)).max())


def triangle_graph(a, root, sign=1):
    nodes = [i for i in range(len(a)) if i != root]
    g = nx.Graph()
    g.add_nodes_from(nodes)
    for i, j in itertools.combinations(nodes, 2):
        if sign * int(a[root, i]) * int(a[i, j]) * int(a[j, root]) == 1:
            g.add_edge(i, j)
    return g


def orbit_map(a, b):
    """Return p,s with b_ij=s_i s_j a_(p_i,p_j), or certify no map."""
    n = len(a)
    if n == 1:
        return [0], [1]
    gb = triangle_graph(b, 0)
    for root in range(n):
        ga = triangle_graph(a, root)
        matcher = nx.algorithms.isomorphism.GraphMatcher(gb, ga)
        for mapping in matcher.isomorphisms_iter():
            p = [root] + [mapping[i] for i in range(1, n)]
            s = [1] + [int(b[0, i]) * int(a[root, p[i]]) for i in range(1, n)]
            assert np.array_equal(b, np.array(s)[:, None] * a[np.ix_(p, p)] * np.array(s)[None, :])
            return p, s
    return None


def inverse(a, max_structures=1000, seconds=60):
    n2 = len(a)
    assert n2 % 2 == 0 and np.array_equal(a, a.T)
    assert np.all(np.diag(a) == 0)
    assert np.all(np.abs(a + np.eye(n2, dtype=int)) == 1)
    values = walsh_values(a)
    result = {"order": n2, "cap": int(np.abs(values).max()),
              "energy_min": int(values.min()), "energy_max": int(values.max()),
              "input_sha256": hashlib.sha256(a.astype(np.int8).tobytes()).hexdigest(),
              "matrix": a.tolist(), "structures_checked": 0, "transversals_checked": 0,
              "found": [], "search_exhaustive": False}
    # Odd spectral trace is a necessary exact invariant for every chiral form.
    odd_trace = int(np.trace(a @ a @ a))
    result["trace_cube"] = odd_trace
    if odd_trace:
        result["search_exhaustive"] = True
        result["no_chiral_reason"] = "nonzero exact trace(A^3)"
        return result
    start = time.monotonic()
    ga = triangle_graph(a, 0)
    for root in range(1, n2):
        gb = triangle_graph(a, root, -1)
        matcher = nx.algorithms.isomorphism.GraphMatcher(ga, gb)
        for mapping in matcher.isomorphisms_iter():
            if time.monotonic() - start > seconds:
                result["stopping_reason"] = "declared time bound"
                return result
            p = [root] + [mapping[i] for i in range(1, n2)]
            if any(p[p[i]] != i or p[i] == i for i in range(n2)):
                continue
            s = [1] + [-int(a[0, i]) * int(a[root, p[i]]) for i in range(1, n2)]
            if any(s[i] * s[p[i]] != -1 for i in range(n2)):
                continue
            assert np.array_equal(np.array(s)[:, None] * a[np.ix_(p, p)] * np.array(s)[None, :], -a)
            result["structures_checked"] += 1
            pairs = [(i, p[i]) for i in range(n2) if i < p[i]]
            n = len(pairs)
            # Global exchange of halves duplicates the last orientation bit.
            for bits in range(1 << max(0, n - 1)):
                first = [pair[(bits >> k) & 1] for k, pair in enumerate(pairs)]
                second = [p[i] for i in first]
                order = first + second
                signs = np.array([1] * n + [s[i] for i in first], dtype=np.int64)
                dmat = signs[:, None] * a[np.ix_(order, order)] * signs[None, :]
                child, c = dmat[:n, :n], dmat[:n, n:]
                assert np.array_equal(c, c.T) and np.array_equal(dmat[n:, n:], -child)
                b = c.copy()
                np.fill_diagonal(b, 0)
                result["transversals_checked"] += 1
                match = orbit_map(child, b)
                if match is not None:
                    pp, ss = match
                    result["found"].append({"child": child.tolist(), "child_cap": cap(child),
                        "bridge_hollow": b.tolist(), "d": np.diag(c).tolist(),
                        "twist_source_indices": pp, "twist_output_signs": ss,
                        "parent_reordering": order, "parent_switching": signs.tolist(),
                        "chiral_permutation": p, "chiral_signs": s,
                        "parent_cap": result["cap"]})
                    result["stopping_reason"] = "first exact twist witness found"
                    return result
            if result["structures_checked"] >= max_structures:
                result["stopping_reason"] = "declared structure bound"
                return result
    result["search_exhaustive"] = True
    result["stopping_reason"] = "all signed chiral structures and transversals exhausted"
    return result


def matrices_in(obj, key=""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k in {"matrix", "parent_matrix", "representative_matrix", "conference_matrix"} and isinstance(v, list):
                yield key + k, np.array(v, dtype=np.int64)
            elif isinstance(v, (dict, list)):
                yield from matrices_in(v, key + k + "/")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            if isinstance(v, dict):
                yield from matrices_in(v, key + str(i) + "/")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+")
    parser.add_argument("--seconds", type=float, default=60)
    parser.add_argument("--max-structures", type=int, default=1000)
    args = parser.parse_args()
    for name in args.inputs:
        payload = json.loads((ROOT / name).read_text())
        for key, a in matrices_in(payload):
            if len(a) % 2 or len(a) > 20:
                continue
            record = inverse(a, args.max_structures, args.seconds)
            record["source"] = name
            record["source_key"] = key
            print(json.dumps(record, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
