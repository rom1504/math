"""Exact structural regression of the logarithmic-scale port counterexample.

The asymptotic theorem is analytic. This finite model has three nonzero
amplitude classes and one zero class, so total energy is 3N^2/4. It tests
actual sign blocks, even factors, simultaneous zero-defect placement, and
the pointwise 1/5 inequality using integer arithmetic. No cap optimum is
claimed or needed.
"""
import json
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def regular_hadamard(depth):
    base = np.ones((4, 4), dtype=np.int64) - 2*np.eye(4, dtype=np.int64)
    h = np.ones((1, 1), dtype=np.int64)
    for _ in range(depth):
        h = np.kron(h, base)
    return h


def even_factors(adjacency, count):
    graph = nx.from_numpy_array(adjacency)
    assert nx.is_eulerian(graph)
    directed = list(nx.eulerian_circuit(graph))
    size = len(adjacency)
    bipartite = nx.Graph()
    bipartite.add_nodes_from(range(2*size))
    bipartite.add_edges_from((i, size+j) for i, j in directed)
    factors = []
    for _ in range(count):
        matching = nx.algorithms.bipartite.hopcroft_karp_matching(
            bipartite, top_nodes=set(range(size)))
        edges = [(i, matching[i]-size) for i in range(size)]
        factor = np.zeros((size, size), dtype=np.int64)
        for i, j in edges:
            assert i != j and not factor[i, j]
            factor[i, j] = factor[j, i] = 1
            bipartite.remove_edge(i, size+j)
        assert np.all(factor.sum(axis=1) == 2)
        factors.append(factor)
    assert np.max(sum(factors)) == 1
    return factors


def main():
    small = regular_hadamard(4)
    s = len(small)
    assert np.array_equal(small @ small, s*np.eye(s, dtype=np.int64))
    assert np.all(small.sum(axis=1) == 16)
    assert np.all(np.diag(small) == 1)
    positive = (-small == 1).astype(np.int64)
    assert np.all(positive.sum(axis=1) == 120)
    factors = even_factors(positive, 16)
    t = np.kron(regular_hadamard(1), small)
    n = len(t)
    np.fill_diagonal(t, 0)
    signing = t.copy()
    ports = np.zeros_like(t)
    degrees = [32, 8, 2, 0]
    relative_amplitudes = [1, 2, 4, 0]
    for g, (degree, amplitude) in enumerate(zip(degrees, relative_amplitudes)):
        interval = slice(g*s, (g+1)*s)
        signing[interval, interval] = -1
        if degree:
            ports[interval, interval] = amplitude*sum(factors[:degree//2])
    np.fill_diagonal(signing, 0)
    assert np.all(np.abs(t + np.eye(n, dtype=np.int64)) == 1)
    assert np.all(np.abs(signing + np.eye(n, dtype=np.int64)) == 1)
    # Physical amplitudes are sqrt(32) times these integers.
    assert np.array_equal((ports*ports).sum(axis=1),
                          np.repeat([32, 32, 32, 0], s))
    assert np.all(ports - t*ports.T == 0)
    total_energy_units = int((ports*ports).sum())
    assert total_energy_units*32 == 3*n*n//4

    # Every possible edge-type comparison, not only sampled arrangements.
    for g, a in enumerate(relative_amplitudes):
        for h, b in enumerate(relative_amplitudes):
            signs = [-1] if g == h else [-1, 1]
            for sign in signs:
                assert 5*(a-sign*b)**2 >= a*a+b*b

    rng = np.random.default_rng(202609071934)
    observed = []
    for _ in range(12):
        shuffled = np.zeros_like(ports)
        for i in range(n):
            columns = np.concatenate((np.arange(i), np.arange(i+1, n)))
            shuffled[i, columns] = rng.permutation(ports[i, columns])
        defect_units = int(np.square(shuffled-signing*shuffled.T).sum()//2)
        assert 5*defect_units >= total_energy_units
        observed.append(defect_units)
    output = {
        "status": "PASS", "order": n, "block_order": s,
        "positive_degree": 120, "factor_degrees": degrees,
        "relative_amplitudes": relative_amplitudes,
        "amplitude_square_multiplier": 32,
        "total_port_energy": total_energy_units*32,
        "exact_T_defect": 0,
        "universal_S_defect_lower_bound": {
            "numerator": total_energy_units*32, "denominator": 5},
        "random_permutation_defects": [x*32 for x in observed],
        "seed": 202609071934,
        "scope": "Exact structural identities; no finite cap optimality or asymptotic numerical inference."
    }
    target = ROOT / "computations/results/principle_director_sharp_port_check_2026_09_07.json"
    target.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "order": n,
                      "exact_T_defect": 0, "sampled_permutations": len(observed)}))


if __name__ == "__main__":
    main()
