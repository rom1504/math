"""Exact bounded screen: one-edge local optimality versus oriented balance."""
import itertools
import json
import argparse
from pathlib import Path

import numpy as np


def scan(n):
    edges = list(itertools.combinations(range(n), 2))
    residual = [e for e in edges if e[0] != 0]
    masks = np.arange(1 << len(residual), dtype=np.uint32)
    signs = np.ones((len(masks), len(edges)), dtype=np.int16)
    for k, edge in enumerate(residual):
        signs[:, edges.index(edge)] = 1 - 2 * ((masks >> k) & 1).astype(np.int16)
    x = np.asarray([(1,) + z for z in itertools.product((-1, 1), repeat=n-1)], dtype=np.int16)
    features = np.asarray([x[:, i] * x[:, j] for i, j in edges], dtype=np.int16)
    energy = signs @ features
    positive = energy.max(axis=1)
    negative = -energy.min(axis=1)
    cap = np.maximum(positive, negative)
    exact_minimum = int(cap.min())
    candidate = np.flatnonzero(positive - negative >= 4)
    original_candidate_count = len(candidate)
    for e in range(len(edges)):
        altered = energy[candidate] - 2 * signs[candidate, e, None] * features[e]
        candidate = candidate[altered.max(axis=1) >= positive[candidate]]
        if not len(candidate):
            break
    result = {"order": n, "root_normalized_signings": len(signs),
              "exact_minimum_from_complete_screen": exact_minimum,
              "minimizer_absolute_gap_values": sorted(set(int(z) for z in
                   abs(positive[cap == exact_minimum] - negative[cap == exact_minimum]))),
              "positive_dominant_gap_at_least_4_signings": original_candidate_count,
              "positive_one_edge_local_minima_with_gap_at_least_4": len(candidate)}
    if len(candidate):
        index = min(candidate, key=lambda j: (int(cap[j]), -int(positive[j]-negative[j]), int(j)))
        a = [[0] * n for _ in range(n)]
        for (i, j), coefficient in zip(edges, signs[index]):
            a[i][j] = a[j][i] = int(coefficient)
        neighbor_caps = []
        neighbor_positive = []
        for e in range(len(edges)):
            altered = energy[index] - 2 * signs[index, e] * features[e]
            neighbor_caps.append(int(abs(altered).max()))
            neighbor_positive.append(int(altered.max()))
        assert min(neighbor_caps) >= int(cap[index])
        assert min(neighbor_positive) >= int(positive[index])
        # Independent direct integer enumeration of this exact witness.
        direct = [sum(a[i][j] * int(z[i]) * int(z[j]) for i, j in edges) for z in x]
        assert direct == energy[index].tolist()
        result["local_counterexample"] = {
            "matrix": a, "P": int(positive[index]), "R": int(negative[index]),
            "Q": int(cap[index]), "gap": int(positive[index]-negative[index]),
            "edge_order": edges, "single_edge_neighbor_caps": neighbor_caps,
            "single_edge_neighbor_P": neighbor_positive,
            "global_optimality": "NOT optimal; exhaustive minimum is strictly lower",
        }
        assert result["local_counterexample"]["Q"] > exact_minimum
    return result


def hadamard_local_trap():
    n = 16
    h = np.ones((4, 4), dtype=np.int16) - 2 * np.eye(4, dtype=np.int16)
    a = np.eye(n, dtype=np.int16) - np.kron(h, h)
    assert np.array_equal(a, a.T) and not np.any(np.diag(a))
    edges = list(itertools.combinations(range(n), 2))
    x = np.ones((1 << (n-1), n), dtype=np.int16)
    masks = np.arange(len(x), dtype=np.uint32)
    for k in range(1, n):
        x[:, k] = 1 - 2 * ((masks >> (k-1)) & 1).astype(np.int16)
    contribution = np.asarray([a[i, j] * x[:, i] * x[:, j] for i, j in edges], dtype=np.int16)
    energy = contribution.sum(axis=0, dtype=np.int16)
    assert (int(energy.max()), int(energy.min())) == (40, -24)
    single = energy[None, :] - 2 * contribution
    assert np.all(single.max(axis=1) == 42)
    assert np.all(abs(single).max(axis=1) == 42)
    pair_hist = {}
    for i in range(len(edges)-1):
        changed = single[i][None, :] - 2 * contribution[i+1:]
        assert np.all(changed.max(axis=1) >= 40)
        values, counts = np.unique(abs(changed).max(axis=1), return_counts=True)
        for value, count in zip(values, counts):
            pair_hist[int(value)] = pair_hist.get(int(value), 0) + int(count)
    assert sum(pair_hist.values()) == 7140 and min(pair_hist) == 40
    escaping_edges = [(0, 5), (0, 10), (0, 15)]
    escaping = energy - 2 * contribution[[edges.index(e) for e in escaping_edges]].sum(axis=0, dtype=np.int16)
    assert (int(escaping.max()), int(escaping.min())) == (38, -30)
    return {
        "matrix": a.tolist(), "formula": "I_16-(J_4-2I_4) tensor (J_4-2I_4)",
        "P": 40, "R": 24, "Q": 40, "orientation_gap": 16,
        "projective_states": len(x), "single_flip_checks": 120,
        "every_single_flip_P_and_Q": 42,
        "double_flip_checks": 7140, "double_flip_Q_histogram": pair_hist,
        "three_flip_escape_edges": escaping_edges,
        "three_flip_escape_P": 38, "three_flip_escape_R": 30,
        "minimum_number_of_edges_needed_to_reduce_P_or_Q": 3,
        "global_optimality": "NOT optimal; independently verified order-16 witness has Q=30",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-order", type=int, default=7, choices=range(3, 9))
    args = parser.parse_args()
    cases = [scan(n) for n in range(3, args.max_order + 1)]
    result = {"scope": f"Complete root-normalized enumeration n=3..{args.max_order}; scalable Hadamard local trap is NOT an exact minimizer", "cases": cases,
              "explicit_hadamard_local_trap": hadamard_local_trap()}
    path = Path(__file__).resolve().parent / "results/transfer_adversary_orientation_local_vs_global_2026_09_06.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
