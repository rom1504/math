"""Finite regression for the hard multicolor degree/sign compiler.

No asymptotic theorem is inferred from this test.  The program verifies
degree preservation, connectedness, and the final signed incidences.
"""

from __future__ import annotations

import argparse
import json
import math

import numpy as np


def compile_colors(m: int, seed: int) -> tuple[np.ndarray, dict]:
    if (m - 1) % 8:
        raise ValueError("Need m=1 modulo 8 for the chosen even target degrees")
    rng = np.random.default_rng(seed)
    probs = np.array([0.25, 0.25, 0.5])
    colors = np.full((m, m), -1, dtype=np.int8)
    ii, jj = np.triu_indices(m, 1)
    values = rng.choice(3, len(ii), p=probs)
    colors[ii, jj] = values
    colors[jj, ii] = values
    targets = ((m - 1) * probs).astype(int)
    before = np.stack([(colors == r).sum(axis=1) for r in range(3)])
    loads = np.zeros(m, dtype=int)
    internal_loads = np.zeros(m, dtype=int)
    cap = m ** 0.75 / 2
    edits = 0
    steps = [0, 0, 0]

    def recolor(a: int, b: int, new: int) -> None:
        nonlocal edits
        assert a != b
        assert colors[a, b] != new
        colors[a, b] = colors[b, a] = new
        loads[a] += 1
        loads[b] += 1
        edits += 1

    for r in (0, 1):
        delta = targets[r] - (colors == r).sum(axis=1)
        while np.any(delta):
            plus = np.flatnonzero(delta > 0)
            minus = np.flatnonzero(delta < 0)
            eligible = internal_loads < cap
            if plus.size and minus.size:
                u, v = int(plus[0]), int(minus[0])
                candidates = np.flatnonzero(
                    eligible & (colors[u] == 2) & (colors[v] == r)
                )
                if not candidates.size:
                    raise RuntimeError(("no transfer mediator", m, r, u, v))
                w = int(candidates[0])
                recolor(u, w, r)
                recolor(v, w, 2)
                internal_loads[w] += 2
                delta[u] -= 1
                delta[v] += 1
                steps[0] += 1
            else:
                sign = 1 if plus.size else -1
                active = plus if plus.size else minus
                u = int(active[0])
                v = int(active[1]) if active.size > 1 else u
                assert u != v or abs(delta[u]) >= 2
                outer, middle = (2, r) if sign == 1 else (r, 2)
                found = None
                for a in np.flatnonzero(eligible & (colors[u] == outer)):
                    a = int(a)
                    if a == v:
                        continue
                    candidates = np.flatnonzero(
                        eligible & (colors[a] == middle) & (colors[v] == outer)
                    )
                    candidates = candidates[(candidates != u) & (candidates != a)]
                    if candidates.size:
                        found = a, int(candidates[0])
                        break
                if found is None:
                    raise RuntimeError(("no same-sign alternating path", m, r, u, v))
                a, b = found
                recolor(u, a, middle)
                recolor(a, b, outer)
                recolor(b, v, middle)
                internal_loads[a] += 2
                internal_loads[b] += 2
                delta[u] -= sign
                delta[v] -= sign
                steps[1 if sign == 1 else 2] += 1
            assert np.array_equal(delta, targets[r] - (colors == r).sum(axis=1))
        for old in range(r + 1):
            assert np.all((colors == old).sum(axis=1) == targets[old])
    for r in range(3):
        assert np.all((colors == r).sum(axis=1) == targets[r])
    assert np.array_equal(colors, colors.T)
    return colors, {
        "m": m,
        "seed": seed,
        "targets": targets.tolist(),
        "initial_max_degree_error": int(np.max(np.abs(before - targets[:, None]))),
        "edge_recolorings": edits,
        "max_total_incident_changes": int(loads.max()),
        "max_internal_incident_changes": int(internal_loads.max()),
        "internal_load_cap": cap,
        "transfer_two_edge_steps": steps[0],
        "positive_three_edge_steps": steps[1],
        "negative_three_edge_steps": steps[2],
    }


def euler_tour(adjacency: np.ndarray) -> list[int]:
    remaining = adjacency.copy()
    stack = [0]
    reverse = []
    while stack:
        v = stack[-1]
        neighbors = np.flatnonzero(remaining[v])
        if neighbors.size:
            w = int(neighbors[0])
            remaining[v, w] = remaining[w, v] = False
            stack.append(w)
        else:
            reverse.append(stack.pop())
    tour = reverse[::-1]
    assert not remaining.any(), "The color graph must be connected"
    assert len(tour) == int(adjacency.sum() // 2) + 1
    assert tour[0] == tour[-1]
    return tour


def compile_signed_incidence(colors: np.ndarray, seed: int) -> dict:
    m = len(colors)
    rng = np.random.default_rng(seed)
    ii, jj = np.triu_indices(m, 1)
    signing = np.ones((m, m), dtype=np.int8)
    signing[ii, jj] = rng.choice([-1, 1], len(ii))
    signing[jj, ii] = signing[ii, jj]
    incidence = np.zeros((m, m), dtype=np.int8)
    parity_edits = []
    for r, amplitude in ((1, 1), (2, 2)):
        adjacency = colors == r
        tour = euler_tour(adjacency)
        local_signing = signing.copy()
        closure = 1
        for a, b in zip(tour, tour[1:]):
            closure *= -int(local_signing[a, b])
        if closure == -1:
            a, b = tour[0], tour[1]
            local_signing[a, b] *= -1
            local_signing[b, a] *= -1
            parity_edits.append([r, a, b])
        outgoing = 1
        for a, b in zip(tour, tour[1:]):
            incidence[a, b] = amplitude * outgoing
            incoming = int(local_signing[a, b]) * outgoing
            incidence[b, a] = amplitude * incoming
            outgoing = -incoming
        assert outgoing == 1
        for vertex in range(m):
            row = incidence[vertex][colors[vertex] == r]
            assert row.sum() == 0
            assert (row == amplitude).sum() == (row == -amplitude).sum()
    assert np.all(incidence[colors == 0] == 0)
    defects = incidence - signing * incidence.T
    defect_sum = int(np.sum(np.square(defects[ii, jj].astype(int))))
    expected = sum(4 * r * r for r, _, _ in parity_edits)
    assert defect_sum == expected
    return {
        "parity_edits": parity_edits,
        "undirected_squared_defect": defect_sum,
        "all_exact_row_types": True,
        "all_color_graphs_connected": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orders", nargs="+", type=int, default=[257, 513])
    parser.add_argument("--seed", type=int, default=20260907)
    args = parser.parse_args()
    reports = []
    for index, m in enumerate(args.orders):
        colors, report = compile_colors(m, args.seed + index)
        report.update(compile_signed_incidence(colors, args.seed + 100 + index))
        reports.append(report)
    print(json.dumps({"status": "finite regression PASS", "reports": reports}, indent=2))


if __name__ == "__main__":
    main()
