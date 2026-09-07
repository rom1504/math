"""Exact finite audit of extremal-pair cut geometry; no asymptotic inference."""
import itertools
import json


def run():
    matrices = tests = 0
    for n in range(2, 7):
        edges = list(itertools.combinations(range(n), 2))
        free = [(i, j) for i, j in edges if i]
        spins = [(1,) + s for s in itertools.product((-1, 1), repeat=n - 1)]
        for signs in itertools.product((-1, 1), repeat=len(free)):
            a = {(0, j): 1 for j in range(1, n)}
            a.update(dict(zip(free, signs)))
            energies = [sum(a[i, j] * x[i] * x[j] for i, j in edges) for x in spins]
            lo, hi = min(energies), max(energies)
            assert (hi + lo) % 2 == (hi - lo) % 2 == 0
            midpoint, width = (hi + lo) // 2, (hi - lo) // 2
            top = spins[energies.index(hi)]
            bottom = spins[energies.index(lo)]
            relative = [u * v for u, v in zip(top, bottom)]
            switched = {(i, j): a[i, j] * top[i] * top[j] for i, j in edges}
            cross = [(i, j) for i, j in edges if relative[i] != relative[j]]
            assert sum(switched[e] for e in cross) == width
            cross_cap = 0
            for x in spins:
                bridge = sum(switched[i, j] * x[i] * x[j] for i, j in cross)
                full = sum(switched[i, j] * x[i] * x[j] for i, j in edges)
                internal = full - bridge
                assert abs(bridge) + abs(internal - midpoint) <= width
                cross_cap = max(cross_cap, abs(bridge))
                tests += 1
            assert cross_cap == width
            matrices += 1
    print(json.dumps({"status": "PASS", "matrices": matrices,
                      "exact_spin_inequalities": tests,
                      "maximum_order": 6,
                      "scope": "all switching classes at orders 2 through 6"}))


if __name__ == "__main__":
    run()
