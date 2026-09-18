#!/usr/bin/env python3
"""Search degree-nonincreasing signed-coordinate restrictions of the seed."""
import importlib.util
import itertools
from pathlib import Path

import numpy as np


def main():
    path = Path(__file__).with_name("bh_boundedness_barrier_2026_09_18_tangent.py")
    spec = importlib.util.spec_from_file_location("exact_tangent", path)
    t = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(t)
    roots = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
    masks = [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 16, 17, 18, 20, 24]
    powers = [0, 2, 0, 5, 1, 4, 2, 4, 3, 1, 4, 1, 2, 2, 3, 2]
    f4 = dict(zip(masks, [roots[j] for j in powers]))
    fpairs = [t.evaluate(f4, x) for x in range(32)]
    root_indices = [roots.index((u // 4, v // 4)) for u, v in fpairs]
    fv = np.exp(1j * np.pi * np.array(root_indices) / 3)
    for n in [1, 2, 3]:
        signs = np.array([[1 - 2 * ((x >> j) & 1) for x in range(1 << n)]
                          for j in range(n)], dtype=int)
        options = np.concatenate([np.ones((1, 1 << n), dtype=int), signs,
                                  -np.ones((1, 1 << n), dtype=int), -signs])
        J = np.zeros((1 << n, 1 << n))
        for x in range(1 << n):
            for j in range(n):
                J[x, x ^ (1 << j)] = 0.5
        best = (-1.0, None, None)
        for mapping in itertools.product(range(2 * (n + 1)), repeat=5):
            inputs = np.zeros(1 << n, dtype=int)
            for j, choice in enumerate(mapping):
                inputs += ((1 - options[choice]) // 2) << j
            values = fv[inputs]
            A = values.conj()[:, None] * J * values[None, :] - J
            top = np.linalg.eigvalsh(A)[-1]
            if top > best[0] + 1e-10:
                best = (float(top), mapping, [root_indices[x] for x in inputs])
        print({"new_variables": n, "best_velocity": best[0],
               "option_indices": best[1], "restricted_sixth_root_powers": best[2]},
              flush=True)


if __name__ == "__main__":
    main()
