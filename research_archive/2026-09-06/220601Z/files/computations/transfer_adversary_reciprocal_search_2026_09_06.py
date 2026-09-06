"""Bounded heuristic lower-defect witnesses, with exact final verification."""

import argparse
import json

import numpy as np

from transfer_adversary_reciprocal_pair_probe_2026_09_06 import hadamard, probe


def search(n, batch, iterations, seed):
    rng = np.random.default_rng(seed)
    h = hadamard(n)
    f = rng.choice([-1, 1], size=(batch, n))
    differing = np.zeros((batch, n), dtype=bool)
    for row in range(batch):
        differing[row, rng.choice(n, n // 8, replace=False)] = True
    g = f * (1 - 2 * differing)
    base = probe(16)["l1_witness"]
    if n >= 16 and n.bit_length() % 2 == 1:
        padding = np.ones(n // 16, dtype=np.int64)
        coords = np.arange(n // 16)
        for shift in range(0, (n // 16).bit_length() - 1, 2):
            padding *= 1 - 2 * (((coords >> shift) & 1) & ((coords >> (shift + 1)) & 1))
        f[:batch // 2] = np.kron(base["f"], padding)
        g[:batch // 2] = np.kron(base["g"], padding)
        differing = f != g
    a = f @ h
    b = g @ h
    values = np.sum(np.abs(a * b - 3 * n // 4), axis=1)
    best = int(np.min(values))
    best_pair = (f[int(np.argmin(values))].copy(), g[int(np.argmin(values))].copy())
    rows = np.arange(batch)
    for step in range(iterations):
        if rng.random() < 0.5:
            positions = rng.integers(n, size=batch)
            next_a = a - 2 * f[rows, positions, None] * h[positions]
            next_b = b - 2 * g[rows, positions, None] * h[positions]
            common = True
        else:
            inside = np.argmax(rng.random((batch, n)) * differing, axis=1)
            outside = np.argmax(rng.random((batch, n)) * ~differing, axis=1)
            next_a = a
            next_b = b - 2 * g[rows, inside, None] * h[inside] - 2 * g[rows, outside, None] * h[outside]
            common = False
        next_values = np.sum(np.abs(next_a * next_b - 3 * n // 4), axis=1)
        temperature = n * (0.2 * (1 - step / iterations) ** 2 + 0.002)
        accepted = (next_values <= values) | (np.log(rng.random(batch)) < (values - next_values) / temperature)
        changed = rows[accepted]
        a[accepted] = next_a[accepted]
        b[accepted] = next_b[accepted]
        values[accepted] = next_values[accepted]
        if common:
            f[changed, positions[accepted]] *= -1
            g[changed, positions[accepted]] *= -1
        else:
            g[changed, inside[accepted]] *= -1
            g[changed, outside[accepted]] *= -1
            differing[changed, inside[accepted]] = False
            differing[changed, outside[accepted]] = True
        idx = int(np.argmin(values))
        if int(values[idx]) < best:
            best = int(values[idx])
            best_pair = (f[idx].copy(), g[idx].copy())
            print(json.dumps({"order": n, "step": step, "integer_l1_numerator": best,
                              "normalized_l1_defect": best / n ** 2}), flush=True)
    f, g = best_pair
    a, b = h @ f, h @ g
    assert int(f @ g) == 3 * n // 4
    assert int(np.sum(np.abs(a * b - 3 * n // 4))) == best
    return {"order": n, "batch": batch, "iterations": iterations, "seed": seed,
            "integer_l1_numerator": best, "normalized_l1_defect": best / n ** 2,
            "f": f.tolist(), "g": g.tolist(), "walsh_f": a.tolist(), "walsh_g": b.tolist()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, default=64)
    parser.add_argument("--batch", type=int, default=256)
    parser.add_argument("--iterations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=260906)
    args = parser.parse_args()
    print(json.dumps(search(args.order, args.batch, args.iterations, args.seed)))
