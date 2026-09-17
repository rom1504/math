#!/usr/bin/env python3
"""Exact switching-quotient cap/Fourier census; no optimization oracle.

Enumerates all signings with first row +1, computes the exact absolute cap,
and performs an integer Walsh transform. Gauge invariance identifies each
chord character with its unique Eulerian completion by star edges. Thus
the resulting polynomial is the actual independent edge-noise response,
not noise only on the gauge-fixed chords. Analytic conclusions are in the
companion artifact; floating polynomial evaluations are diagnostics.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np


def walsh(values: np.ndarray) -> np.ndarray:
    out = values.astype(np.int64, copy=True)
    width = 1
    while width < len(out):
        blocks = out.reshape(-1, 2 * width)
        left, right = blocks[:, :width].copy(), blocks[:, width:].copy()
        blocks[:, :width], blocks[:, width:] = left + right, left - right
        width *= 2
    return out


def census(n: int, batch: int) -> dict:
    chords = list(itertools.combinations(range(1, n), 2))
    q = len(chords)
    size = 1 << q
    spins = np.array(list(itertools.product((-1, 1), repeat=n - 1)), dtype=np.int16)
    witnesses = np.column_stack((np.ones(len(spins), dtype=np.int16), spins))
    feature = np.array([witnesses[:, i] * witnesses[:, j] for i, j in chords], dtype=np.int16)
    star = spins.sum(axis=1)
    caps = np.empty(size, dtype=np.int16)
    degrees = np.empty(size, dtype=np.uint8)
    for lo in range(0, size, batch):
        indices = np.arange(lo, min(lo + batch, size), dtype=np.uint64)
        bits = ((indices[:, None] >> np.arange(q, dtype=np.uint64)) & 1).astype(np.int16)
        signs = 1 - 2 * bits
        caps[lo:lo + len(indices)] = np.max(np.abs(signs @ feature + star), axis=1)
        parity = np.zeros((len(indices), n - 1), dtype=np.int16)
        for e, (i, j) in enumerate(chords):
            parity[:, i - 1] += bits[:, e]
            parity[:, j - 1] += bits[:, e]
        degrees[lo:lo + len(indices)] = bits.sum(axis=1) + (parity & 1).sum(axis=1)
    transform = walsh(caps)
    assert np.all(transform[degrees % 2 == 1] == 0)
    assert np.all(transform[(degrees > 0) & (degrees < 4)] == 0)
    c4 = transform[degrees == 4]
    assert len(c4) == 3 * (n * (n - 1) * (n - 2) * (n - 3) // 24)
    assert np.all(c4 == c4[0])
    ground = np.flatnonzero(caps == caps.min())
    # The four-cycle sum is another exact Walsh transform, evaluated on signings.
    cycles = walsh((degrees == 4).astype(np.int64))
    selected = sorted(set([0, int(ground[0]), int(ground[np.argmax(cycles[ground])]),
                           int(ground[np.argmin(cycles[ground])])]))
    polynomials = []
    max_degree = int(degrees.max())
    for index in selected:
        characters = np.ones(size, dtype=np.int64)
        mask = index
        bit = 0
        all_indices = np.arange(size, dtype=np.uint64)
        while mask:
            if mask & 1:
                characters *= 1 - 2 * ((all_indices >> bit) & 1).astype(np.int64)
            mask >>= 1
            bit += 1
        numerators = [int(np.sum(transform[degrees == k] * characters[degrees == k]))
                      for k in range(max_degree + 1)]
        assert sum(numerators) == int(caps[index]) * size
        ts = np.linspace(0, 1, 501)
        vals = np.polynomial.polynomial.polyval(ts, np.array(numerators) / size)
        polynomials.append({"chord_index": index, "cap": int(caps[index]),
                            "cycle4_sum": int(cycles[index]),
                            "coefficient_numerators": numerators,
                            "denominator": size,
                            "largest_sampled_increase": float(np.max(np.diff(vals))),
                            "sampled_min": float(vals.min()), "sampled_max": float(vals.max())})
    return {"n": n, "gauge_classes": size, "exact_min_cap": int(caps.min()),
            "minimizer_classes": len(ground), "random_cap_sum": int(caps.sum()),
            "random_cap_denominator": size, "cycle4_fourier_numerator": int(c4[0]),
            "cycle4_fourier_denominator": size,
            "minimizer_cycle4_min": int(cycles[ground].min()),
            "minimizer_cycle4_max": int(cycles[ground].max()),
            "symmetry_support_checks_pass": True, "selected_noise_polynomials": polynomials}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=7)
    parser.add_argument("--batch", type=int, default=2048)
    parser.add_argument("--output", type=Path,
                        default=Path("tmp/paper_portfolio_2026_09_17/director/noise_census.json"))
    args = parser.parse_args()
    result = {"status": "exact integer census; polynomial-grid monotonicity is diagnostic",
              "orders": []}
    for n in range(4, args.max_n + 1):
        row = census(n, args.batch)
        result["orders"].append(row)
        print(json.dumps(row), flush=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
