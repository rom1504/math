#!/usr/bin/env python3
"""Exact exhaustive classification of nonzero perfect dyadic pairs at N<=16.

For each f up to common reversal and each positive c in 4,8,...,N, invert
the nonzero integer Walsh spectrum and test whether the inverse is Boolean.
Negative c follows by reversing g. No files are written.
"""

import argparse
import json

import numpy as np


def classify(n):
    if n not in (4, 8, 16):
        raise ValueError("Only orders 4,8,16 are deliberately authorized")
    w = np.array([[1 if bin(i & j).count("1") % 2 == 0 else -1
                   for j in range(n)] for i in range(n)], dtype=np.int64)
    counts = {c: 0 for c in range(4, n + 1, 4)}
    examples = {}
    nonzero_spectra = 0
    for base in range(0, 2 ** (n - 1), 2048):
        masks = np.arange(base, min(base + 2048, 2 ** (n - 1)), dtype=np.int64)
        f = np.ones((len(masks), n), dtype=np.int64)
        f[:, 1:] = 1 - 2 * ((masks[:, None] >> np.arange(n - 1)) & 1)
        spectra = f @ w
        good = np.all(spectra != 0, axis=1)
        f, spectra = f[good], spectra[good]
        nonzero_spectra += len(f)
        for c in counts:
            divides = np.all(c % spectra == 0, axis=1)
            cand_f, cand_spectra = f[divides], spectra[divides]
            gspectra = c // cand_spectra
            numerator_g = gspectra @ w
            flat = np.all(np.abs(numerator_g) == n, axis=1)
            counts[c] += int(np.sum(flat))
            if np.any(flat) and c not in examples:
                j = int(np.flatnonzero(flat)[0])
                ff = cand_f[j]
                gg = numerator_g[j] // n
                assert np.all((ff @ w) * (gg @ w) == c)
                assert sum(int(ff[x]) * int(gg[x]) for x in range(n)) == c
                examples[c] = {"f": ff.tolist(), "g": gg.tolist()}
    return {"order": n, "f_representatives": 2 ** (n - 1),
            "nonzero_spectra": nonzero_spectra,
            "perfect_pair_counts_by_positive_c": counts,
            "examples": examples, "all_checks_integer": True}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, choices=[4, 8, 16], nargs="+",
                        default=[4, 8, 16])
    args = parser.parse_args()
    print(json.dumps([classify(n) for n in args.order], sort_keys=True))
