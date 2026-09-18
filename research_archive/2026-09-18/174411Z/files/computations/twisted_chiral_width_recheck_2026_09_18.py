#!/usr/bin/env python3
"""Direct independent profile-width enumeration without automorphism quotient."""

import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def run(n, class_index):
    source = ROOT / ("computations/results/m%d_minimizer_orbits.json" % n)
    data = json.loads(source.read_text())
    a = np.asarray(data["classes"][class_index]["representative_matrix"], dtype=np.int64)
    codes = np.arange(1 << (n-1), dtype=np.int64)
    spins = np.column_stack((np.ones(len(codes), dtype=np.int64),
                            1 - 2 * ((codes[:, None] >> np.arange(n-1)) & 1)))
    h = np.einsum("bi,ij,bj->b", spins, a, spins) // 2
    xor = codes[:, None] ^ codes[None, :]
    rooted = {}
    for p in itertools.permutations(range(n)):
        ap = a[np.ix_(p, p)]
        s = np.r_[1, ap[0, 1:]]
        b = ap * s[:, None] * s[None, :]
        rooted[b.astype(np.int8).tobytes()] = b
    hist = Counter()
    representative = {}
    digest = hashlib.sha256()
    for key in sorted(rooted):
        root = rooted[key]
        for s in spins:
            b = root * s[:, None] * s[None, :]
            e = h[:, None] - h[None, :] + spins @ b @ spins.T
            grouped = e[codes[:, None], xor]
            low, high = grouped.min(axis=0), grouped.max(axis=0)
            widths = high - low
            assert np.all(widths % 2 == 0)
            radius = int(widths.max() // 2)
            hist[radius] += 1
            digest.update(b.astype(np.int8).tobytes())
            digest.update(radius.to_bytes(2, "little"))
            if radius not in representative:
                t = int(np.argmax(widths))
                xp, xm = int(np.argmax(grouped[:, t])), int(np.argmin(grouped[:, t]))
                representative[radius] = {
                    "B": b.tolist(), "relative_spin_t": spins[t].tolist(),
                    "high": int(high[t]), "low": int(low[t]),
                    "high_x": spins[xp].tolist(), "high_y": spins[xp ^ t].tolist(),
                    "low_x": spins[xm].tolist(), "low_y": spins[xm ^ t].tolist()}
    return {"order": n, "class_index": class_index, "source": str(source.relative_to(ROOT)),
            "method": "all rooted permutation copies, every switch; direct integer matrix energies; no group quotient",
            "rooted_matrices": len(rooted), "distinct_B": sum(hist.values()),
            "unconstrained_profile_center_radius_histogram": dict(hist),
            "representatives": representative, "enumeration_sha256": digest.hexdigest()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=7)
    parser.add_argument("--class-index", type=int, default=2)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.n, args.class_index)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "representatives"}))


if __name__ == "__main__":
    main()
