#!/usr/bin/env python3
"""Independent exhaustive core-cap audit of even signed-matching twists."""

from collections import Counter
import itertools
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def projective_cube(n):
    masks = np.arange(1 << (n - 1), dtype=np.uint32)
    return np.column_stack((np.ones(len(masks), dtype=np.int16),
        1 - 2 * ((masks[:, None] >> np.arange(n - 1)) & 1))).astype(np.int16)


def perfect_matchings(items):
    if not items:
        yield []
        return
    i = items[0]
    for at in range(1, len(items)):
        j = items[at]
        remaining = items[1:at] + items[at+1:]
        for rest in perfect_matchings(remaining):
            yield [(i, j)] + rest


def signed_complex_structures(n):
    for matching in perfect_matchings(list(range(n))):
        # G and -G give the same conjugate; fix the first pair orientation.
        for tail in itertools.product((-1, 1), repeat=n//2-1):
            orientations = (1,) + tail
            g = np.zeros((n, n), dtype=np.int16)
            for (i, j), sign in zip(matching, orientations):
                g[j, i] = sign
                g[i, j] = -sign
            yield g


def invariants(a, b, d=None):
    n = len(a)
    x = projective_cube(n)
    h = np.sum((x @ a) * x, axis=1) // 2
    bridge = (x @ b) @ x.T
    diff = h[:, None] - h[None, :]
    core = int((np.abs(diff) + np.abs(bridge)).max())
    values = diff + bridge
    assert core == int(np.abs(values).max())
    masks = np.arange(len(x), dtype=np.uint32)
    profiles = values[masks[None, :], np.bitwise_xor(masks[:, None], masks[None, :])]
    low, high = profiles.min(axis=1), profiles.max(axis=1)
    radius = int((high-low).max() // 2)
    result = {"child_cap": int(np.abs(h).max()), "core_cap": core,
              "unconstrained_center_radius": radius}
    if d is not None:
        filled = values + (x * d) @ x.T
        result["full_cap"] = int(np.abs(filled).max())
    return result


def main():
    records = []
    for n in (4, 6, 8):
        source = ROOT / f"computations/results/m{n}_minimizer_orbits.json"
        data = json.loads(source.read_text())
        for item in data["classes"]:
            a = np.asarray(item["representative_matrix"], dtype=np.int16)
            seen = set()
            histogram = Counter()
            radius_histogram = Counter()
            best = None
            total_g = 0
            for g in signed_complex_structures(n):
                total_g += 1
                assert np.array_equal(g.T, -g)
                assert np.array_equal(g @ g, -np.eye(n, dtype=np.int16))
                b = g.T @ a @ g
                code = b.astype(np.int8).tobytes()
                if code in seen:
                    continue
                seen.add(code)
                stats = invariants(a, b)
                histogram[stats["core_cap"]] += 1
                radius_histogram[stats["unconstrained_center_radius"]] += 1
                if best is None or stats["core_cap"] < best["core_cap"]:
                    best = dict(stats, g=g.tolist(), bridge=b.tolist())
            record = {"label": f"m{n}_class{item['class']}", "order": n,
                      "source": str(source.relative_to(ROOT)), "matrix": a.tolist(),
                      "exhaustive_skew_matching_subfamily": True,
                      "complex_structures_mod_global_sign": total_g,
                      "distinct_bridges": len(seen), "core_cap_histogram": dict(histogram),
                      "unconstrained_center_radius_histogram": dict(radius_histogram), "best": best}
            records.append(record)
            print(record["label"], len(seen), dict(histogram), dict(radius_histogram), flush=True)
    source = ROOT / "computations/results/twisted_chiral_2026_09_18_skew_matching.json"
    verified = []
    if source.exists():
        for item in json.loads(source.read_text())["records"]:
            a = np.asarray(item["child_matrix"], dtype=np.int16)
            b = np.asarray(item["bridge_hollow_matrix"], dtype=np.int16)
            d = np.asarray(item["d"], dtype=np.int16)
            g = np.zeros_like(a)
            g[np.arange(len(a)), item["p"]] = item["s"]
            assert np.array_equal(b, g @ a @ g.T)
            if len(a) % 2 == 0:
                assert np.array_equal(g.T, -g)
            stats = invariants(a, b, d)
            assert stats["full_cap"] == item["cap"]
            verified.append(dict(stats, label=item["label"], order=len(a),
                                 skew_symmetric=bool(np.array_equal(g.T, -g))))
    result = {"records": records, "independent_best_completion_checks": verified}
    output = ROOT / "computations/results/twisted_chiral_matching_core_audit_2026_09_18.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print("BEST COMPLETION CHECKS", verified, flush=True)


if __name__ == "__main__":
    main()
