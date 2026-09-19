#!/usr/bin/env python3
"""Exact finite test: paired restrictions of the stored order-26 conference.

This is a search for NONOPTIMAL child presentations, not a claim that
conference seeds solve the asymptotic problem. All caps use a source-built
integer Gray evaluator and winners receive an independent full Walsh check.
"""
import itertools
import json
import subprocess
import time
from collections import Counter
from pathlib import Path

import numpy as np

from conference_prime_square import PrimeSquare, evaluate
from twisted_chiral_director_inverse_2026_09_18 import cap, orbit_map, walsh_values

ROOT = Path(__file__).resolve().parents[1]
SCRATCH = ROOT / "tmp/twisted_chiral_2026_09_18"
OUTPUT = ROOT / "computations/results/twisted_chiral_paired_restrictions_2026_09_18.json"


def main():
    start = time.monotonic()
    field = PrimeSquare(5)
    a = np.asarray(json.loads((ROOT / "computations/results/conference_order26_gf25.json").read_text())["conference_matrix"], dtype=np.int64)
    assert np.array_equal(a, field.conference())
    # Projective involution x -> c/x, c a quadratic nonresidue; infinity<->0.
    c = next(x for x in field.elements if field.character(x) == -1)
    index = {x: i + 1 for i, x in enumerate(field.elements)}
    p = [1, 0] + [index[field.multiply(c, field.power(x, 23))] for x in field.elements[1:]]
    s = np.array([1, -1] + [field.character(x) for x in field.elements[1:]], dtype=np.int64)
    assert all(p[p[i]] == i and p[i] != i and s[i] * s[p[i]] == -1 for i in range(26))
    assert np.array_equal(s[:, None] * a[np.ix_(p, p)] * s[None, :], -a)
    pairs = [(i, p[i]) for i in range(26) if i < p[i]]
    exe = SCRATCH / "verify"
    source = ROOT / "computations/exact_fixed_signing_gray.cpp"
    if not exe.exists() or exe.stat().st_mtime < source.stat().st_mtime:
        subprocess.run(["g++", "-O3", "-std=c++17", str(source), "-o", str(exe)], check=True)
    records = []
    histogram = Counter()
    best = 1000
    winners = []
    for removed in itertools.combinations(range(13), 3):
        remaining_pairs = [pair for k, pair in enumerate(pairs) if k not in removed]
        keep = sorted(v for pair in remaining_pairs for v in pair)
        matrix = a[np.ix_(keep, keep)]
        profile = evaluate(exe, matrix)
        histogram[profile["cap"]] += 1
        records.append({"removed_pair_indices": removed, "kept_vertices": keep, "profile": profile})
        if profile["cap"] < best:
            best = profile["cap"]
            winners = [(remaining_pairs, keep, matrix, profile)]
        elif profile["cap"] == best:
            winners.append((remaining_pairs, keep, matrix, profile))
    winner = winners[0]
    values = walsh_values(winner[2])
    assert int(np.abs(values).max()) == best
    # Examine every transversal of one best signed chiral involution.
    found = []
    checked = 0
    for pairs_kept, keep, matrix, profile in winners:
        for bits in range(1 << 9):
            first = [pair[(bits >> k) & 1] for k, pair in enumerate(pairs_kept)]
            second = [p[i] for i in first]
            signs = np.array([1] * 10 + [s[i] for i in first])
            dmat = signs[:, None] * a[np.ix_(first + second, first + second)] * signs[None, :]
            child, bridge = dmat[:10, :10], dmat[:10, 10:].copy()
            diagonal = np.diag(bridge).copy()
            np.fill_diagonal(bridge, 0)
            assert np.array_equal(dmat[10:, 10:], -child)
            assert np.array_equal(bridge, bridge.T)
            match = orbit_map(child, bridge)
            checked += 1
            if match:
                found.append({"matrix": dmat.tolist(), "child": child.tolist(),
                    "child_cap": cap(child), "parent_cap": best,
                    "p": match[0], "s": match[1], "d": diagonal.tolist(),
                    "kept_vertices": keep, "parent_order": first + second,
                    "parent_signs": signs.tolist()})
                break
        if found or time.monotonic() - start > 180:
            break
    # Separately recertify the earlier, generally nonchiral cap-42 witness.
    old = json.loads((ROOT / "computations/results/conference_order26_greedy_chain_22_to19.json").read_text())
    old_record = next(r for r in old["records"] if r["order"] == 20)
    old_a = a[np.ix_(old_record["kept_vertices"], old_record["kept_vertices"])]
    old_values = walsh_values(old_a)
    old_hist = Counter(map(int, old_values))
    result = {"status": "exact finite witnesses/search; no global optimality assertion",
        "input_source": "computations/results/conference_order26_gf25.json",
        "field_nonresidue": c, "chiral_permutation": p, "chiral_signs": s.tolist(),
        "pairs": pairs, "all_pair_deletions_complete": len(records) == 286,
        "cap_histogram": dict(sorted(histogram.items())), "minimum_cap": best,
        "records": records, "best_parent_matrix": winner[2].tolist(),
        "independent_walsh_assignments": len(values), "transversals_checked": checked,
        "twisted_presentations_found": found,
        "old_order20_witness": {"matrix": old_a.tolist(), "cap": int(np.abs(old_values).max()),
            "full_spin_count": len(old_values), "energy_histogram": dict(sorted(old_hist.items())),
            "chiral_obstruction": "histogram not invariant under energy negation" if any(old_hist[t] != old_hist[-t] for t in old_hist) else None},
        "elapsed_seconds": time.monotonic() - start}
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in ["cap_histogram", "minimum_cap", "transversals_checked", "elapsed_seconds"]}))
    print("twisted presentations", len(found), "old cap", result["old_order20_witness"]["cap"])


if __name__ == "__main__":
    main()
