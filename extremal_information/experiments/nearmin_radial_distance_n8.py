#!/usr/bin/env python3
"""Exact distance of saved order-8 cap-12 signings to all minimizers.

The two certified minimizer orbit representatives are closed under vertex
permutation, global sign, and switching.  Distances use the full labelled
edge Hamming metric; thus no symmetry is omitted.
"""

import itertools
import json
from pathlib import Path
from collections import Counter
import numpy as np


ROOT = Path(__file__).resolve().parents[2]


def edge_index(n):
    es = [(i, j) for i in range(n) for j in range(i + 1, n)]
    return es, {e: k for k, e in enumerate(es)}


def mask_of_matrix(a):
    n = len(a)
    es, _ = edge_index(n)
    return sum((a[i][j] < 0) << k for k, (i, j) in enumerate(es))


def root_gauge_permutation(a, perm, global_sign):
    n = len(a)
    b = [[global_sign * a[perm[i]][perm[j]] for j in range(n)] for i in range(n)]
    d = [1] + [b[0][i] for i in range(1, n)]
    c = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            c[i][j] = c[j][i] = b[i][j] * d[i] * d[j]
    return mask_of_matrix(c)


def cut_masks(n):
    es, _ = edge_index(n)
    ans = []
    for word in range(1 << (n - 1)):
        c = 0
        for k, (i, j) in enumerate(es):
            bi = 0 if i == 0 else (word >> (i - 1)) & 1
            bj = 0 if j == 0 else (word >> (j - 1)) & 1
            c |= (bi ^ bj) << k
        ans.append(c)
    return ans


def main():
    orbit_data = json.loads((ROOT / "computations/results/m8_minimizer_orbits.json").read_text())
    root_min = set()
    for item in orbit_data["classes"]:
        a = item["representative_matrix"]
        for perm in itertools.permutations(range(8)):
            root_min.add(root_gauge_permutation(a, perm, 1))
            root_min.add(root_gauge_permutation(a, perm, -1))
    assert len(root_min) == orbit_data["minimizing_signing_count"] == 4200

    physical = np.array(sorted({a ^ c for a in root_min for c in cut_masks(8)}), dtype=np.uint32)
    assert len(physical) == 4200 * 128
    lut = np.array([bin(i).count("1") for i in range(1 << 16)], dtype=np.uint8)

    audit = json.loads((ROOT / "extremal_information/experiments/nearmin_blind_structural_results.json").read_text())
    unique = {}
    for source, records in audit.items():
        if not isinstance(records, list):
            continue
        for rec in records:
            obs = rec.get("observables", {}) if isinstance(rec, dict) else {}
            if obs.get("n") == 8 and obs.get("cap") == 12:
                unique.setdefault(obs["matrix_sha256"], {"matrix": rec["matrix"], "sources": []})
                unique[obs["matrix_sha256"]]["sources"].append(source)

    rows = []
    for sha, rec in sorted(unique.items()):
        m = np.uint32(mask_of_matrix(rec["matrix"]))
        xor = np.bitwise_xor(physical, m)
        ds = lut[xor & np.uint32(0xFFFF)] + lut[xor >> np.uint32(16)]
        rows.append({
            "matrix_sha256": sha,
            "sources": rec["sources"],
            "cap": 12,
            "excess": 2,
            "distance_to_switching_closed_minimizer_set": int(ds.min()),
        })

    # Certify the response-metric distance for one maximally Hamming-distant
    # witness.  If two signings differ on r edges, orthogonality of the
    # quadratic Boolean characters gives max_x |sum_F c_e x_i x_j| >= sqrt(r).
    # Hence response distance below 6 is possible only for r <= 4.  We check
    # every exact minimizer in that finite radius explicitly.
    witness_row = next(row for row in rows if row["distance_to_switching_closed_minimizer_set"] == 3)
    witness_matrix = unique[witness_row["matrix_sha256"]]["matrix"]
    witness_mask = int(mask_of_matrix(witness_matrix))
    xor = np.bitwise_xor(physical, np.uint32(witness_mask))
    hamming = lut[xor & np.uint32(0xFFFF)] + lut[xor >> np.uint32(16)]
    near_indices = np.where(hamming <= 4)[0]
    sparse_response = []
    es, _ = edge_index(8)
    for idx in near_indices:
        diff = witness_mask ^ int(physical[idx])
        vals = []
        for c in cut_masks(8):
            val = 0
            for k, _edge in enumerate(es):
                if (diff >> k) & 1:
                    coeff = -1 if (witness_mask >> k) & 1 else 1
                    character = -1 if (c >> k) & 1 else 1
                    val += coeff * character
            vals.append(abs(val))
        sparse_response.append(2 * max(vals))
    response_distance = min(sparse_response)
    assert response_distance == 6

    out = {
        "status": "exact finite computation",
        "n": 8,
        "minimizer_root_gauged_count": len(root_min),
        "minimizer_switching_closed_count": len(physical),
        "saved_unique_cap12_count": len(rows),
        "distance_distribution": dict(sorted(Counter(
            row["distance_to_switching_closed_minimizer_set"] for row in rows
        ).items())),
        "certified_response_distance_witness": {
            "matrix_sha256": witness_row["matrix_sha256"],
            "upper_triangle_encoding": "/".join(
                "".join("+" if witness_matrix[i][j] > 0 else "-"
                        for j in range(i + 1, 8))
                for i in range(7)
            ),
            "distance_to_minimizers_in_edge_Hamming_metric": 3,
            "distance_to_minimizers_in_d_square": response_distance,
            "proof_note": (
                "All exact minimizers at Hamming radius <=4 were checked. "
                "Beyond radius 4, character orthogonality gives d_square "
                ">=2*sqrt(5), and even integrality makes d_square>=6."
            ),
        },
        "records": rows,
    }
    dest = Path(__file__).with_suffix(".json")
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
