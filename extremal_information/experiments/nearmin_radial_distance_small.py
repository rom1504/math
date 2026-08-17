#!/usr/bin/env python3
"""Exact small-order audit of Hamming distance to the minimizer set.

For n <= 7 this enumerates every root-gauged signing, expands the exact
minimizers by every switching, and performs a multi-source BFS in the full
edge Hamming cube.  It tests the proposed radial error bound against the
actual distance to the full (switching-closed) minimizer set.
"""

from collections import defaultdict, deque
import json
from pathlib import Path


M = {3: 3, 4: 4, 5: 4, 6: 5, 7: 9}


def edges(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def cut_masks(n, es):
    ans = []
    for spin_bits in range(1 << (n - 1)):
        # x_0=+1; bit i-1 records x_i=-1.
        mask = 0
        for k, (i, j) in enumerate(es):
            bi = 0 if i == 0 else (spin_bits >> (i - 1)) & 1
            bj = 0 if j == 0 else (spin_bits >> (j - 1)) & 1
            if bi ^ bj:
                mask |= 1 << k
        ans.append(mask)
    return ans


def cap(mask, N, cuts):
    # Negative edge bits.  H=N-2*d(mask, cut); absolute value includes the
    # global-negative augmented codeword automatically.
    return max(abs(N - 2 * bin(mask ^ c).count("1")) for c in cuts)


def audit(n):
    es = edges(n)
    N = len(es)
    cuts = cut_masks(n, es)
    root_edge_bits = [k for k, (i, _j) in enumerate(es) if i != 0]
    root_masks = []
    cap_by_root = {}
    minimizer_roots = []
    for word in range(1 << len(root_edge_bits)):
        mask = 0
        for j, k in enumerate(root_edge_bits):
            if (word >> j) & 1:
                mask |= 1 << k
        q = cap(mask, N, cuts)
        root_masks.append(mask)
        cap_by_root[mask] = q
        if q == M[n]:
            minimizer_roots.append(mask)

    total = 1 << N
    dist = bytearray([255]) * total
    queue = deque()
    for a in minimizer_roots:
        for c in cuts:
            b = a ^ c
            if dist[b] == 255:
                dist[b] = 0
                queue.append(b)
    while queue:
        a = queue.popleft()
        d = dist[a] + 1
        for k in range(N):
            b = a ^ (1 << k)
            if dist[b] == 255:
                dist[b] = d
                queue.append(b)

    by_excess = defaultdict(list)
    witnesses = {}
    for a in root_masks:
        excess = cap_by_root[a] - M[n]
        by_excess[excess].append(dist[a])
        if excess not in witnesses or dist[a] > witnesses[excess][0]:
            witnesses[excess] = (dist[a], a)
    return {
        "n": n,
        "N": N,
        "root_gauged_signings": len(root_masks),
        "root_gauged_minimizers": len(minimizer_roots),
        "switching_closed_minimizers": sum(d == 0 for d in dist),
        "levels": {
            str(e): {
                "count": len(ds),
                "min_distance": min(ds),
                "max_distance": max(ds),
                "mean_distance": sum(ds) / len(ds),
                "max_ratio_distance_over_half_excess": (
                    None if e == 0 else 2 * max(ds) / e
                ),
                "witness_negative_edge_mask": witnesses[e][1],
            }
            for e, ds in sorted(by_excess.items())
        },
    }


def main():
    results = [audit(n) for n in range(3, 8)]
    out = Path(__file__).with_name("nearmin_radial_distance_small.json")
    out.write_text(json.dumps(results, indent=2) + "\n")
    for row in results:
        first = row["levels"].get("2")
        print(row["n"], row["root_gauged_minimizers"], first)


if __name__ == "__main__":
    main()
