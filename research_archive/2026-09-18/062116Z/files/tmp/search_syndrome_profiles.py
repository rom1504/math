#!/usr/bin/env python3
"""Search rooted syndrome-profile collisions invisible to outer spectra."""

from collections import defaultdict, deque
from itertools import combinations


def dist_profile(w, gens):
    size = 1 << w
    dist = [99] * size
    dist[0] = 0
    queue = deque([0])
    while queue:
        x = queue.popleft()
        for g in gens:
            y = x ^ g
            if dist[y] > dist[x] + 1:
                dist[y] = dist[x] + 1
                queue.append(y)
    return tuple(dist)


def histogram(profile):
    out = [0] * (max(profile) + 1)
    for d in profile:
        out[d] += 1
    return tuple(out)


def kernel_weight_enum(gens):
    out = [0] * (len(gens) + 1)
    for mask in range(1 << len(gens)):
        synd = 0
        for i, g in enumerate(gens):
            if (mask >> i) & 1:
                synd ^= g
        if synd == 0:
            out[bin(mask).count("1")] += 1
    return tuple(out)


def main():
    for w in range(3, 6):
        universe = list(range(1, 1 << w))
        records = []
        for k in range(w + 1, min(len(universe), w + 5) + 1):
            for gens in combinations(universe, k):
                profile = dist_profile(w, gens)
                if max(profile) >= 99:
                    continue
                records.append((gens, profile, histogram(profile)))
        buckets = defaultdict(list)
        for rec in records:
            buckets[(len(rec[0]), rec[2])].append(rec)
        for key, bucket in buckets.items():
            if len(bucket) < 2:
                continue
            for i, (s, ps, hs) in enumerate(bucket):
                for t, pt, ht in bucket[i + 1:]:
                    if ps == pt:
                        continue
                    if kernel_weight_enum(s) == kernel_weight_enum(t):
                        continue
                    # Seek a single fixed appended generator set E exposing
                    # a different covering radius after syndrome cancellation.
                    for e_size in range(1, min(4, w + 1)):
                        for e in combinations(universe, e_size):
                            rs = max(dist_profile(w, tuple(set(s) | set(e))))
                            rt = max(dist_profile(w, tuple(set(t) | set(e))))
                            if rs != rt:
                                print("w", w, "k", key[0], "hist", hs)
                                print("S", s, "profile", ps)
                                print("T", t, "profile", pt)
                                print("kernel enums", kernel_weight_enum(s), kernel_weight_enum(t))
                                print("E", e, "radii", rs, rt)
                                return
        print("no collision at w", w, "records", len(records))


if __name__ == "__main__":
    main()
