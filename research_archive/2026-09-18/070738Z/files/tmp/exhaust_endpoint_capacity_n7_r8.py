#!/usr/bin/env python3
"""Exhaust switching classes at order seven for Wave-8 endpoint inequalities."""

from functools import lru_cache
from itertools import combinations, product


def energy(A, x):
    return 2 * sum(A[i][j] * x[i] * x[j]
                   for i in range(len(A)) for j in range(i + 1, len(A)))


@lru_cache(None)
def extrema(A):
    n = len(A)
    xs = ((1,) + t for t in product((-1, 1), repeat=max(0, n - 1)))
    rows = tuple((energy(A, x), x) for x in xs)
    P = max(v for v, _ in rows) if rows else 0
    lo = min(v for v, _ in rows) if rows else 0
    return P, -lo, tuple(x for v, x in rows if v == P), tuple(x for v, x in rows if v == lo)


def induced(A, ids):
    return tuple(tuple(A[i][j] for j in ids) for i in ids)


def matrix7(code):
    A = [[0] * 7 for _ in range(7)]
    for j in range(1, 7):
        A[0][j] = A[j][0] = 1
    for k, (i, j) in enumerate(combinations(range(1, 7), 2)):
        A[i][j] = A[j][i] = 1 if code >> k & 1 else -1
    return tuple(tuple(row) for row in A)


def main():
    min_main = (10**9, None)
    min_strong = (10**9, None)
    pairs = 0
    for code in range(1 << 15):
        A = matrix7(code)
        P, N, ps, ns = extrema(A)
        R = P + N
        for p in ps:
            for n in ns:
                same = tuple(i for i in range(7) if p[i] == n[i])
                diff = tuple(i for i in range(7) if p[i] != n[i])
                if not same or not diff:
                    continue
                info = []
                for ids, other in ((same, diff), (diff, same)):
                    AX = induced(A, ids)
                    PX, NX, _, _ = extrema(AX)
                    QX = max(PX, NX)
                    px = tuple(p[i] for i in ids)
                    h = energy(AX, px)
                    L = sum(abs(sum(A[j][i] * p[i] for i in ids)) for j in other)
                    info.append((QX, h, L))
                H = sum(z[1] for z in info)
                M = sum(2 * z[2] - z[0] + abs(z[1]) for z in info)
                main_slack = M - abs(P - N)
                strong_slack = R - abs(H) - sum(z[0] for z in info)
                data = (code, P, N, p, n, info)
                min_main = min(min_main, (main_slack, data), key=lambda z: z[0])
                min_strong = min(min_strong, (strong_slack, data), key=lambda z: z[0])
                pairs += 1
                if main_slack < 0 or strong_slack < 0 or any(4*z[2] != R for z in info):
                    print("FAIL", main_slack, strong_slack, data)
                    return 1
        if code % 4096 == 4095:
            print("progress", code + 1, "pairs", pairs,
                  "min", min_main[0], min_strong[0], flush=True)
    print("PASS matrices", 1 << 15, "pairs", pairs)
    print("min_main", min_main)
    print("min_strong", min_strong)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
