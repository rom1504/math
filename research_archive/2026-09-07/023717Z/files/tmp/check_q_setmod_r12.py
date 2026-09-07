#!/usr/bin/env python3
"""Find small sub/supermodularity failures for S -> Q(A[S])."""

from itertools import combinations, product


def q_subset(a, mask):
    verts = [i for i in range(len(a)) if mask >> i & 1]
    if len(verts) < 2:
        return 0
    best = 0
    for tail in product((-1, 1), repeat=len(verts) - 1):
        x = (1,) + tail
        val = 0
        for ii, i in enumerate(verts):
            for jj in range(ii + 1, len(verts)):
                j = verts[jj]
                val += 2 * a[i][j] * x[ii] * x[jj]
        best = max(best, abs(val))
    return best


def main():
    for n in range(3, 7):
        edges = list(combinations(range(n), 2))
        for signs in product((-1, 1), repeat=len(edges)):
            a = [[0] * n for _ in range(n)]
            for (i, j), z in zip(edges, signs):
                a[i][j] = a[j][i] = z
            f = [q_subset(a, s) for s in range(1 << n)]
            sub = sup = None
            for s in range(1 << n):
                for t in range(1 << n):
                    defect = f[s] + f[t] - f[s | t] - f[s & t]
                    if defect < 0 and sub is None:
                        sub = (s, t, defect)
                    if defect > 0 and sup is None:
                        sup = (s, t, defect)
            if sub or sup:
                print("n", n, "signs", signs, "sub_failure", sub, "super_failure", sup)
                if sub and sup:
                    return


if __name__ == "__main__":
    main()
