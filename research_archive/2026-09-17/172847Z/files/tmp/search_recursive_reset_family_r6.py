#!/usr/bin/env python3
"""Search the alternating nested-clique recurrence for forced-reset towers.

C_new = [[-(J-I), J], [J, -C_old]].  A profile records exact min/max
energy at each magnetization, so orders far beyond brute-force 2^n are exact.
"""

from math import isqrt


def base_one_negative_triangle():
    # Exact profile by direct enumeration for the triangle with edge 01 negative.
    A = [[0, -1, 1], [-1, 0, 1], [1, 1, 0]]
    prof = {}
    for mask in range(8):
        x = [1 if (mask >> i) & 1 else -1 for i in range(3)]
        M = sum(x)
        H = sum(A[i][j] * x[i] * x[j] for i in range(3) for j in range(3))
        lo, hi = prof.get(M, (H, H))
        prof[M] = (min(lo, H), max(hi, H))
    return prof


def extrema(prof):
    P = max(hi for lo, hi in prof.values())
    mn = min(lo for lo, hi in prof.values())
    n = max(prof)
    Hall = prof[n][0]
    assert prof[n][0] == prof[n][1]
    return n, P, -mn, Hall


def extend(prof, s):
    m = max(prof)
    out = {}
    for U in range(-s, s + 1, 2):
        for V, (lo, hi) in prof.items():
            W = U + V
            const = s - U * U + 2 * U * V
            # Core is -C_old, so its extrema at V are -hi and -lo.
            vals = (const - hi, const - lo)
            if W in out:
                olo, ohi = out[W]
                out[W] = (min(olo, vals[0]), max(ohi, vals[1]))
            else:
                out[W] = vals
    return out


def valid_extension(prof, s, verbose=False):
    m, P0, N0, H0 = extrema(prof)
    assert H0 == P0 and N0 > P0
    nxt = extend(prof, s)
    n, P, N, Hall = extrema(nxt)
    # all-one state is required to be the positive (subdominant) ground.
    ok = Hall == P and N > P
    # The desired negative ground flips exactly the new top block.
    Hflip = s - s*s - H0 - 2*s*m
    ok = ok and Hflip == -N
    if verbose or ok:
        print("m,s,n old(P,N)", m, s, n, (P0, N0),
              "new(P,N,H1,Hflip)", (P, N, Hall, Hflip), "OK", ok)
    return ok, nxt


def greedy(steps=20):
    prof = base_one_negative_triangle()
    print("BASE", extrema(prof))
    sizes = []
    reset_sum = 0
    for k in range(steps):
        m = max(prof)
        good = []
        for s in range(1, m + 2):
            ok, nxt = valid_extension(prof, s)
            if ok:
                n, P, N, Hall = extrema(nxt)
                # reset cost at this level with inherited all-one and old +.
                a = N + P
                good.append((s, nxt, P, N, a))
        if not good:
            print("NO EXTENSION at", m)
            break
        # Smallest added block maximizes depth.
        s, prof, P, N, a = good[0]
        sizes.append(s)
        reset_sum += a
        print("CHOOSE", s, "n", max(prof), "P,N", P, N,
              "cum reset/Q", reset_sum, reset_sum / N)
    print("sizes outerward", sizes)


if __name__ == "__main__":
    greedy(30)
