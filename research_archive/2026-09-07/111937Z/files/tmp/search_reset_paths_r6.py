#!/usr/bin/env python3
"""Exact DP search for compatible nested paths maximizing reset costs.

All energies use doubled normalization x^T A x.  At state (R,y,tau),
choose an absolute endpoint x with tau_new H_R(x)=Q(R), peel its
disagreement D={x!=y}, retain E={x=y}, and carry tau_new.  Ties are
preserved: tau_new != tau is allowed only when the opposite endpoint
strictly dominates.  At a forced reset, h=tau*H(y), q=Q=N(tau*A),
the augmented gap is g=q-h, while the reset cost in (10.390) is
a=q+h=2q-g.
"""

from functools import lru_cache
from itertools import combinations, product
import random


def energy(A, inds, bits):
    return 2 * sum(A[i][j] * bits[k] * bits[l]
                   for k, i in enumerate(inds)
                   for l, j in enumerate(inds) if k < l)


def endpoint_data(A, inds):
    vals = []
    for bits in product((-1, 1), repeat=len(inds)):
        vals.append((energy(A, inds, bits), bits))
    P = max(e for e, _ in vals) if vals else 0
    mn = min(e for e, _ in vals) if vals else 0
    N = -mn
    Q = max(P, N)
    eps = [(tau, bits) for e, bits in vals for tau in (-1, 1)
           if tau * e == Q]
    return P, N, Q, eps


def solve(A, y0, tau0=1):
    n = len(A)

    @lru_cache(None)
    def ed(inds):
        return endpoint_data(A, inds)

    @lru_cache(None)
    def dp(inds, tau):
        if not inds:
            return (0, ())
        P, N, Q, eps = ed(inds)
        y = tuple(y0[i] for i in inds)
        Hy = energy(A, inds, y)
        # Endpoint magnitudes relative to the carried orientation tau.
        same = P if tau == 1 else N
        opp = N if tau == 1 else P
        best = (0, ())
        for tau_new, x in eps:
            reset = tau_new != tau
            # Never reset on a tie, and never use the subdominant orientation.
            if reset and not (opp > same):
                continue
            if not reset and same < opp:
                continue
            augmented_gap = Q - tau * Hy
            # A zero augmented gap has no leader discrepancy to descend.
            if augmented_gap == 0:
                continue
            nxt = tuple(i for i, xb, yb in zip(inds, x, y) if xb == yb)
            # An orientation-only reset (D empty) is allowed once; after it the
            # inherited state has zero ordinary cut gap, so terminate.
            if len(nxt) == len(inds):
                tail, path = (0, ())
            else:
                tail, path = dp(nxt, tau_new)
            cost = (Q + tau * Hy) if reset else 0
            cand = (cost + tail,
                    ((inds, tau, tau_new, same, opp, Q, Hy,
                      augmented_gap, cost, x, nxt),) + path)
            if cand[0] > best[0]:
                best = cand
        return best

    P, N, Q, _ = ed(tuple(range(n)))
    val, path = dp(tuple(range(n)), tau0)
    return P, N, Q, val, path


def matrix_from_mask(n, mask):
    A = [[0] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            a = -1 if (mask >> k) & 1 else 1
            A[i][j] = A[j][i] = a
            k += 1
    return A


def show(A, y, tau=1):
    P, N, Q, val, path = solve(A, y, tau)
    print("A=", A)
    print("y=", y, "tau=", tau, "P,N,Q,S,ratios=", P, N, Q, val,
          val / Q if Q else None, val / (P + N) if P + N else None)
    for row in path:
        print(row)


def exhaustive(n):
    E = n * (n - 1) // 2
    best = (-1, None)
    for mask in range(1 << E):
        A = matrix_from_mask(n, mask)
        for y in product((-1, 1), repeat=n):
            for tau in (-1, 1):
                P, N, Q, val, path = solve(A, y, tau)
                key = val / Q if Q else 0
                if key > best[0]:
                    best = (key, (A, y, tau, P, N, Q, val, path))
    print("EXHAUSTIVE", n, "BEST", best[0])
    A, y, tau, *_ = best[1]
    show(A, y, tau)


def random_search(n, trials, seed=1):
    rng = random.Random(seed)
    best = (-1, None)
    for z in range(trials):
        E = n * (n - 1) // 2
        A = matrix_from_mask(n, rng.randrange(1 << E))
        y = tuple(rng.choice((-1, 1)) for _ in range(n))
        tau = rng.choice((-1, 1))
        P, N, Q, val, path = solve(A, y, tau)
        key = val / Q if Q else 0
        if key > best[0]:
            best = (key, (A, y, tau))
            print("RANDOM", n, z, "BEST", key, "S,Q", val, Q)
    show(*best[1])


if __name__ == "__main__":
    exhaustive(3)
    exhaustive(4)
    # n=5 full signing/y exhaustion is intentionally omitted.
    random_search(5, 3000, 7)
    random_search(6, 3000, 8)
    random_search(7, 2000, 9)
