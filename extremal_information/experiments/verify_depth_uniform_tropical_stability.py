#!/usr/bin/env python3
"""Exact finite checks for depth-uniform tropical stability identities.

The script checks formulas rather than replacing their proofs.  It uses only
integer arithmetic and the standard library.
"""

from __future__ import annotations

import itertools
import random


def minplus(A, B):
    return tuple(
        tuple(min(A[i][k] + B[k][j] for k in range(len(B))) for j in range(len(B[0])))
        for i in range(len(A))
    )


def maxplus_apply(u, S):
    return tuple(max(u[a] + S[a][b] for a in range(len(u))) for b in range(len(S[0])))


def osc(v):
    return max(v) - min(v)


def hilbert(u, v):
    return osc(tuple(a - b for a, b in zip(u, v))) / 2


def line_metric(q):
    return tuple(tuple(abs(i - j) for j in range(q)) for i in range(q))


def shell(d, lam, perm):
    q = len(d)
    return tuple(tuple(lam * d[b][perm[a]] for b in range(q)) for a in range(q))


def compose_perm(left, right):
    """right o left."""
    return tuple(right[left[i]] for i in range(len(left)))


def add_gauge(D, source, target, c):
    return tuple(
        tuple(D[a][b] + source[a] - target[b] + c for b in range(len(D)))
        for a in range(len(D))
    )


def rectangular_defect(E):
    q = len(E)
    return max(
        abs(E[a][b] + E[aa][bb] - E[a][bb] - E[aa][b])
        for a, aa, b, bb in itertools.product(range(q), repeat=4)
    )


def interface_defect(E, F):
    q = len(E)
    return max(
        abs(E[a][b] - E[a][bb] + F[b][c] - F[bb][c])
        for a, b, bb, c in itertools.product(range(q), repeat=4)
    )


def projective_diameter_formula(S):
    r, s = len(S), len(S[0])
    twice = max(
        max(S[a][b] - S[a][c] for a in range(r))
        - min(S[a][b] - S[a][c] for a in range(r))
        for b, c in itertools.product(range(s), repeat=2)
    )
    return twice / 2


def row_diameter(S):
    return max(hilbert(S[a], S[aa]) for a in range(len(S)) for aa in range(len(S)))


def check_gauge_shells(rng):
    checks = 0
    for q in range(2, 7):
        d = line_metric(q)
        for _ in range(120):
            T = rng.randrange(1, 7)
            perms = []
            lambdas = []
            phis = [tuple(rng.randrange(-5, 6) for _ in range(q)) for _ in range(T + 1)]
            constants = []
            factors = []
            for t in range(T):
                # Line isometries are identity and reversal.
                p = tuple(range(q)) if rng.randrange(2) == 0 else tuple(reversed(range(q)))
                lam = rng.randrange(0, 7)
                c = rng.randrange(-4, 5)
                perms.append(p)
                lambdas.append(lam)
                constants.append(c)
                factors.append(add_gauge(shell(d, lam, p), phis[t], phis[t + 1], c))

            actual = factors[0]
            for factor in factors[1:]:
                actual = minplus(actual, factor)
            total_perm = tuple(range(q))
            for p in perms:
                total_perm = compose_perm(total_perm, p)
            expected = add_gauge(
                shell(d, min(lambdas), total_perm), phis[0], phis[-1], sum(constants)
            )
            assert actual == expected

            residuals = []
            for t in range(T):
                D = shell(d, lambdas[t], perms[t])
                residuals.append(
                    tuple(tuple(factors[t][a][b] - D[a][b] for b in range(q)) for a in range(q))
                )
                assert rectangular_defect(residuals[-1]) == 0
            for E, F in zip(residuals, residuals[1:]):
                assert interface_defect(E, F) == 0
            checks += 1
    return checks


def check_projective_diameter(rng):
    formula_checks = 0
    dominance_checks = 0
    local_isometry_checks = 0
    for r in range(2, 7):
        for s in range(2, 7):
            for _ in range(80):
                S = tuple(tuple(rng.randrange(-8, 9) for _ in range(s)) for _ in range(r))
                delta = projective_diameter_formula(S)
                assert delta == row_diameter(S)
                formula_checks += 1

                # A very large source coordinate makes its row dominate.
                images = []
                for a in range(r):
                    u = [-1000] * r
                    u[a] = 1000
                    images.append(maxplus_apply(tuple(u), S))
                    assert hilbert(images[-1], S[a]) == 0
                    dominance_checks += 1
                assert max(hilbert(x, y) for x in images for y in images) == delta

                # On a unique-winner cell, mixed selected rows preserve a
                # projective perturbation with ratio one.
                for _ in range(30):
                    u = tuple(rng.randrange(-30, 31) for _ in range(r))
                    winners = []
                    unique = True
                    margins = []
                    for b in range(s):
                        vals = sorted(((u[a] + S[a][b], a) for a in range(r)), reverse=True)
                        if vals[0][0] == vals[1][0]:
                            unique = False
                            break
                        winners.append(vals[0][1])
                        margins.append(vals[0][0] - vals[1][0])
                    if unique and len(set(winners)) >= 2:
                        i, j = winners[0], next(x for x in winners if x != winners[0])
                        h = [0] * r
                        h[i], h[j] = 1, -1
                        # Scale down if needed; integer margin >=1 allows 1/4.
                        v = tuple(u[a] + h[a] / 4 for a in range(r))
                        ratio = hilbert(maxplus_apply(u, S), maxplus_apply(v, S)) / hilbert(u, v)
                        assert abs(ratio - 1) < 1e-12
                        local_isometry_checks += 1
                        break
    return formula_checks, dominance_checks, local_isometry_checks


def check_reset_bound(rng):
    checks = 0
    for q in range(2, 7):
        for _ in range(120):
            L = rng.randrange(2, 8)
            eps = 1
            # Every Lth map is rank one plus a row-dependent term, hence an
            # exact projective reset (diameter zero).
            maps = []
            for t in range(6 * L):
                if (t + 1) % L == 0:
                    alpha = [rng.randrange(-4, 5) for _ in range(q)]
                    beta = [rng.randrange(-4, 5) for _ in range(q)]
                    S = tuple(tuple(alpha[a] + beta[b] for b in range(q)) for a in range(q))
                    assert projective_diameter_formula(S) == 0
                else:
                    S = tuple(tuple(rng.randrange(-5, 6) for _ in range(q)) for _ in range(q))
                maps.append(S)

            clean = tuple(0 for _ in range(q))
            noisy = clean
            for t, S in enumerate(maps, start=1):
                clean = maxplus_apply(clean, S)
                base = maxplus_apply(noisy, S)
                noise = tuple(rng.randrange(-1, 2) for _ in range(q))
                assert hilbert(noise, tuple(0 for _ in range(q))) <= eps
                noisy = tuple(base[i] + noise[i] for i in range(q))
                if t >= L:
                    assert hilbert(clean, noisy) <= 2 * L * eps
                    checks += 1
    return checks


def check_ising():
    checks = 0
    spins = (-1, 1)
    for J in range(-20, 21):
        S = tuple(tuple(J * s * t for t in spins) for s in spins)
        assert projective_diameter_formula(S) == 2 * abs(J)
        checks += 1
    return checks


def main():
    rng = random.Random(20260816)
    gauge = check_gauge_shells(rng)
    diameter, dominance, local = check_projective_diameter(rng)
    reset = check_reset_bound(rng)
    ising = check_ising()
    print(f"gauge-shell compositions/recognition checks: {gauge}")
    print(f"projective diameter formula checks: {diameter}")
    print(f"dominant-row realization checks: {dominance}")
    print(f"mixed-cell unit Lipschitz checks: {local}")
    print(f"depth-uniform reset checks: {reset}")
    print(f"Ising diameter checks: {ising}")


if __name__ == "__main__":
    main()
