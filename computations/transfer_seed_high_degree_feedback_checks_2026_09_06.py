#!/usr/bin/env python3
"""Finite algebra checks for the high-degree ACTUAL feedback probe.

No asymptotic Boolean comparison is certified by this script. It checks
normalizations, the signed tensor-rank inequality, the graph dichotomy,
and the exact finite feasible spin-ascent identity used after comparison.
Run from repository root with .venv/bin/python. No inputs or output files.
"""

from itertools import product
import math

import numpy as np
from scipy.special import eval_hermitenorm, ndtr


def rank_checks(rng):
    count = 0
    for n in range(2, 13):
        for dim in range(1, min(5, n) + 1):
            v = rng.normal(size=(n, dim)) / math.sqrt(dim)
            t = v @ v.T
            q = np.diag(t)
            d = rng.normal(size=n)  # BOTH SIGNS deliberately retained.
            for k in (1, 3, 5, 7):
                lhs = float(np.sum(d[:, None] * d[None, :] * t ** (k + 1)))
                rhs = float(np.sum(d * q ** ((k + 1) // 2)) ** 2 / n)
                assert lhs + 1e-8 * max(1, abs(lhs)) >= rhs
                count += 1
    return count


def gaussian_mixed_checks(rng):
    # A bounded smooth response permits exponentially accurate quadrature.
    # Independent literal coherent shift has a non-Gaussian two-atom law.
    nodes, weights = np.polynomial.hermite_e.hermegauss(96)
    weights /= math.sqrt(2 * math.pi)
    count = 0
    for _ in range(30):
        z = rng.normal(size=(6, 4)) / 4
        t = z @ z.T
        extra = rng.normal(size=(6, 3)) / 5
        total = t + extra @ extra.T
        shift = rng.uniform(0, 1.5, size=6)
        mask = rng.uniform(0.1, 1, size=6)
        for k in (1, 3, 5):
            for i in range(6):
                sigma2 = total[i, i]
                # psi(x)=2Phi(x)-1. Gaussian noise smooths variance to 1+sigma2.
                scale = math.sqrt(1 + sigma2)
                x = shift[i] / scale
                d = (
                    2 * mask[i] * math.exp(-x*x/2) / math.sqrt(2*math.pi)
                    * eval_hermitenorm(k-1, x)
                    / (scale**k * math.sqrt(math.factorial(k)))
                )
                # For odd k, the (-1)^(k-1) derivative sign is positive.
                for j in range(6):
                    cov = t[i, j]
                    # E[U_k(Z_Pj)|Z_total,i=sigma N]
                    # = (cov/sigma)^k h_k(N), exactly with our normalization.
                    response = mask[i] * (
                        ndtr(shift[i] + math.sqrt(sigma2)*nodes)
                        + ndtr(-shift[i] + math.sqrt(sigma2)*nodes) - 1
                    )
                    conditional_u = (cov/math.sqrt(sigma2))**k * (
                        eval_hermitenorm(k, nodes)/math.sqrt(math.factorial(k))
                    )
                    actual = float(weights @ (response * conditional_u))
                    expected = d * cov**k
                    assert abs(actual-expected) < 2e-11, (k, actual, expected)
                    count += 1
    return count


def graph_checks():
    count = 0
    # Every connected bipartite graph is an edge, a star, or contains P4.
    for left in range(1, 5):
        right = 3
        for bits in range(1, 1 << (left * right)):
            adjacency = [[] for _ in range(left + right)]
            for a in range(left):
                for b in range(right):
                    if bits >> (a*right+b) & 1:
                        adjacency[a].append(left+b)
                        adjacency[left+b].append(a)
            if any(not x for x in adjacency):
                continue
            unseen = set(range(left + right))
            while unseen:
                stack = [next(iter(unseen))]
                component = set()
                while stack:
                    a = stack.pop()
                    if a in component:
                        continue
                    component.add(a)
                    stack.extend(adjacency[a])
                unseen -= component
                edges = sum(len(adjacency[a]) for a in component)//2
                star = any(len(adjacency[a]) == len(component)-1 for a in component)
                if edges != 1 and not star:
                    p4 = any(
                        len({a, b, c, d}) == 4
                        for a in component for b in adjacency[a]
                        for c in adjacency[b] for d in adjacency[c]
                    )
                    assert p4
            count += 1
    return count


def spin_checks(rng):
    count = 0
    strict = 0
    for n in range(4, 13):
        spins = np.array(list(product((-1.0, 1.0), repeat=n)))
        for _ in range(8):
            a = rng.choice((-1.0, 1.0), size=(n, n))
            a = np.triu(a, 1)
            a += a.T
            b = a / math.sqrt(n-1)
            op = float(np.max(np.abs(np.linalg.eigvalsh(b))))
            fields = spins @ b
            alpha_old = rng.uniform(0.3, 1.3)
            f = np.sign(fields) * (np.abs(fields) > alpha_old)
            h = 1-f*f
            bf = f @ b
            c = h * np.where(bf >= 0, 1, -1)
            first_gain = float(np.mean(h*np.abs(bf)))
            candidates = []
            for endpoint in (f+c, -f+c):
                energy = np.einsum('ij,ij->i', endpoint, endpoint@b)/(2*n)
                for orientation in (-1, 1):
                    candidates.append((orientation*float(np.mean(energy)), endpoint, orientation))
            baseline, endpoint, orientation = max(candidates, key=lambda x: x[0])
            assert baseline >= first_gain - 1e-12
            local = orientation * endpoint * (endpoint @ b)
            bad = local < 0
            gain = float(np.mean(-local * bad))
            step = gain/op
            assert 0 <= step <= 1 + 1e-12
            improved = endpoint * (1-step*bad)
            actual = orientation * float(np.mean(np.einsum('ij,ij->i', improved, improved@b)))/(2*n)
            lower = baseline + gain*gain/(2*op)
            assert actual >= lower-1e-12
            cap = float(np.max(np.abs(np.einsum('ij,ij->i', spins, fields))))/(2*n)
            assert cap >= actual-1e-12
            strict += gain > 1e-10
            count += 1
    return count, strict


def main():
    rng = np.random.default_rng(20260906)
    print('signed_rank_cases', rank_checks(rng))
    print('gaussian_mixed_normalization_cases', gaussian_mixed_checks(rng))
    print('bipartite_graph_cases', graph_checks())
    cases, strict = spin_checks(rng)
    print('actual_spin_ascent_cases', cases, 'strict_finite_cases', strict)
    print('PASS: finite algebra only; no asymptotic mixed closure certified by computation')


if __name__ == '__main__':
    main()
