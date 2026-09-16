#!/usr/bin/env python3
"""Finite checks for the Wave 21 shared-prior channel calculations."""

import itertools
import math


def h2(p):
    if p in (0.0, 1.0):
        return 0.0
    return -p * math.log(p) - (1.0 - p) * math.log(1.0 - p)


def projective_noise_distribution(m, p):
    probs = {}
    for bits in itertools.product((0, 1), repeat=m):
        comp = tuple(1 - b for b in bits)
        key = min(bits, comp)
        k = sum(bits)
        probs[key] = probs.get(key, 0.0) + p**k * (1.0 - p) ** (m - k)
    return probs


def projective_entropy(m, p):
    probs = projective_noise_distribution(m, p)
    assert len(probs) == 2 ** (m - 1)
    assert abs(sum(probs.values()) - 1.0) < 1e-12
    return -sum(z * math.log(z) for z in probs.values() if z)


def projective_entropy_formula(m, p):
    ans = 0.0
    for k in range(m + 1):
        z = p**k * (1.0 - p) ** (m - k) + p ** (m - k) * (1.0 - p) ** k
        if z:
            ans -= 0.5 * math.comb(m, k) * z * math.log(z)
    return ans


def oriented_projective_states(m):
    for sigma in (-1, 1):
        for tail in itertools.product((-1, 1), repeat=m - 1):
            yield sigma, (1,) + tail


def ordered_score(A, sigma, x):
    return sum(A[i][j] * sigma * x[i] * x[j]
               for i in range(len(A)) for j in range(len(A)) if i != j)


def noisy_ground_check():
    A = [
        [0, 1, 1, -1],
        [1, 0, -1, 1],
        [1, -1, 0, 1],
        [-1, 1, 1, 0],
    ]
    states = list(oriented_projective_states(4))
    vals = [(ordered_score(A, s, x), s, x) for s, x in states]
    Q, sigma0, x0 = max(vals)
    for p, q in ((0.13, 0.0), (0.23, 0.17), (0.5, 0.5)):
        mean = 0.0
        for flips in itertools.product((-1, 1), repeat=4):
            k = sum(z == -1 for z in flips)
            pz = p**k * (1 - p) ** (4 - k)
            y = tuple(x0[i] * flips[i] for i in range(4))
            for orient_flip, pr in ((1, 1 - q), (-1, q)):
                mean += pz * pr * ordered_score(A, sigma0 * orient_flip, y)
        retained = (1 - 2 * q) * (1 - 2 * p) ** 2
        assert abs(mean - retained * Q) < 1e-10


def gibbs_check():
    A = [
        [0, 1, 1, -1],
        [1, 0, -1, 1],
        [1, -1, 0, 1],
        [-1, 1, 1, 0],
    ]
    scores = [ordered_score(A, s, x) for s, x in oriented_projective_states(4)]
    Q = max(scores)
    for beta in (0.0, 0.2, 1.3):
        weights = [math.exp(beta * z) for z in scores]
        Z = sum(weights)
        probs = [z / Z for z in weights]
        mean = sum(p * z for p, z in zip(probs, scores))
        entropy = -sum(p * math.log(p) for p in probs)
        delta = Q - mean
        deficit_Z = sum(math.exp(-beta * (Q - z)) for z in scores)
        assert abs(entropy - (beta * delta + math.log(deficit_Z))) < 1e-11
        assert delta <= entropy / beta + 1e-10 if beta else True


def f(a):
    return math.log(2.0) - h2((1.0 - math.sqrt(a)) / 2.0)


def main():
    for m in range(2, 11):
        for p in (0.0, 0.03, 0.17, 0.37, 0.5):
            H = projective_entropy(m, p)
            assert abs(H - projective_entropy_formula(m, p)) < 1e-11
            assert (m - 1) * h2(p) <= H + 1e-11
            assert H <= m * h2(p) + 1e-11
    noisy_ground_check()
    gibbs_check()

    # f(a)/a is increasing; the finite grid checks the analytic power-series proof.
    ratios = [f(k / 10000) / (k / 10000) for k in range(1, 10001)]
    assert all(x <= y + 1e-12 for x, y in zip(ratios, ratios[1:]))

    gamma = 0.672986728862
    print("max fixed-ratio coefficient budget", 27 / 256, "at rho=9/16")
    print("rho b_rho a_worst f(a_worst) a_from_known_bounds f(a_bound)")
    for rho in (0.5, 9 / 16, 0.75, 0.9):
        b = rho**1.5 - rho**2
        a_worst = 1 - b
        a_bound = max(0.0, 1 - (1 - math.sqrt(rho)) / gamma)
        print(f"{rho:.6f} {b:.12f} {a_worst:.12f} {f(a_worst):.12f} "
              f"{a_bound:.12f} {f(a_bound):.12f}")
    print("all entropy, distortion, and Gibbs identities passed")


if __name__ == "__main__":
    main()
