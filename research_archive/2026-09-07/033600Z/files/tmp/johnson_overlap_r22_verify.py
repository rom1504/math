#!/usr/bin/env python3
"""Exhaustive A9 audit for the Wave 22 Johnson/Gibbs identities."""

from itertools import combinations, product
from math import exp, log, sqrt


A9 = (
    (0, 1, -1, 1, 1, 1, 1, -1, -1),
    (1, 0, 1, -1, -1, -1, 1, -1, -1),
    (-1, 1, 0, 1, 1, 1, 1, 1, -1),
    (1, -1, 1, 0, 1, 1, 1, -1, 1),
    (1, -1, 1, 1, 0, 1, -1, 1, -1),
    (1, -1, 1, 1, 1, 0, -1, -1, -1),
    (1, 1, 1, 1, -1, -1, 0, 1, 1),
    (-1, -1, 1, -1, 1, -1, 1, 0, -1),
    (-1, -1, -1, 1, -1, -1, 1, -1, 0),
)


def energy(a, x, sigma, vertices):
    return 2 * sigma * sum(
        a[i][j] * x[i] * x[j] for i, j in combinations(vertices, 2)
    )


def gibbs_rows(beta):
    n = 9
    selectors = [tuple(j for j in range(n) if j != i) for i in range(n)]
    states = [((1,) + tail, sigma) for tail in product((-1, 1), repeat=8)
              for sigma in (1, -1)]
    energies = [[energy(A9, x, sigma, s) for x, sigma in states]
                for s in selectors]
    laws = []
    for row in energies:
        logits = [beta * value for value in row]
        top = max(logits)
        weights = [exp(value - top) for value in logits]
        total = sum(weights)
        laws.append([value / total for value in weights])
    return selectors, energies, laws


def audit(beta):
    selectors, energies, laws = gibbs_rows(beta)
    n = len(selectors)
    mixture = [sum(laws[s][k] for s in range(n)) / n
               for k in range(len(laws[0]))]
    information = sum(
        laws[s][k] * log(laws[s][k] / mixture[k])
        for s in range(n) for k in range(len(mixture))
    ) / n
    mean_energy = sum(
        laws[s][k] * energies[s][k]
        for s in range(n) for k in range(len(mixture))
    ) / n

    kls = []
    js = []
    hs = []
    for s, t in combinations(range(n), 2):
        kl_st = sum(p * log(p / q) for p, q in zip(laws[s], laws[t]))
        kl_ts = sum(q * log(q / p) for p, q in zip(laws[s], laws[t]))
        bc = sum(sqrt(p * q) for p, q in zip(laws[s], laws[t]))
        kls.extend((kl_st, kl_ts))
        js.append(kl_st + kl_ts)
        hs.append(1 - bc)

    mean_kl = sum(kls) / len(kls)
    mean_j = sum(js) / len(js)
    mean_h = sum(hs) / len(hs)
    expected_kl = 2 * beta * mean_energy / 8
    expected_j = 4 * beta * mean_energy / 8
    assert abs(mean_kl - expected_kl) < 2e-12
    assert abs(mean_j - expected_j) < 4e-12
    # Johnson MLSI with tau_mls <= 1; d/(2n)=8/2=4.
    assert information <= 4 * mean_j + 2e-12
    return 24 - mean_energy, information, mean_h, mean_j


def main():
    print(" beta       distortion       information       mean_H2       mean_Jeffreys")
    for beta in (0.10, 0.50, 1.00, 2.00):
        values = audit(beta)
        print(f"{beta:5.2f} " + " ".join(f"{value:16.9f}" for value in values))
    print("A9 Johnson/Gibbs identities pass")


if __name__ == "__main__":
    main()
