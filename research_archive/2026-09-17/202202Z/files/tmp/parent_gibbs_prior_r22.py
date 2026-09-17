#!/usr/bin/env python3
"""Exact finite verifier for the Wave 22 parent-Gibbs identities."""

import itertools
import math


A9 = [
    [0, 1, -1, 1, 1, 1, 1, -1, -1],
    [1, 0, 1, -1, -1, -1, 1, -1, -1],
    [-1, 1, 0, 1, 1, 1, 1, 1, -1],
    [1, -1, 1, 0, 1, 1, 1, -1, 1],
    [1, -1, 1, 1, 0, 1, -1, 1, -1],
    [1, -1, 1, 1, 1, 0, -1, -1, -1],
    [1, 1, 1, 1, -1, -1, 0, 1, 1],
    [-1, -1, 1, -1, 1, -1, 1, 0, -1],
    [-1, -1, -1, 1, -1, -1, 1, -1, 0],
]

PI_NUM = [4, 2, 4, 0, 4, 4, 2, 5, 0]
PI = [z / 25 for z in PI_NUM]


def logsumexp(values):
    top = max(values)
    return top + math.log(sum(math.exp(z - top) for z in values))


def enumerate_states(A):
    n = len(A)
    ans = []
    for sigma in (-1, 1):
        for tail in itertools.product((-1, 1), repeat=n - 1):
            x = (1,) + tail
            d = [[0] * n for _ in range(n)]
            for i in range(n):
                for j in range(i + 1, n):
                    d[i][j] = d[j][i] = sigma * x[i] * x[j]
            energy = sum(A[i][j] * d[i][j] for i in range(n) for j in range(n) if i != j)
            ans.append((energy, d))
    return ans


def restriction_key(d, keep):
    return tuple(d[i][j] for pos, i in enumerate(keep) for j in keep[pos + 1 :])


def child_profile(states, omitted):
    n = len(states[0][1])
    keep = [v for v in range(n) if v != omitted]
    groups = {}
    child_score = {}
    for index, (energy, d) in enumerate(states):
        key = restriction_key(d, keep)
        c = 2 * sum(A9[i][j] * d[i][j]
                    for pos, i in enumerate(keep) for j in keep[pos + 1 :])
        groups.setdefault(key, []).append(index)
        child_score[key] = c
    assert len(groups) == 2 ** len(keep)
    assert all(len(v) == 2 for v in groups.values())
    return groups, child_score


def distributions(states, profile, beta, gamma):
    groups, child_score = profile
    log_ZA = logsumexp([beta * e for e, _ in states])
    log_Zg = logsumexp([gamma * c for c in child_score.values()])
    mu = {y: math.exp(gamma * c - log_Zg) for y, c in child_score.items()}
    K = {}
    for y, inds in groups.items():
        c = child_score[y]
        K[y] = sum(math.exp(beta * (states[k][0] - c)) for k in inds)
    parent_marginal = {
        y: math.exp(beta * child_score[y]) * K[y] / math.exp(log_ZA)
        for y in groups
    }
    assert abs(sum(mu.values()) - 1) < 1e-11
    assert abs(sum(parent_marginal.values()) - 1) < 1e-11
    return log_ZA, log_Zg, mu, K, parent_marginal


def check_two_temperature(states, profiles):
    for beta, gamma in ((0.0, 0.4), (0.2, 0.7), (0.7, 0.7), (1.3, 0.5)):
        for profile in profiles:
            log_ZA, log_Zg, mu, K, parent = distributions(states, profile, beta, gamma)
            groups, score = profile
            direct = sum(mu[y] * math.log(mu[y] / parent[y]) for y in mu)
            formula = (log_ZA - log_Zg
                       + (gamma - beta) * sum(mu[y] * score[y] for y in mu)
                       - sum(mu[y] * math.log(K[y]) for y in mu))
            L = {y: math.exp((beta - gamma) * score[y]) * K[y] for y in mu}
            jensen = math.log(sum(mu[y] * L[y] for y in mu)) - sum(
                mu[y] * math.log(L[y]) for y in mu
            )
            assert abs(direct - formula) < 1e-10
            assert abs(direct - jensen) < 1e-10
            assert direct >= -1e-12
            # Conditional outside Jensen: there are two extensions for m=8.
            assert min(K.values()) >= 2 - 1e-12


def hard_event_probabilities(states, tolerance, beta):
    log_ZA = logsumexp([beta * e for e, _ in states])
    probs = []
    for omitted in range(9):
        keep = [v for v in range(9) if v != omitted]
        child_scores = []
        for energy, d in states:
            c = 2 * sum(A9[i][j] * d[i][j]
                        for pos, i in enumerate(keep) for j in keep[pos + 1 :])
            child_scores.append(c)
        Q = max(child_scores)
        event_logs = [beta * states[k][0] for k, c in enumerate(child_scores) if Q - c <= tolerance]
        probs.append(math.exp(logsumexp(event_logs) - log_ZA))
    return probs


def event_level_counts(states, tolerance):
    levels = sorted(set(e for e, _ in states))
    totals = {e: sum(E == e for E, _ in states) for e in levels}
    all_counts = []
    for omitted in range(9):
        keep = [v for v in range(9) if v != omitted]
        values = []
        for energy, d in states:
            c = 2 * sum(A9[i][j] * d[i][j]
                        for pos, i in enumerate(keep) for j in keep[pos + 1 :])
            values.append((energy, c))
        Q = max(c for _, c in values)
        counts = [sum(E == level and Q - c <= tolerance for E, c in values) for level in levels]
        fractions = [a / totals[level] for a, level in zip(counts, levels)]
        assert all(a <= b + 1e-15 for a, b in zip(fractions, fractions[1:]))
        all_counts.append(counts)
    return levels, totals, all_counts


def common_prior_chain(states, tolerance, beta):
    weights = [math.exp(beta * e) for e, _ in states]
    Z = sum(weights)
    nu = [z / Z for z in weights]
    conditionals = []
    avg_kl = 0.0
    for omitted in range(9):
        keep = [v for v in range(9) if v != omitted]
        scores = [2 * sum(A9[i][j] * d[i][j]
                          for pos, i in enumerate(keep) for j in keep[pos + 1 :])
                  for _, d in states]
        Q = max(scores)
        event = [Q - z <= tolerance for z in scores]
        mass = sum(p for p, yes in zip(nu, event) if yes)
        P = [(p / mass if yes else 0.0) for p, yes in zip(nu, event)]
        conditionals.append(P)
        if PI[omitted]:
            avg_kl += PI[omitted] * (-math.log(mass))
    output = [sum(PI[i] * conditionals[i][k] for i in range(9)) for k in range(len(states))]
    mutual = 0.0
    for i in range(9):
        if not PI[i]:
            continue
        for p, out in zip(conditionals[i], output):
            if p:
                mutual += PI[i] * p * math.log(p / out)
    marginal_kl = sum(p * math.log(p / q) for p, q in zip(output, nu) if p)
    assert abs(avg_kl - mutual - marginal_kl) < 1e-10
    return avg_kl, mutual, marginal_kl


def main():
    states = enumerate_states(A9)
    assert len(states) == 2 ** 9
    assert max(e for e, _ in states) == 24
    assert sum(e == 24 for e, _ in states) == 25
    profiles = [child_profile(states, i) for i in range(9)]
    assert all(max(score.values()) == 24 for _, score in profiles)
    check_two_temperature(states, profiles)

    expected_ground_counts = {
        0: [4, 8, 4, 8, 4, 4, 8, 8, 8],
        4: [12, 16, 12, 18, 12, 12, 16, 16, 18],
        8: [20, 22, 20, 24, 20, 20, 22, 22, 24],
        12: [24, 24, 24, 25, 24, 24, 24, 24, 25],
    }
    for tolerance, expected in expected_ground_counts.items():
        levels, totals, counts = event_level_counts(states, tolerance)
        assert levels[-1] == 24 and [row[-1] for row in counts] == expected
        # Monotone energy-level fractions prove beta -> event probability increases.
        limiting_cost = -sum(w * math.log(k / 25) for w, k in zip(PI_NUM, expected)) / 25
        finite_cost = -sum(w * math.log(p) for w, p in zip(PI, hard_event_probabilities(states, tolerance, 1.0)))
        assert finite_cost >= limiting_cost - 1e-12
        print("tolerance", tolerance, "beta=infinity average KL", limiting_cost,
              "ground counts", expected)

    avg_kl, mutual, marginal = common_prior_chain(states, 0, 1.0)
    print("A9 beta=1 exact-event avgKL, mutual information, marginal KL",
          avg_kl, mutual, marginal)
    print("all conditional-partition, Jensen-gap, and common-prior identities passed")


if __name__ == "__main__":
    main()
