#!/usr/bin/env python3
"""Exact A9 audit for Wave 23 conditional-free-energy identities."""

import itertools
import math
import sys

sys.path.insert(0, "/home/math/quadra/tmp")
import parent_gibbs_prior_r22 as old


A = old.A9
PI = old.PI
N = len(A)


def lse(xs):
    return old.logsumexp(xs)


def key(d, keep):
    return old.restriction_key(d, keep)


def profile(states, keep):
    groups = {}
    score = {}
    for k, (_, d) in enumerate(states):
        y = key(d, keep)
        c = 2 * sum(
            A[i][j] * d[i][j]
            for a, i in enumerate(keep)
            for j in keep[a + 1 :]
        )
        groups.setdefault(y, []).append(k)
        score[y] = c
    assert len(groups) == 2 ** len(keep)
    assert all(len(v) == 2 ** (N - len(keep)) for v in groups.values())
    return groups, score


def full_key(d):
    return tuple(d[i][j] for i in range(N) for j in range(i + 1, N))


def flipped_key(d, keep, vertex):
    return tuple(
        -d[i][j] if vertex in (i, j) else d[i][j]
        for a, i in enumerate(keep)
        for j in keep[a + 1 :]
    )


def flipped_full_key(d, vertex):
    return tuple(
        -d[i][j] if vertex in (i, j) else d[i][j]
        for i in range(N)
        for j in range(i + 1, N)
    )


def data(states, keep, beta, gamma=None):
    if gamma is None:
        gamma = beta
    groups, score = profile(states, keep)
    logza = lse([beta * e for e, _ in states])
    logzs = lse([gamma * c for c in score.values()])
    mu = {y: math.exp(gamma * c - logzs) for y, c in score.items()}
    kval = {
        y: sum(math.exp(beta * (states[k][0] - score[y])) for k in inds)
        for y, inds in groups.items()
    }
    marginal = {
        y: math.exp(beta * score[y] + math.log(kval[y]) - logza)
        for y in groups
    }
    return groups, score, mu, kval, marginal, logza, logzs


def parent_probabilities(states, beta):
    logz = lse([beta * e for e, _ in states])
    return [math.exp(beta * e - logz) for e, _ in states]


def row_fields(d):
    return [sum(A[i][j] * d[i][j] for j in range(N) if j != i) for i in range(N)]


def audit_child_flip(states):
    index = {full_key(d): k for k, (_, d) in enumerate(states)}
    for beta in (0.1, 0.5, 1.0):
        nu = parent_probabilities(states, beta)
        mean_energy = sum(p * e for p, (e, _) in zip(nu, states))
        for omitted in range(N):
            keep = [i for i in range(N) if i != omitted]
            groups, score, _, kval, marginal, _, _ = data(states, keep, beta)
            fval = {y: math.log(z) for y, z in kval.items()}
            outside = {omitted}
            for y, inds in groups.items():
                d0 = states[inds[0]][1]
                for i in keep:
                    yi = flipped_key(d0, keep, i)
                    probs = [
                        math.exp(beta * (states[k][0] - score[y]) - fval[y])
                        for k in inds
                    ]
                    lhs = fval[yi] - fval[y]
                    rhs = math.log(
                        sum(
                            p
                            * math.exp(
                                -4
                                * beta
                                * sum(A[i][v] * states[k][1][i][v] for v in outside)
                            )
                            for p, k in zip(probs, inds)
                        )
                    )
                    assert abs(lhs - rhs) < 2e-11

            # Full-to-marginal KL chain for every retained coordinate.
            for i in keep:
                marginal_kl = 0.0
                conditional_kl = 0.0
                for y, inds in groups.items():
                    d0 = states[inds[0]][1]
                    yi = flipped_key(d0, keep, i)
                    marginal_kl += marginal[y] * math.log(marginal[y] / marginal[yi])
                    py = {
                        k: math.exp(beta * states[k][0])
                        / sum(math.exp(beta * states[j][0]) for j in inds)
                        for k in inds
                    }
                    inds_i = groups[yi]
                    denom_i = sum(math.exp(beta * states[j][0]) for j in inds_i)
                    local = 0.0
                    for k in inds:
                        ki = index[flipped_full_key(states[k][1], i)]
                        assert ki in inds_i
                        qi = math.exp(beta * states[ki][0]) / denom_i
                        local += py[k] * math.log(py[k] / qi)
                    conditional_kl += marginal[y] * local
                mean_hi = sum(p * row_fields(d)[i] for p, (_, d) in zip(nu, states))
                full_kl = 4 * beta * mean_hi
                assert abs(full_kl - marginal_kl - conditional_kl) < 3e-10
                assert conditional_kl >= -1e-12 and marginal_kl >= -1e-12
        # Coordinate-averaged version of the chain has budget 4 beta E/n.
        mean_h = [
            sum(p * row_fields(d)[i] for p, (_, d) in zip(nu, states))
            for i in range(N)
        ]
        assert abs(sum(mean_h) - mean_energy) < 1e-10
        assert min(mean_h) >= -1e-12


def audit_fixed_slice_identity(states):
    # Exact finite-population identity for adjacent child-energy differences.
    for _, d in states:
        h = row_fields(d)
        for u in range(N):
            for v in range(N):
                if u == v:
                    continue
                rest = [i for i in range(N) if i not in (u, v)]
                z = [A[u][i] * d[u][i] - A[v][i] * d[v][i] for i in rest]
                assert sum(z) == h[u] - h[v]
                pop = len(rest)
                for k in (1, 3, pop):
                    direct = 0.0
                    count = 0
                    for C in itertools.combinations(rest, k):
                        dc = 2 * sum(
                            A[u][i] * d[u][i] - A[v][i] * d[v][i] for i in C
                        )
                        direct += dc * dc
                        count += 1
                    direct /= count
                    formula = 4 * (
                        k * (pop - k) / (pop * (pop - 1)) * sum(a * a for a in z)
                        + k
                        * (k - 1)
                        / (pop * (pop - 1))
                        * (h[u] - h[v]) ** 2
                    )
                    assert abs(direct - formula) < 1e-10


def selector_statistics(states, beta):
    nu = parent_probabilities(states, beta)
    omit_data = {}
    for omitted in range(N):
        keep = [i for i in range(N) if i != omitted]
        omit_data[omitted] = (keep, data(states, keep, beta))
    core_data = {}
    for u in range(N):
        for v in range(u + 1, N):
            keep = [i for i in range(N) if i not in (u, v)]
            core_data[(u, v)] = (keep, data(states, keep, beta))

    eg2 = ef2 = ec2 = 0.0
    pairs = 0
    for ou in range(N):
        for ov in range(N):
            if ou == ov:
                continue
            # S omits ou and T omits ov.
            keep_s, ds = omit_data[ou]
            keep_t, dt = omit_data[ov]
            core_keep, dc = core_data[tuple(sorted((ou, ov)))]
            gs, cs, _, ks, ps, _, _ = ds
            gt, ct, _, kt, pt, _, _ = dt
            gc, cc, _, kc, pc, _, _ = dc
            for prob, (_, d) in zip(nu, states):
                ys, yt, yc = key(d, keep_s), key(d, keep_t), key(d, core_keep)
                Fs, Ft = math.log(ks[ys]), math.log(kt[yt])
                Gs, Gt = beta * cs[ys] + Fs, beta * ct[yt] + Ft
                # Both ratios are binary conditional probabilities over the core.
                rhs = math.log(ps[ys] / pc[yc]) - math.log(pt[yt] / pc[yc])
                assert abs((Gs - Gt) - rhs) < 3e-10
                assert abs((Fs - Ft) - ((Gs - Gt) - beta * (cs[ys] - ct[yt]))) < 1e-12
                eg2 += prob * (Gs - Gt) ** 2
                ef2 += prob * (Fs - Ft) ** 2
                ec2 += prob * (cs[ys] - ct[yt]) ** 2
            pairs += 1
    eg2 /= pairs
    ef2 /= pairs
    ec2 /= pairs
    assert eg2 <= 32 / math.e**2 + 1e-11

    # Gibbs row-field moment estimate used in the analytic bound.
    mean_energy = sum(prob * e for prob, (e, _) in zip(nu, states))
    mean_h2 = sum(
        prob * sum(z * z for z in row_fields(d)) for prob, (_, d) in zip(nu, states)
    )
    cap = (N - 1) / math.tanh(2 * beta * (N - 1))
    assert mean_h2 <= cap * mean_energy + 2e-10
    return eg2, ef2, ec2, mean_energy, mean_h2, cap


def gap_statistics(states, beta):
    avg_gap = avg_var = 0.0
    for omitted, weight in enumerate(PI):
        if not weight:
            continue
        keep = [i for i in range(N) if i != omitted]
        _, _, mu, kval, _, _, _ = data(states, keep, beta)
        f = {y: math.log(z) for y, z in kval.items()}
        ef = sum(mu[y] * f[y] for y in mu)
        avg_gap += weight * (math.log(sum(mu[y] * kval[y] for y in mu)) - ef)
        avg_var += weight * sum(mu[y] * (f[y] - ef) ** 2 for y in mu)
    return avg_gap, avg_var


def limiting_gap(states):
    parent_ground_count = sum(e == 24 for e, _ in states)
    assert parent_ground_count == 25
    answer = 0.0
    counts = []
    for omitted, weight in enumerate(PI):
        keep = [i for i in range(N) if i != omitted]
        groups, score = profile(states, keep)
        child_grounds = [y for y, c in score.items() if c == 24]
        assert all(
            [states[k][0] for k in groups[y]] == [24, 24] for y in child_grounds
        )
        qualifying = 2 * len(child_grounds)
        counts.append(qualifying)
        if weight:
            answer += weight * math.log(parent_ground_count / qualifying)
    assert counts == [4, 8, 4, 8, 4, 4, 8, 8, 8]
    return answer


def main():
    states = old.enumerate_states(A)
    audit_child_flip(states)
    audit_fixed_slice_identity(states)
    print("selector beta, E[dG^2], E[dF^2], E[dc^2], E[parent energy], E[sum h_i^2], cap")
    for beta in (0.1, 0.5, 1.0, 2.0):
        print(beta, selector_statistics(states, beta))
    print("matched beta, Jensen gap, base variance")
    for beta in (0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0):
        print(beta, *gap_statistics(states, beta))
    limit = limiting_gap(states)
    assert abs(limit - 1.5830484787467298) < 1e-13
    print("beta=infinity Jensen-gap limit", limit, "while base variance tends to zero")
    print("all Wave 23 conditional-free-energy identities passed")


if __name__ == "__main__":
    main()
