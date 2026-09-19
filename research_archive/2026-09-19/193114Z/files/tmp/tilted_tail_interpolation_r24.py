#!/usr/bin/env python3
"""Wave 24 audit of selector collapse and the full joint interpolation."""

import itertools
import math
import sys
from collections import defaultdict
from fractions import Fraction

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, "/home/math/quadra/tmp")
import parent_gibbs_prior_r22 as old


def local_data(states, omitted, beta, gamma):
    prof = old.child_profile(states, omitted)
    logza, logzs, mu, kval, parent = old.distributions(
        states, prof, beta, gamma
    )
    score = prof[1]
    ell = {
        y: math.log(kval[y]) + (beta - gamma) * score[y] for y in mu
    }
    jensen = math.log(sum(mu[y] * math.exp(ell[y]) for y in mu)) - sum(
        mu[y] * ell[y] for y in mu
    )
    direct = sum(mu[y] * math.log(mu[y] / parent[y]) for y in mu)
    assert abs(jensen - direct) < 2e-11
    return prof, logza, logzs, mu, kval, parent, ell, jensen


def audit_soft_selector_and_joint_interpolation(states, beta, gamma):
    rows = [local_data(states, i, beta, gamma) for i in range(9)]
    js = [r[-1] for r in rows]

    # Exact optimization over selector laws for the local Jensen costs.
    soft_norm = sum(math.exp(-j) for j in js) / 9
    pi_star = [math.exp(-j) / (9 * soft_norm) for j in js]
    objective = sum(
        p * j + p * math.log(9 * p) for p, j in zip(pi_star, js)
    )
    assert abs(objective + math.log(soft_norm)) < 2e-12

    # One global interpolation on the disjoint union {(S,y)}.
    records = []
    logzs = []
    for omitted, row in enumerate(rows):
        prof, _, logzs_i, mu, kval, parent, ell, _ = row
        logzs.append(logzs_i)
        for y, c in prof[1].items():
            records.append((omitted, y, c, math.log(kval[y]), ell[y]))

    def logzt(t):
        alpha = (1 - t) * gamma + t * beta
        return old.logsumexp(
            [-math.log(9) + alpha * c + t * f for _, _, c, f, _ in records]
        )

    def moments(t):
        lz = logzt(t)
        probs = [
            math.exp(
                -math.log(9)
                + ((1 - t) * gamma + t * beta) * c
                + t * f
                - lz
            )
            for _, _, c, f, _ in records
        ]
        mean = sum(p * ell for p, (_, _, _, _, ell) in zip(probs, records))
        var = sum(
            p * (ell - mean) ** 2
            for p, (_, _, _, _, ell) in zip(probs, records)
        )
        return mean, var

    logz0, logz1 = logzt(0), logzt(1)
    logza = rows[0][1]
    assert abs(logz1 - logza) < 2e-12
    e0, _ = moments(0)
    global_kl = logz1 - logz0 - e0
    integral = quad(lambda t: (1 - t) * moments(t)[1], 0, 1, epsabs=2e-11)[0]
    assert abs(global_kl - integral) < 2e-9

    # M_0 selector is proportional to U(S) Z_S(gamma).
    top = max(logzs)
    weights = [math.exp(z - top) for z in logzs]
    pi0 = [w / sum(weights) for w in weights]
    chain = sum(p * math.log(9 * p) for p in pi0) + sum(
        p * j for p, j in zip(pi0, js)
    )
    assert abs(chain - global_kl) < 3e-11
    return js, -math.log(soft_norm), global_kl, integral


def audit_hard_selector_collapse(states, beta, tolerance):
    weights = [math.exp(beta * e) for e, _ in states]
    z = sum(weights)
    nu = [w / z for w in weights]
    omegas = []
    for omitted in range(9):
        keep = [i for i in range(9) if i != omitted]
        values = []
        for _, d in states:
            c = 2 * sum(
                old.A9[i][j] * d[i][j]
                for a, i in enumerate(keep)
                for j in keep[a + 1 :]
            )
            values.append(c)
        q = max(values)
        omegas.append(
            sum(p for p, c in zip(nu, values) if q - c <= tolerance)
        )
    bar = sum(omegas) / 9
    pi = [o / (9 * bar) for o in omegas]
    objective = sum(
        p * (-math.log(o)) + p * math.log(9 * p)
        for p, o in zip(pi, omegas)
        if p
    )
    assert abs(objective + math.log(bar)) < 2e-12
    return omegas, bar, objective


def audit_mean_deficit_identity(states, beta):
    weights = [math.exp(beta * e) for e, _ in states]
    z = sum(weights)
    nu = [w / z for w in weights]
    mean_parent = sum(p * e for p, (e, _) in zip(nu, states))
    direct = 0.0
    mean_q = 0.0
    for omitted in range(9):
        keep = [i for i in range(9) if i != omitted]
        child = [
            2
            * sum(
                old.A9[i][j] * d[i][j]
                for a, i in enumerate(keep)
                for j in keep[a + 1 :]
            )
            for _, d in states
        ]
        q = max(child)
        mean_q += q / 9
        direct += sum(p * (q - c) for p, c in zip(nu, child)) / 9
    p2 = 8 * 7 / (9 * 8)
    assert abs(direct - (mean_q - p2 * mean_parent)) < 2e-11
    thermal_gap = 24 - mean_parent
    assert thermal_gap <= 9 * math.log(2) / beta + 1e-12
    return direct, mean_q, mean_parent


def generic_nonmonotone_example():
    n = 5
    edges = list(itertools.combinations(range(n), 2))
    signs = [1, 1, 1, 1, 1, 1, 1, 1, -1, -1]
    states = []
    for sigma in (-1, 1):
        for tail in itertools.product((-1, 1), repeat=n - 1):
            x = (1,) + tail
            d = [2 * sigma * x[i] * x[j] for i, j in edges]
            states.append(d)
    energy = [sum(a * d for a, d in zip(signs, row)) for row in states]
    subset = {0, 1, 4}
    inside = [k for k, e in enumerate(edges) if set(e) <= subset]
    child = [sum(signs[k] * row[k] for k in inside) for row in states]
    event = [c >= max(child) - 4 for c in child]
    beta = 0.5
    w = [math.exp(beta * (e - max(energy))) for e in energy]
    w = [x / sum(w) for x in w]
    prob = sum(p for p, yes in zip(w, event) if yes)
    mean_e = sum(p * e for p, e in zip(w, energy))
    covariance = sum(p * e for p, e, yes in zip(w, energy, event) if yes) - prob * mean_e

    # Rigorous sign certificate.  At beta=1/2 all weights are integer powers
    # of e because every ordered energy is even.  The covariance numerator is
    # an integer Laurent polynomial P(e).  Bound e by a rational Taylor
    # interval and evaluate P interval-wise.
    n_count = defaultdict(int)
    a_count = defaultdict(int)
    for e_value, yes in zip(energy, event):
        exponent = e_value // 2
        n_count[exponent] += 1
        if yes:
            a_count[exponent] += 1
    polynomial = defaultdict(int)
    for k, ak in a_count.items():
        for ell, nell in n_count.items():
            polynomial[k + ell] += k * ak * nell - ak * ell * nell
    shift = -min(polynomial)
    terms = {k + shift: value for k, value in polynomial.items() if value}
    order = 32
    e_lo = sum(Fraction(1, math.factorial(j)) for j in range(order + 1))
    e_hi = e_lo + Fraction(1, math.factorial(order + 1)) * Fraction(order + 2, order + 1)
    p_lo = Fraction(0)
    p_hi = Fraction(0)
    for exponent, coefficient in terms.items():
        if coefficient >= 0:
            p_lo += coefficient * e_lo**exponent
            p_hi += coefficient * e_hi**exponent
        else:
            p_lo += coefficient * e_hi**exponent
            p_hi += coefficient * e_lo**exponent
    assert p_hi < 0, (p_lo, p_hi)

    assert max(energy) == 12  # Not an exact n=5 minimizer: q_5=8.
    assert abs(prob - 0.6791552231529845) < 2e-14
    assert abs(covariance - (-0.08675337519956727)) < 2e-14
    return prob, covariance, float(p_hi)


def exhaustive_small_minimizer_monotonicity():
    """Finite audit only: no negative covariance for exact n<=6 minimizers."""
    report = []
    for n in (4, 5, 6):
        edges = list(itertools.combinations(range(n), 2))
        rows = []
        for sigma in (-1, 1):
            for tail in itertools.product((-1, 1), repeat=n - 1):
                x = (1,) + tail
                rows.append([2 * sigma * x[i] * x[j] for i, j in edges])
        cuts = np.asarray(rows, dtype=np.int16)
        qbest = 10**9
        minimizers = []
        for mask in range(1 << len(edges)):
            signs = np.asarray(
                [1 if (mask >> k) & 1 else -1 for k in range(len(edges))],
                dtype=np.int16,
            )
            q = int(np.max(cuts @ signs))
            if q < qbest:
                qbest, minimizers = q, [signs]
            elif q == qbest:
                minimizers.append(signs)
        for signs in minimizers:
            full = cuts @ signs
            for m in range(3, n):
                for subset_tuple in itertools.combinations(range(n), m):
                    subset = set(subset_tuple)
                    inside = [k for k, e in enumerate(edges) if set(e) <= subset]
                    child = cuts[:, inside] @ signs[inside]
                    qchild = int(np.max(child))
                    for tolerance in sorted(set((qchild - child).tolist()))[:-1]:
                        event = child >= qchild - tolerance
                        for beta in (0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0):
                            w = np.exp(beta * (full - np.max(full)))
                            w /= np.sum(w)
                            p = float(w @ event)
                            covariance = float(w @ (full * event) - (w @ full) * p)
                            assert covariance >= -2e-10
        report.append((n, qbest, len(minimizers)))
    return report


def ground_conditioned_row_moments(states):
    grounds = [d for e, d in states if e == 24]
    all_r2 = []
    conditioned = {}
    for d in grounds:
        h = [
            sum(old.A9[i][j] * d[i][j] for j in range(9) if j != i)
            for i in range(9)
        ]
        all_r2.append(sum(z * z for z in h))
    for tolerance in (0, 4, 8, 12):
        vals = []
        for omitted in range(9):
            keep = [i for i in range(9) if i != omitted]
            for d in grounds:
                c = 2 * sum(
                    old.A9[i][j] * d[i][j]
                    for a, i in enumerate(keep)
                    for j in keep[a + 1 :]
                )
                if 24 - c <= tolerance:
                    h = [
                        sum(old.A9[i][j] * d[i][j] for j in range(9) if j != i)
                        for i in range(9)
                    ]
                    vals.append(sum(z * z for z in h))
        conditioned[tolerance] = sum(vals) / len(vals)
    assert abs(sum(all_r2) / len(all_r2) - 104.32) < 1e-12
    assert abs(conditioned[0] - 107.42857142857143) < 1e-12
    return sum(all_r2) / len(all_r2), conditioned


def main():
    states = old.enumerate_states(old.A9)
    print("soft beta: local costs; optimized selector cost; joint full-path cost")
    for beta in (0.1, 0.5, 1.0, 2.0, 4.0):
        js, opt, joint, integral = audit_soft_selector_and_joint_interpolation(
            states, beta, beta
        )
        print(beta, [round(x, 9) for x in js], opt, joint, integral)
    print("hard beta=infinity: tolerance, global probability, collapsed cost")
    counts = {
        0: [4, 8, 4, 8, 4, 4, 8, 8, 8],
        4: [12, 16, 12, 18, 12, 12, 16, 16, 18],
        8: [20, 22, 20, 24, 20, 20, 22, 22, 24],
        12: [24, 24, 24, 25, 24, 24, 24, 24, 25],
    }
    for tolerance, row in counts.items():
        bar = sum(row) / (9 * 25)
        print(tolerance, bar, -math.log(bar))
    # Finite-beta hard chain and the mean-deficit identity.
    print("hard beta=1", audit_hard_selector_collapse(states, 1.0, 0)[1:])
    print("mean deficit beta=1", audit_mean_deficit_identity(states, 1.0))
    print("ground row moments", ground_conditioned_row_moments(states))
    print("generic nonmonotone beta example", generic_nonmonotone_example())
    print("exact small minimizer audit", exhaustive_small_minimizer_monotonicity())
    print("all Wave 24 selector-collapse and interpolation identities passed")


if __name__ == "__main__":
    main()
