#!/usr/bin/env python3
"""Exact finite audit for Wave 43 colored-codegree/code-length package."""

from __future__ import annotations

import itertools
import math
import sys
from collections import Counter, defaultdict
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from anchored_conflict_r40_check import A5
from complement_consensus_r42_check import build, complement
from envelope_block_cover_r27 import A6, A8, A9


def entropy(probabilities) -> float:
    return -sum(float(p) * math.log(float(p)) for p in probabilities if p)


def map_entropy(choice) -> float:
    counts = Counter(choice)
    return entropy(Fraction(c, len(choice)) for c in counts.values())


def conflict(choice, selectors, states, n) -> Fraction:
    count = len(selectors)
    numerator = 0
    common_degree = None
    for i in range(1, n):
        plus = sum(bool(states[d][1][i] == 1)
                   for s0, d in zip(selectors, choice) if i in s0)
        minus = sum(bool(states[d][1][i] == -1)
                    for s0, d in zip(selectors, choice) if i in s0)
        if plus + minus:
            if common_degree is None:
                common_degree = plus + minus
            assert common_degree == plus + minus
            numerator += plus * minus
    return Fraction(2 * numerator, count * common_degree)


def precedence(choice, selectors, fibers):
    support = set(choice)
    graph = {d: set() for d in support}
    for s0, selected in zip(selectors, choice):
        graph[selected] |= (support & fibers[s0]) - {selected}
    return graph


def acyclic(graph) -> bool:
    indegree = {d: 0 for d in graph}
    for d in graph:
        for e in graph[d]:
            indegree[e] += 1
    stack = [d for d in graph if indegree[d] == 0]
    seen = 0
    while stack:
        d = stack.pop()
        seen += 1
        for e in graph[d]:
            indegree[e] -= 1
            if indegree[e] == 0:
                stack.append(e)
    return seen == len(graph)


def incidence_entropy_bound(choice, selectors, fibers) -> tuple[float, float]:
    degree = defaultdict(int)
    for s0 in selectors:
        for d in fibers[s0]:
            degree[d] += 1
    lhs = map_entropy(choice)
    rhs = math.log(len(selectors)) - sum(
        math.log(degree[d]) for d in choice
    ) / len(selectors)
    assert lhs + 1e-12 >= rhs
    return lhs, rhs


def audit_matrix(a: np.ndarray, m: int, name: str, priority_seed):
    q, selectors, states, fibers = build(a, m, False)
    choices = itertools.product(*(tuple(fibers[s0]) for s0 in selectors))
    zero = []
    rationalizable_zero = []
    seed_graph = precedence(priority_seed, selectors, fibers)
    assert acyclic(seed_graph)
    best_rational_conflict = conflict(priority_seed, selectors, states, len(a))
    best_rational = []
    total_maps = 0

    for choice in choices:
        total_maps += 1
        c = conflict(choice, selectors, states, len(a))
        row = Fraction(sum(states[d][2] for d in choice), len(choice))
        if c == 0:
            h = map_entropy(choice)
            graph = precedence(choice, selectors, fibers)
            is_rationalizable = acyclic(graph)
            zero.append((choice, h, row, graph))
            if is_rationalizable:
                rationalizable_zero.append((choice, h, row, graph))
        if c <= best_rational_conflict:
            graph = precedence(choice, selectors, fibers)
            is_rationalizable = acyclic(graph)
        else:
            is_rationalizable = False
        if is_rationalizable:
            h = map_entropy(choice)
            if best_rational_conflict is None or c < best_rational_conflict:
                best_rational_conflict = c
                best_rational = [(choice, h, row)]
            elif c == best_rational_conflict:
                best_rational.append((choice, h, row))

    best_zero = min(zero, key=lambda item: (item[1], item[2]))
    choice, h, row, _ = best_zero
    incidence_entropy_bound(choice, selectors, fibers)
    incidence_entropy_bound(best_rational[0][0], selectors, fibers)
    witness = []
    for s0, d in zip(selectors, choice):
        sigma, x, r2 = states[d]
        value = sigma * int(x @ complement(a, s0) @ x)
        assert value >= q
        witness.append((s0, d, sigma, tuple(map(int, x)), r2, value))

    result = {
        "name": name,
        "maps": total_maps,
        "zero_maps": len(zero),
        "rationalizable_zero_maps": len(rationalizable_zero),
        "best_zero_choice": choice,
        "best_zero_entropy": h,
        "best_zero_row": row,
        "best_priority_conflict": best_rational_conflict,
        "best_priority_count": len(best_rational),
        "witness": witness,
    }
    print(result)
    return result, (selectors, states, fibers)


def randomized_chain_rule(selectors, maps):
    """Verify E_J I(S;D|J)=I(S;D)+I(S;J|D) on a nontrivial mixture."""
    # J has probabilities 1/3,2/3 and is independent of uniform S.
    pj = (Fraction(1, 3), Fraction(2, 3))
    ns = len(selectors)
    joint_sd = defaultdict(Fraction)
    joint_sjd = defaultdict(Fraction)
    joint_d = defaultdict(Fraction)
    joint_jd = defaultdict(Fraction)
    conditional_information = 0.0
    for j, (prob_j, choice) in enumerate(zip(pj, maps)):
        conditional_information += float(prob_j) * map_entropy(choice)
        for si, d in enumerate(choice):
            mass = prob_j / ns
            joint_sd[(si, d)] += mass
            joint_sjd[(si, j, d)] += mass
            joint_d[d] += mass
            joint_jd[(j, d)] += mass

    # I(S;D), with S uniform.
    i_sd = 0.0
    for (si, d), mass in joint_sd.items():
        i_sd += float(mass) * math.log(float(mass / (Fraction(1, ns) * joint_d[d])))

    # I(S;J|D)=sum p(s,j,d) log[p(s,j|d)/(p(s|d)p(j|d))].
    i_sj_d = 0.0
    for (si, j, d), mass in joint_sjd.items():
        p_sd = joint_sd[(si, d)]
        p_jd = joint_jd[(j, d)]
        i_sj_d += float(mass) * math.log(float(mass * joint_d[d] / (p_sd * p_jd)))
    assert abs(conditional_information - i_sd - i_sj_d) < 1e-12
    print({"conditional_I": conditional_information,
           "marginal_I": i_sd, "extra_I": i_sj_d})


def incidence_collapse_audit(a: np.ndarray, m: int, name: str):
    """Verify conditional-KL support collapse on a nonuniform joint law."""
    q, selectors, states, fibers = build(a, m, False)
    ns = len(selectors)
    active = sorted(set().union(*fibers.values()))
    columns = {d: [si for si, s0 in enumerate(selectors) if d in fibers[s0]]
               for d in active}
    hcol = {d: math.log(Fraction(ns, len(columns[d]))) for d in active}

    # Arbitrary positive P_D and deliberately nonuniform P(S|d).
    raw_d = {d: (d % 7) + 1 for d in active}
    norm_d = sum(raw_d.values())
    pd = {d: Fraction(raw_d[d], norm_d) for d in active}
    cond = {}
    joint = defaultdict(Fraction)
    ps = defaultdict(Fraction)
    for d in active:
        raw_s = {si: (si % 5) + 1 for si in columns[d]}
        norm_s = sum(raw_s.values())
        cond[d] = {si: Fraction(raw_s[si], norm_s) for si in columns[d]}
        for si, mass in cond[d].items():
            joint[(si, d)] = pd[d] * mass
            ps[si] += joint[(si, d)]

    htot_cond = sum(
        float(pd[d] * mass) * math.log(float(mass * ns))
        for d in active for mass in cond[d].values()
    )
    kl_selector = sum(float(mass) * math.log(float(mass * ns))
                      for mass in ps.values())
    mutual = sum(
        float(mass) * math.log(float(mass / (ps[si] * pd[d])))
        for (si, d), mass in joint.items()
    )
    assert abs(htot_cond - kl_selector - mutual) < 1e-11

    mean_hcol = sum(float(pd[d]) * hcol[d] for d in active)
    assert mean_hcol <= htot_cond + 1e-12
    mean_row = sum(float(pd[d]) * states[d][2] for d in active)
    selected = min(active, key=lambda d: hcol[d] / htot_cond
                   + states[d][2] / mean_row)
    assert hcol[selected] <= 2 * htot_cond + 1e-12
    assert states[selected][2] <= 2 * mean_row + 1e-12

    # Every actual incidence obeys (10.1044), including orientation.
    p2 = Fraction(m * (m - 1), len(a) * (len(a) - 1))
    assert p2 >= Fraction(1, 2)
    for si, s0 in enumerate(selectors):
        idx = np.asarray(s0)
        child = a[np.ix_(idx, idx)]
        child_cap = max(abs(int(z @ child @ z))
                        for z in itertools.product((-1, 1), repeat=m))
        assert child_cap <= q
        for d in fibers[s0]:
            sigma, x, _ = states[d]
            energy = sigma * int(x @ a @ x)
            delta = q - energy
            c_s = sigma * int(x[idx] @ child @ x[idx])
            assert 2 * c_s - energy >= q
            assert 2 * c_s >= 2 * q - delta

    # Exact scalarization: the lower bound is attained by a constant D and
    # uniform conditional law on its incidence column.
    scalar = {}
    for lam in (0.0, 0.01, 0.1, 1.0):
        dstar = min(active, key=lambda d: hcol[d] + lam * states[d][2])
        value = hcol[dstar] + lam * states[dstar][2]
        scalar[lam] = (dstar, value)
        uniform_kl = math.log(Fraction(ns, len(columns[dstar])))
        assert abs(value - (uniform_kl + lam * states[dstar][2])) < 1e-12

    result = {
        "name": name,
        "selectors": ns,
        "active_states": len(active),
        "Htot": htot_cond,
        "mean_column_cost": mean_hcol,
        "mean_row": mean_row,
        "extracted": (selected, hcol[selected], states[selected][2]),
        "scalar_minima": scalar,
    }
    print(result)
    return result


def main():
    a5, data5 = audit_matrix(np.asarray(A5), 4, "A5", (24, 24, 25, 25))
    a6, data6 = audit_matrix(np.asarray(A6), 5, "A6", (2, 2, 2, 11, 42))

    assert a5["zero_maps"] == 57
    assert a5["rationalizable_zero_maps"] == 27
    assert a5["best_zero_choice"] == (24, 24, 25, 25)
    assert abs(a5["best_zero_entropy"] - math.log(2)) < 1e-12
    assert a5["best_zero_row"] == 16
    assert a5["best_priority_conflict"] == 0

    assert a6["zero_maps"] == 22
    assert a6["rationalizable_zero_maps"] == 0
    assert a6["best_zero_choice"] == (2, 4, 8, 16, 32)
    assert abs(a6["best_zero_entropy"] - math.log(5)) < 1e-12
    assert a6["best_zero_row"] == 30
    assert a6["best_priority_conflict"] == Fraction(2, 5)
    assert a6["best_priority_count"] == 40
    assert all(record[-1] == 10 for record in a6["witness"])

    # A genuine two-cycle in the displayed A6 zero-conflict map.
    selectors, _, fibers = data6
    graph = precedence(a6["best_zero_choice"], selectors, fibers)
    assert 4 in graph[2] and 2 in graph[4]

    # Use two distinct A5 incidence maps for the shared-latent chain rule.
    rational_maps = []
    for choice in itertools.product(*(tuple(data5[2][s]) for s in data5[0])):
        if acyclic(precedence(choice, data5[0], data5[2])):
            rational_maps.append(choice)
            if len(rational_maps) == 2:
                break
    randomized_chain_rule(data5[0], rational_maps)
    for name, a, m in (("A5", np.asarray(A5), 4), ("A6", np.asarray(A6), 5),
                       ("A8", np.asarray(A8), 6), ("A9", np.asarray(A9), 7)):
        incidence_collapse_audit(a, m, name)
    print("PASS colored_codegree_r43_check")


if __name__ == "__main__":
    main()
