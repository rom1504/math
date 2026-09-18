#!/usr/bin/env python3
"""Independent finite checks for the Wave 44 replica/common-mode memo."""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9
from high_replica_pressure_r44_explore import geometry, qnorm


def full_spins(n: int) -> np.ndarray:
    return np.asarray(list(itertools.product((-1, 1), repeat=n)), dtype=np.int64)


def logmeanexp(x: np.ndarray) -> float:
    xmax = float(np.max(x))
    return xmax + math.log(float(np.mean(np.exp(x - xmax))))


def logaddexp(a: float, b: float) -> float:
    return float(np.logaddexp(a, b))


def selector_data(a: np.ndarray, m: int):
    n = len(a)
    z = full_spins(n)
    selectors = list(itertools.combinations(range(n), m))
    q_s = np.empty(len(selectors), dtype=np.int64)
    c = np.empty((len(z), len(selectors)), dtype=np.int64)
    for j, ss in enumerate(selectors):
        s = list(ss)
        child = a[np.ix_(s, s)]
        q_s[j] = qnorm(child)
        c[:, j] = np.einsum("bi,ij,bj->b", z[:, s], child, z[:, s])
    return z, selectors, q_s, c


def consistency_check(a: np.ndarray, m: int, theta: float) -> float:
    """Check the normalized soft-label stitching identity on three views."""
    n = len(a)
    z, selectors, q_s, c = selector_data(a, m)
    row = np.einsum("bi,ij,bj->b", z, a @ a, z)
    c2 = row <= 2 * n * (n - 1)
    c2_density = float(np.mean(c2))
    js = (0, len(selectors) // 2, len(selectors) - 1)
    chosen = [selectors[j] for j in js]
    union = sorted(set().union(*map(set, chosen)))
    outside = [i for i in range(n) if i not in union]

    zeta = []
    mu_tables = []
    for j, ss in zip(js, chosen):
        y = full_spins(m)
        child = a[np.ix_(list(ss), list(ss))]
        energy = np.einsum("bi,ij,bj->b", y, child, y)
        w = np.exp(-theta * (q_s[j] - np.abs(energy)))
        zeta.append(float(np.mean(w)))
        mu_tables.append({tuple(map(int, yy)): float(ww / np.sum(w)) for yy, ww in zip(y, w)})

    # LHS under the old hard C2 law.
    product = np.ones(len(z))
    for j in js:
        product *= np.exp(-theta * (q_s[j] - np.abs(c[:, j])))
    lhs = float(np.mean(product[c2]))

    # Probability that independent child Gibbs labels stitch to one word on
    # the union and that a uniform outside extension is in C2.
    prob_consistent_low = 0.0
    for word in full_spins(len(union)):
        assignment = dict(zip(union, map(int, word)))
        mass = 1.0
        for ss, table in zip(chosen, mu_tables):
            mass *= table[tuple(assignment[i] for i in ss)]
        low_extensions = 0
        for outword in full_spins(len(outside)):
            full = np.empty(n, dtype=np.int64)
            for i in union:
                full[i] = assignment[i]
            for i, value in zip(outside, outword):
                full[i] = value
            low_extensions += int(full @ (a @ a) @ full <= 2 * n * (n - 1))
        eta = low_extensions / (2 ** len(outside))
        prob_consistent_low += mass * eta
    rhs = (
        2 ** (len(js) * m - len(union))
        * math.prod(zeta)
        * prob_consistent_low
        / c2_density
    )
    assert abs(lhs - rhs) <= 2e-12 * max(1.0, lhs)
    return lhs


def audit(name: str, a: np.ndarray, m: int) -> dict[str, float | int | str]:
    n = len(a)
    z, selectors, q_s, cmat = selector_data(a, m)
    qn = qnorm(a)
    psel = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    h = (psel**1.5 - p2) * qn
    theta = 1.0 / h
    energy = np.einsum("bi,ij,bj->b", z, a, z)
    row = np.einsum("bi,ij,bj->b", z, a @ a, z)
    c2 = row <= 2 * n * (n - 1)
    deficit = q_s[None, :] - np.abs(cmat)
    qabs = np.mean(np.exp(-theta * deficit), axis=1)
    exact_k = 1.0 + np.log(qabs)

    # Exact finite Holder inequality.  epsilon=1/2 is only a checker choice;
    # the proof uses epsilon=C*n^{-c}.
    epsilon = 0.5
    holder_p = 1.0 + epsilon
    holder_s = holder_p / epsilon
    holder_bounds = np.empty(len(z))
    for i in range(len(z)):
        logs = []
        centered = cmat[i] - p2 * energy[i]
        for sigma in (-1, 1):
            bound = theta * sigma * p2 * energy[i]
            bound += logmeanexp(-holder_p * theta * q_s) / holder_p
            bound += logmeanexp(holder_s * theta * sigma * centered) / holder_s
            logs.append(bound)
            actual_oriented = logmeanexp(-theta * q_s + theta * sigma * cmat[i])
            assert actual_oriented <= bound + 2e-12
        holder_bounds[i] = 1.0 + logaddexp(logs[0], logs[1])
        assert exact_k[i] <= holder_bounds[i] + 2e-12

    mean_q = float(np.mean(q_s))
    gap = mean_q - psel**1.5 * qn
    principal_pressure = logmeanexp(-theta * (q_s - mean_q))
    circular_core = -theta * gap + principal_pressure  # t=0

    # Exact power expansion for r=3, evaluated both ways.
    w = np.exp(-theta * deficit)
    lhs_power = float(np.mean(qabs[c2] ** 3))
    rhs_power = 0.0
    for j1, j2, j3 in itertools.product(range(len(selectors)), repeat=3):
        rhs_power += float(np.mean(w[c2, j1] * w[c2, j2] * w[c2, j3]))
    rhs_power /= len(selectors) ** 3
    assert abs(lhs_power - rhs_power) < 2e-12

    # Empty hard-fiber pair count.  This scopes, but does not falsify, the
    # strictly positive soft partition.
    favorable = (deficit <= h + 1e-12) & c2[:, None]
    empty_pairs = 0
    for j, k in itertools.combinations(range(len(selectors)), 2):
        empty_pairs += int(not np.any(favorable[:, j] & favorable[:, k]))
    empty_triples = 0
    for j, k, ell in itertools.combinations(range(len(selectors)), 3):
        empty_triples += int(not np.any(favorable[:, j] & favorable[:, k] & favorable[:, ell]))

    stitched = consistency_check(a, m, theta)
    return {
        "name": name,
        "n": n,
        "m": m,
        "selectors": len(selectors),
        "c2_projective": int(np.sum(c2) // 2),
        "projective_total": len(z) // 2,
        "max_row": int(np.max(row)),
        "max_exact_K": float(np.max(exact_k[c2])),
        "max_holder_K": float(np.max(holder_bounds[c2])),
        "G": gap,
        "Pi_Q": principal_pressure,
        "circular_core": circular_core,
        "power_r3": lhs_power,
        "stitched_tuple": stitched,
        "empty_pairs": empty_pairs,
        "all_pairs": math.comb(len(selectors), 2),
        "empty_triples": empty_triples,
        "all_triples": math.comb(len(selectors), 3),
    }


def main() -> None:
    rows = [audit("A6", A6, 5), audit("A8", A8, 6), audit("A9", A9, 7)]
    for x in rows:
        print(
            "{name}: n={n} m={m} selectors={selectors} C2={c2_projective}/{projective_total} "
            "maxR={max_row} maxK={max_exact_K:.9f} HolderK={max_holder_K:.9f} "
            "G={G:.9f} PiQ={Pi_Q:.9f} core={circular_core:.9f} "
            "E_C2[q^3]={power_r3:.9f} stitched={stitched_tuple:.9f} "
            "empty-hard-pairs={empty_pairs}/{all_pairs} "
            "empty-hard-triples={empty_triples}/{all_triples}".format(**x)
        )
    assert rows[1]["empty_pairs"] == 22
    assert rows[2]["empty_pairs"] == 362
    assert rows[1]["empty_triples"] == 1700
    assert rows[2]["empty_triples"] == 6070
    print("PASS: replica expansion, C2 stitching normalization, and Holder common-mode bound")


if __name__ == "__main__":
    main()
