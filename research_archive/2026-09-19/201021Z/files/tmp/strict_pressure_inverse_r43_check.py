#!/usr/bin/env python3
"""Finite audit for the Wave 43 absolute mixed-pressure extraction.

This checks the exact tilt/KL identities, the Markov extraction step, the
absolute-to-oriented sign split, and the finite strength of the global
center average on the certified minimizers A6, A8, and A9.  The normalized
choice theta * H = 1 is only a finite diagnostic, not an asymptotic test.
"""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9, projective_spins


def qnorm(a: np.ndarray) -> int:
    xs = np.asarray(list(projective_spins(len(a))), dtype=np.int64)
    return int(np.max(np.abs(np.einsum("bi,ij,bj->b", xs, a, xs))))


def audit(name: str, a: np.ndarray, m: int) -> dict[str, float | int | str]:
    n = len(a)
    zs = np.asarray(list(projective_spins(n)), dtype=np.int64)
    selectors = list(itertools.combinations(range(n), m))
    qn = qnorm(a)
    p = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    h = (p**1.5 - p2) * qn
    assert h > 0
    theta = 1.0 / h

    deficit = np.empty((len(zs), len(selectors)), dtype=np.int64)
    signed_energy = np.empty_like(deficit)
    child_norms = np.empty(len(selectors), dtype=np.int64)
    for j, selector in enumerate(selectors):
        s = list(selector)
        child = a[np.ix_(s, s)]
        q_s = qnorm(child)
        child_norms[j] = q_s
        c_s = np.einsum("bi,ij,bj->b", zs[:, s], child, zs[:, s])
        signed_energy[:, j] = c_s
        deficit[:, j] = q_s - np.abs(c_s)
    assert np.min(deficit) >= 0

    weights = np.exp(-theta * deficit)
    q_z = np.mean(weights, axis=1)
    partition = float(np.mean(q_z))

    # Joint tilt against U_z x U_m and its center marginal.
    joint = weights / np.sum(weights)
    joint_kl = float(np.sum(joint * np.log(joint * joint.size)))
    tilted_mean_deficit = float(np.sum(joint * deficit))
    assert abs(joint_kl - (-theta * tilted_mean_deficit - math.log(partition))) < 2e-12
    assert joint_kl <= -math.log(partition) + 2e-12

    marginal = q_z / np.sum(q_z)
    marginal_kl = float(np.sum(marginal * np.log(marginal * len(marginal))))
    assert marginal_kl <= joint_kl + 2e-12  # KL data processing.

    row = np.einsum("bi,ij,bj->b", zs, a @ a, zs)
    mean_row_tilt = float(marginal @ row)
    good = row <= 2.0 * mean_row_tilt + 1e-12
    good_mass = float(np.sum(marginal[good]))
    assert good_mass >= 0.5 - 2e-12
    extracted_q = float(np.max(q_z[good]))
    assert extracted_q >= partition / 2.0 - 2e-12

    # The best center diagnoses how much the global average can lose.  Check
    # soft-to-hard and then split the absolute event into one fixed sign.
    best = int(np.argmax(q_z))
    pressure = 1.0 + math.log(float(q_z[best]))
    absolute_event = deficit[best] <= h + 1e-12
    absolute_mass = float(np.mean(absolute_event))
    if pressure > 0:
        soft_hard = (math.exp(pressure) - 1.0) / (math.e - 1.0)
        assert absolute_mass + 2e-12 >= soft_hard
    c_best = signed_energy[best]
    plus = absolute_event & (c_best >= 0)
    minus = absolute_event & (c_best < 0)
    oriented_mass = max(float(np.mean(plus)), float(np.mean(minus)))
    assert oriented_mass + 2e-12 >= absolute_mass / 2.0

    # Every finite center in these examples happens to satisfy the old C2
    # cap.  That is finite evidence only; the proof extraction uses the
    # weaker project-scale O(n^(9/4-c)) cap.
    c2_cap = 2 * n * (n - 1)
    assert int(np.max(row)) <= c2_cap

    return {
        "name": name,
        "n": n,
        "m": m,
        "q_n": qn,
        "H": h,
        "theta": theta,
        "Z": partition,
        "minus_log_Z": -math.log(partition),
        "max_pressure": pressure,
        "best_absolute_tail": absolute_mass,
        "best_oriented_tail": oriented_mass,
        "joint_KL": joint_kl,
        "center_KL": marginal_kl,
        "tilted_mean_row": mean_row_tilt,
        "good_center_mass": good_mass,
        "max_row": int(np.max(row)),
        "C2_cap": c2_cap,
    }


def main() -> None:
    rows = [audit("A6", A6, 5), audit("A8", A8, 6), audit("A9", A9, 7)]
    for row in rows:
        print(
            "{name}: n={n} m={m} q={q_n} H={H:.9f} theta={theta:.9f} "
            "Z={Z:.9f} -logZ={minus_log_Z:.9f} maxK={max_pressure:.9f} "
            "tails(abs,or)=({best_absolute_tail:.9f},{best_oriented_tail:.9f}) "
            "KL(joint,center)=({joint_KL:.9f},{center_KL:.9f}) "
            "EtiltR={tilted_mean_row:.6f} Ptilt(G)={good_center_mass:.6f} "
            "maxR/C2={max_row}/{C2_cap}".format(**row)
        )
    assert rows[1]["minus_log_Z"] > 1.0 and rows[1]["max_pressure"] > 0.0
    assert rows[2]["minus_log_Z"] > 1.0 and rows[2]["max_pressure"] > 0.0
    print("PASS: tilt/KL/extraction/sign-split identities and finite average-pressure wall")


if __name__ == "__main__":
    main()
