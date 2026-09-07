#!/usr/bin/env python3
"""Finite checks for the Wave 35 row-penalized rare-center objective.

The favorable fibers used here are the exact projective child grounds.  This
is a finite diagnostic, not the asymptotic favorable fiber with its full
allowance.  All cube states, selectors, child grounds, and flip blocks are
enumerated.  Exponentials are evaluated in double precision.
"""

from __future__ import annotations

import itertools
import math
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9, projective_spins  # noqa: E402


def canonical(x: np.ndarray) -> tuple[int, ...]:
    if x[0] < 0:
        x = -x
    return tuple(int(v) for v in x)


def projective_distance(x: np.ndarray, y: np.ndarray) -> int:
    h = int(np.count_nonzero(x != y))
    return min(h, len(x) - h)


def exact_fibers(a: np.ndarray, m: int):
    selectors = list(itertools.combinations(range(len(a)), m))
    fibers = []
    for selector in selectors:
        child = a[np.ix_(selector, selector)]
        words = list(projective_spins(m))
        energies = [abs(int(y @ child @ y)) for y in words]
        optimum = max(energies)
        fibers.append([y for y, e in zip(words, energies) if e == optimum])
    return selectors, fibers


def nearest_block(z: np.ndarray, selector, fiber) -> tuple[int, ...]:
    zs = z[list(selector)]
    candidates = []
    for y in fiber:
        for oriented in (y, -y):
            mismatch = tuple(
                selector[j] for j in range(len(selector)) if zs[j] != oriented[j]
            )
            candidates.append(mismatch)
    return min(candidates, key=lambda b: (len(b), b))


def build(a: np.ndarray, m: int, lam: float):
    n = len(a)
    selectors, fibers = exact_fibers(a, m)
    states = {}
    for z0 in projective_spins(n):
        z = np.array(z0, dtype=np.int64)
        distances = []
        blocks = []
        for selector, fiber in zip(selectors, fibers):
            block = nearest_block(z, selector, fiber)
            blocks.append(block)
            distances.append(len(block))
        q = sum(math.exp(-lam * d) for d in distances) / len(selectors)
        row = int((a @ z) @ (a @ z))
        states[canonical(z.copy())] = {
            "z": z,
            "distances": distances,
            "blocks": blocks,
            "q": q,
            "F": -math.log(q),
            "row": row,
        }
    return selectors, fibers, states


def flipped_key(z: np.ndarray, block) -> tuple[int, ...]:
    out = z.copy()
    out[list(block)] *= -1
    return canonical(out)


def audit_minimizer(
    a: np.ndarray,
    selectors,
    fibers,
    states,
    lam: float,
    gamma: float,
):
    n = len(a)
    m = len(selectors[0])
    value = min(rec["F"] + gamma * rec["row"] for rec in states.values())
    minimizers = [
        rec
        for rec in states.values()
        if abs(rec["F"] + gamma * rec["row"] - value) < 1e-11
    ]
    rec = minimizers[0]
    z = rec["z"]
    q = rec["q"]

    # Exact one-bit free-energy and row-square formulas, and Euler inequality.
    tilted = np.array(
        [math.exp(-lam * d) for d in rec["distances"]], dtype=float
    )
    tilted /= tilted.sum()
    beta_sum = 0.0
    mean_distance = float(np.dot(tilted, rec["distances"]))
    log_mgf_sum = 0.0
    delta_row_sum = 0
    for i in range(n):
        other = states[flipped_key(z, (i,))]
        deltas = np.array(other["distances"], dtype=int) - np.array(
            rec["distances"], dtype=int
        )
        assert set(deltas).issubset({-1, 0, 1})
        mgf = float(np.dot(tilted, np.exp(-lam * deltas)))
        assert abs(mgf - other["q"] / q) < 2e-12
        delta_row = other["row"] - rec["row"]
        formula = 4 * ((n - 1) - int(z[i] * (a @ a @ z)[i]))
        assert delta_row == formula
        assert math.log(mgf) <= gamma * delta_row + 2e-10
        beta_sum += float(np.dot(tilted, deltas == -1))
        log_mgf_sum += math.log(mgf)
        delta_row_sum += delta_row
    assert delta_row_sum == 4 * (n * (n - 1) - rec["row"])
    assert beta_sum + 1e-12 >= mean_distance
    assert (
        lam * (2 * mean_distance - m)
        <= log_mgf_sum + 2e-10
        <= gamma * delta_row_sum + 4e-9
    )

    # Every projective flip block satisfies the block Euler inequality.
    for size in range(n + 1):
        for block in itertools.combinations(range(n), size):
            other = states[flipped_key(z, block)]
            assert (
                math.log(other["q"] / q)
                <= gamma * (other["row"] - rec["row"]) + 2e-10
            )

    # Group selectors by a deterministic nearest mismatch block.  Flipping
    # that block lands exactly in the corresponding favorable fiber.
    groups = Counter(rec["blocks"])
    for block, count in groups.items():
        other = states[flipped_key(z, block)]
        assert count / len(selectors) <= q * math.exp(
            gamma * (other["row"] - rec["row"])
        ) + 2e-12
    partition_bound = 1.0 / sum(
        math.exp(
            gamma * (states[flipped_key(z, block)]["row"] - rec["row"])
        )
        for block in groups
    )
    assert q + 2e-12 >= partition_bound

    return {
        "number_minimizers": len(minimizers),
        "row": rec["row"],
        "F": rec["F"],
        "q": q,
        "hist": tuple(sorted(Counter(rec["distances"]).items())),
        "groups": len(groups),
        "partition_bound": partition_bound,
        "tilted_mean_distance": mean_distance,
        "beta_sum": beta_sum,
    }


def run_case(name: str, a: np.ndarray, m: int, lam: float = 5.0):
    selectors, fibers, states = build(a, m, lam)
    n = len(a)
    all_rows = [rec["row"] for rec in states.values()]
    assert max(all_rows) <= 2 * n * (n - 1)
    gamma_hard = 1.01 * lam / (4 * (n - 1))
    results = {
        "unpenalized": audit_minimizer(
            a, selectors, fibers, states, lam, gamma=0.0
        ),
        "hard_scale": audit_minimizer(
            a, selectors, fibers, states, lam, gamma=gamma_hard
        ),
    }
    print(
        name,
        "n=", n,
        "m=", m,
        "selectors=", len(selectors),
        "row_range=", (min(all_rows), max(all_rows)),
        "C2=", 2 * n * (n - 1),
        "gamma_hard=", gamma_hard,
    )
    for label, result in results.items():
        print(" ", label, result)
    return results


def main():
    # The n-1 slice gives the clearest finite pressure on exact child fibers.
    out = {
        "A6": run_case("A6", A6, 5),
        "A8": run_case("A8", A8, 7),
        "A9": run_case("A9", A9, 8),
    }
    assert out["A6"]["unpenalized"]["row"] == 30
    assert out["A8"]["unpenalized"]["row"] == 64
    assert out["A9"]["unpenalized"]["row"] == 128
    # At the generic coefficient sufficient for the C2 conversion, the A9
    # objective moves away from its high-soft-degree center even though every
    # A9 center already lies in C2.
    assert out["A9"]["hard_scale"]["row"] < 128
    assert out["A9"]["hard_scale"]["q"] < out["A9"]["unpenalized"]["q"]
    print("PASS: exact cube identities and finite row-penalty diagnostics")


if __name__ == "__main__":
    main()
