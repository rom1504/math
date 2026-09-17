#!/usr/bin/env python3
"""Symbolic and finite audit of the bounded-cap critical-block embedding.

Finite caps are exact integer enumerations. Spectral norms are numerical
diagnostics, not proof inputs. The actual weak-spike existence argument is
proved by a net bound in the accompanying artifact, not by these small cases.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
import sympy as sp


def words(n: int) -> np.ndarray:
    return np.array(list(itertools.product((-1, 1), repeat=n)), dtype=np.int64)


def energies(a: np.ndarray, x: np.ndarray) -> np.ndarray:
    values = np.einsum("bi,ij,bj->b", x, a, x)
    assert np.all(values % 2 == 0)
    return values // 2


def symbolic_checks() -> dict:
    c, ca, slack = sp.symbols("C C_A slack", nonnegative=True)
    lam = 16 * (c + ca + 1) + slack
    threshold_slack = sp.expand(3 * lam / 16 - c - 2 * ca)
    assert threshold_slack == 2 * c + ca + 3 + 3 * slack / 16
    pinning_slack = sp.expand(lam - 4 * c)
    assert pinning_slack == 12 * c + 16 * ca + 16 + slack
    aa, bb, zz = sp.symbols("a b z", positive=True)
    identity = sp.expand(bb * (zz - aa / (2 * bb)) ** 2
                         - (bb * zz ** 2 - aa * zz + aa ** 2 / (4 * bb)))
    assert identity == 0
    kk = sp.symbols("K", positive=True)
    coefficient = sp.simplify(256 * (1 + kk) / kk ** 2)
    assert sp.simplify(coefficient - 256 * (kk + 1) / kk ** 2) == 0
    # This rational lower bound proves e^3>18 and therefore
    # log(18)+1<4, as used for the constant 32 in the spectral lemma.
    exp3_lower = sum(sp.Rational(3) ** j / sp.factorial(j) for j in range(6))
    assert exp3_lower > 18
    return {
        "negative_threshold_slack": str(threshold_slack),
        "pinning_hypothesis_slack": str(pinning_slack),
        "completed_square_identity": "exactly zero",
        "random_error_coefficient": str(coefficient),
        "rational_lower_bound_for_exp3": str(exp3_lower),
        "child_net_failure_bound_at_m1": 18 * math.exp(-4),
    }


def finite_case(n: int, m: int, rho: float, flipped_right_edge: bool,
                rng: np.random.Generator, trials: int) -> dict:
    x, y = words(n), words(m)
    a = np.triu(rng.choice((-1, 1), size=(n, n)), 1)
    a = a + a.T
    d = np.ones((m, m), dtype=np.int64) - np.eye(m, dtype=np.int64)
    if flipped_right_edge:
        d[0, 1] = d[1, 0] = -1
    ha, hd = energies(a, x), energies(d, y)
    qa, beta = int(np.abs(ha).max()), int(d.sum() // 2)
    r = d - rho * (np.ones((m, m)) - np.eye(m))
    cc = float(np.linalg.norm(r, 2) / math.sqrt(m))
    lam = rho * math.sqrt(m)
    threshold = beta - 2 * qa - rho * m / 2 - cc * m ** 1.5 / 2
    assert lam >= 4 * cc - 1e-12
    assert threshold >= lam * m ** 1.5 / 16 - 1e-12
    worst = {"lower": math.inf, "positive_upper": math.inf,
             "negative_gap_upper": math.inf, "absolute_upper": math.inf}
    for trial in range(trials):
        if trial == 0:
            b = np.ones((n, m), dtype=np.int64)
        elif trial == 1:
            b = np.fromfunction(lambda i, j: 1 - 2 * ((i + j) % 2), (n, m), dtype=int)
        else:
            b = rng.choice((-1, 1), size=(n, m))
        total = ha[:, None] + x @ b @ y.T + hd[None, :]
        positive, negative = int(total.max()), int(-total.min())
        cap = max(positive, negative)
        center = beta + int(np.abs(b.sum(axis=1)).sum())
        op2 = float(np.linalg.norm(b, 2) ** 2)
        err = 4 * cc ** 2 * m ** 1.5 / lam + 4 * n * op2 / (lam * math.sqrt(m))
        negerr = 4 * n * op2 / (lam * math.sqrt(m))
        slacks = {
            "lower": cap - (center - qa),
            "positive_upper": center + qa + err - positive,
            "negative_gap_upper": negerr - max(negative - positive, 0),
            "absolute_upper": center + qa + err + negerr - cap,
        }
        for key, value in slacks.items():
            assert value >= -1e-8, (n, m, trial, key, value)
            worst[key] = min(worst[key], float(value))
        if not flipped_right_edge:
            row_sums = b.sum(axis=1)
            exact_formula = beta + np.max(ha + x @ row_sums)
            assert cap == exact_formula
    return {"n": n, "m": m, "rho": rho, "right_edge_flipped": flipped_right_edge,
            "trials": trials, "exact_left_cap": qa, "exact_right_baseline": beta,
            "numerical_C": cc, "lambda": lam,
            "actual_negative_threshold": threshold,
            "required_negative_threshold": lam * m ** 1.5 / 16,
            "minimum_inequality_slacks": worst}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=60)
    parser.add_argument("--output", type=Path,
                        default=Path("tmp/paper_portfolio_2026_09_17/localization/bounded_cap_audit.json"))
    args = parser.parse_args()
    rng = np.random.default_rng(2026091703)
    cases = [finite_case(n, m, rho, flipped, rng, args.trials)
             for n, m, rho, flipped in [(2, 8, .9, False), (3, 12, .95, True),
                                         (4, 12, .95, False), (4, 12, .95, True)]]
    output = {"status": "PASS", "seed": 2026091703,
              "symbolic_checks": symbolic_checks(), "finite_checks": cases,
              "scope": "Finite tests audit pointwise inequalities; weak-spike asymptotics follow from the proof."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2), flush=True)


if __name__ == "__main__":
    main()
