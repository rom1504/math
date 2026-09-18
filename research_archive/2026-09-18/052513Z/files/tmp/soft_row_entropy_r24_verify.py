from __future__ import annotations

import math
import random
import importlib.util
import itertools
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


ROOT = Path(__file__).resolve().parents[1]


def check_reference_identity(seed: int = 2401) -> None:
    rng = np.random.default_rng(seed)
    for k in range(2, 9):
        for _ in range(30):
            p = rng.dirichlet(np.ones(k))
            lam = rng.normal(size=k)

            def objective(z: np.ndarray) -> float:
                nu = np.exp(z - np.max(z))
                nu /= nu.sum()
                return float(np.dot(p, -np.log(nu)) + math.log(np.dot(nu, np.exp(lam))))

            out = minimize(objective, np.zeros(k), method="BFGS", tol=1e-11)
            exact = float(-np.dot(p, np.log(p)) + np.dot(p, lam))
            assert abs(out.fun - exact) < 2e-8, (out.fun, exact)

            nu = p * np.exp(-lam)
            nu /= nu.sum()
            direct = float(np.dot(p, -np.log(nu)) + math.log(np.dot(nu, np.exp(lam))))
            assert abs(direct - exact) < 1e-11

    # Stochastic-channel form: inf_nu D(P || U x nu) + log E_nu exp(Lambda)
    # equals I(S;D) + E_{P_D} Lambda.
    for ns, nd in [(3, 4), (5, 3), (6, 7)]:
        joint = rng.dirichlet(np.ones(ns * nd)).reshape(ns, nd)
        psel = joint.sum(axis=1)
        pout = joint.sum(axis=0)
        lam = rng.normal(size=nd)
        mutual = float(np.sum(joint * np.log(joint / (psel[:, None] * pout[None, :]))))
        selector_kl = float(np.dot(psel, np.log(psel * ns)))
        exact = mutual + selector_kl + float(np.dot(pout, lam))

        def stochastic_objective(z: np.ndarray) -> float:
            nu = np.exp(z - np.max(z))
            nu /= nu.sum()
            divergence = float(np.sum(joint * np.log(joint / ((1 / ns) * nu[None, :]))))
            return divergence + math.log(float(np.dot(nu, np.exp(lam))))

        out = minimize(stochastic_objective, np.zeros(nd), method="BFGS", tol=1e-11)
        assert abs(out.fun - exact) < 2e-8, (out.fun, exact)


def check_hoeffding(seed: int = 2402) -> None:
    rng = random.Random(seed)
    for n in range(2, 11):
        for _ in range(100):
            p = rng.uniform(0.05, 0.95)
            r = np.array([rng.uniform(-5, 5) for _ in range(n)])
            u = rng.uniform(-3, 3)
            # Independence factorizes the exact Bernoulli log-mgf.
            exact = sum(
                math.log((1 - p) * math.exp(u * p * ri * (-p))
                         + p * math.exp(u * p * ri * (1 - p)))
                for ri in r
            )
            bound = u * u * p * p * float(np.dot(r, r)) / 8
            assert exact <= bound + 1e-11, (exact, bound)


def audit_a9_global_overlap() -> None:
    spec = importlib.util.spec_from_file_location(
        "selector_codebook_r22", ROOT / "tmp" / "selector_codebook_r22.py"
    )
    sc = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(sc)

    A = sc.PC.A9
    n = len(A)
    parent_q = sc.PC.matrix_norm(A)[0]
    child_q = [
        sc.PC.matrix_norm(sc.PC.principal(A, (i,))[0])[0]
        for i in range(n)
    ]
    grounds = []
    for orientation in (-1, 1):
        for tail in itertools.product((-1, 1), repeat=n - 1):
            x = (1,) + tail
            if sc.score(A, orientation, x) != parent_q:
                continue
            rows = tuple(
                orientation * x[i]
                * sum(A[i][j] * x[j] for j in range(n) if j != i)
                for i in range(n)
            )
            r2 = sum(v * v for v in rows)
            losses = tuple(
                2 * (child_q[i] - sc.score(A, orientation, x, (i,)))
                for i in range(n)
            )
            grounds.append((r2, losses))
    assert len(grounds) == 25

    for tolerance in (0, 4, 8):
        good = [r2 for r2, losses in grounds for loss in losses if loss <= tolerance]
        z_num = len(good)
        z_den = len(grounds) * n
        conditional_r2 = Fraction(sum(good), len(good))
        print(
            f"A9 ledger_t={tolerance}: global_overlap={z_num}/{z_den}, "
            f"conditional_R2={conditional_r2} ({float(conditional_r2):.12g})"
        )


if __name__ == "__main__":
    check_reference_identity()
    check_hoeffding()
    audit_a9_global_overlap()
    print("PASS: output-reference identity and domain-free linear Hoeffding bound")
