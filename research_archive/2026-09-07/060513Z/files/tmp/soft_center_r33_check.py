"""Finite checks for the Wave-33 soft-center identities."""

from itertools import product

import numpy as np
from scipy.optimize import minimize


def check_instance(seed: int) -> None:
    rng = np.random.default_rng(seed)
    nz, ns, s = 6, 4, 3
    lam, D = 0.73, 4
    costs = rng.integers(0, 6, size=(nz, ns))
    K = np.exp(-lam * costs)

    def moment(w: np.ndarray) -> float:
        return float(np.mean((K @ w) ** s))

    result = minimize(
        moment,
        np.full(ns, 1.0 / ns),
        method="SLSQP",
        bounds=[(0.0, 1.0)] * ns,
        constraints={"type": "eq", "fun": lambda w: np.sum(w) - 1.0},
        options={"ftol": 1e-13, "maxiter": 3000},
    )
    assert result.success, result.message
    w = result.x
    L = moment(w)
    V = L ** (1.0 / s)

    # The explicit norming functional is dual feasible and satisfies KKT.
    u = K @ w
    h = u ** (s - 1) / V ** (s - 1)
    q = s / (s - 1)
    assert abs(np.mean(h**q) - 1.0) < 2e-7
    overlaps = np.mean(h[:, None] * K, axis=0)
    assert np.min(overlaps) >= V - 2e-7
    assert np.max(np.abs(overlaps[w > 2e-5] - V)) < 2e-6

    # Enumerate the iid selector tuples and verify the exact transform and
    # both sides of the hard/soft sandwich.
    P = 0.0
    L_enum = 0.0
    for z in range(nz):
        for tup in product(range(ns), repeat=s):
            prob = np.prod(w[list(tup)]) / nz
            X = int(sum(costs[z, j] for j in tup))
            P += prob * (X <= D)
            L_enum += prob * np.exp(-lam * X)
    assert abs(L_enum - L) < 2e-12
    b = np.exp(-lam * (D + 1))
    lower = max(0.0, (L - b) / (1.0 - b))
    upper = np.exp(lam * D) * L
    assert P >= lower - 2e-12
    assert P <= upper + 2e-12


def check_shared_core() -> None:
    # This indicator matrix is the hard-kernel / lambda -> infinity toy used
    # only to demonstrate the strict normalized-degree versus L^s separation.
    nz, ns, s, core = 100, 7, 4, 2
    K = np.zeros((nz, ns))
    K[:core, :] = 1.0
    delta = core / nz
    w = np.full(ns, 1.0 / ns)
    V = np.mean((K @ w) ** s) ** (1.0 / s)
    h = np.zeros(nz)
    h[:core] = delta ** (-(s - 1) / s)
    q = s / (s - 1)
    overlaps = np.mean(h[:, None] * K, axis=0)
    assert abs(V - delta ** (1.0 / s)) < 1e-13
    assert abs(np.mean(h**q) - 1.0) < 1e-13
    assert np.max(np.abs(overlaps - V)) < 1e-13
    assert delta < V  # strict gain over normalized degrees for s>1.


def check_sharp_positive_branch() -> None:
    lam, D = 0.4, 3
    b = np.exp(-lam * (D + 1))
    for mass_at_zero in (0.0, 0.1, 0.7, 1.0):
        L = mass_at_zero + (1.0 - mass_at_zero) * b
        lower = max(0.0, (L - b) / (1.0 - b))
        assert abs(lower - mass_at_zero) < 1e-13


for trial_seed in range(10):
    check_instance(trial_seed)
check_shared_core()
check_sharp_positive_branch()
print("soft-center hard/soft, KKT dual, and shared-core checks passed")
